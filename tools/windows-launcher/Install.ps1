#requires -Version 5.1
# Installs (or updates) the DELFIN app for the current Windows user. No
# administrator rights. The parameters exist for the isolated install test
# (Test-InstallRoundTrip.ps1); a normal install passes none of them.
[CmdletBinding()]
param(
    [string]$Destination = (Join-Path $env:LOCALAPPDATA 'DELFIN Launcher'),
    [string]$Desktop = [Environment]::GetFolderPath('Desktop'),
    [string]$Programs = [Environment]::GetFolderPath('Programs'),
    [string]$RegistryKey = 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\DELFINLauncher',
    [string]$SshFolder = (Join-Path $env:USERPROFILE '.ssh'),
    [switch]$CloseRunning
)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$destination = $Destination
if ($destination.Contains('"')) { throw 'Invalid installation path.' }
if (-not $Programs) { throw 'Windows Start Menu Programs folder is unavailable.' }

function Get-RunningApp([string]$folder) {
    # The app's own PowerShell hosts (main window and SSH consoles) run
    # DELFIN.ps1 from the installed folder; earlier versions ran DELFIN.exe.
    $script = (Join-Path $folder 'DELFIN.ps1').ToLowerInvariant()
    $exe = (Join-Path $folder 'DELFIN.exe').ToLowerInvariant()
    @(Get-CimInstance Win32_Process | Where-Object {
        ($_.CommandLine -and $_.CommandLine.ToLowerInvariant().Contains($script)) -or
        ($_.ExecutablePath -and $_.ExecutablePath.ToLowerInvariant() -eq $exe)
    })
}
function Stop-RunningApp($processes) {
    # The app and the ssh/scp it started: closing the app ends its connections
    # just as its own Disconnect does. A kept server session stays on the server.
    $ids = @($processes | ForEach-Object { $_.ProcessId })
    $children = @(Get-CimInstance Win32_Process | Where-Object { $ids -contains $_.ParentProcessId -and $_.Name -match '^(ssh|scp|powershell)\.exe$' })
    foreach ($process in @($children) + @($processes)) {
        try { Stop-Process -Id $process.ProcessId -Force -ErrorAction Stop } catch {}
    }
    Start-Sleep -Milliseconds 500
}

# An update must not replace the files of an app that is running: Windows
# refuses to remove the old DELFIN.exe while it runs, and a half-copied app
# breaks the connection that is open. Ask, and stop before anything changed.
$running = @(Get-RunningApp $destination)
if ($running.Count -gt 0) {
    Write-Host ''
    Write-Host 'DELFIN is open. It has to be closed for the update.' -ForegroundColor Yellow
    Write-Host 'Open connections end; a session kept on the server stays there.'
    $answer = if ($CloseRunning) { 'y' } else { Read-Host 'Close DELFIN now and continue? (y/N)' }
    if ($answer -ne 'y' -and $answer -ne 'Y') { throw 'DELFIN is still open. Close it, then run Install.cmd again.' }
    Stop-RunningApp $running
    if (@(Get-RunningApp $destination).Count -gt 0) { throw 'DELFIN could not be closed. Close its windows, then run Install.cmd again.' }
}

[void][IO.Directory]::CreateDirectory($destination)
foreach ($file in @('Install.cmd','Uninstall.cmd','DELFIN.ps1','Install.ps1','Uninstall.ps1','remote_launcher.py','remote_bootstrap.sh','Test-Launcher.ps1','Setup-Wsl.ps1','wsl_setup.sh','README.md','DELFIN_logo.png','DELFIN.ico')) {
    $source = Join-Path $PSScriptRoot $file
    $target = Join-Path $destination $file
    if ([IO.Path]::GetFullPath($source) -ne [IO.Path]::GetFullPath($target)) { Copy-Item -LiteralPath $source -Destination $target -Force }
}
# Files of earlier versions, which compiled a starter program at install time.
foreach ($old in @('DELFIN.exe','Starter.cs','Build-Starter.ps1','Taskbar.cs')) {
    $path = Join-Path $destination $old
    if (Test-Path -LiteralPath $path) { Remove-Item -LiteralPath $path -Force }
}
# Leftovers of connections that ended without cleaning up (a crash, a killed
# window): they only ever mattered while that connection ran.
Get-ChildItem -LiteralPath $destination -Filter 'connection-*' -File -ErrorAction SilentlyContinue | Remove-Item -Force -ErrorAction SilentlyContinue
# A standard Windows shortcut (WScript.Shell) starts Windows PowerShell on the
# installed DELFIN.ps1.
#
# -WindowStyle Hidden, and only here: the app is a window of its own, and the
# PowerShell console hosting it carried nothing a person could read or use.
# Started minimized it was still an empty black entry in the taskbar that
# people took for a second, broken application. A .lnk cannot express
# "hidden" (its WindowStyle is normal/minimized/maximized only), so the
# switch goes on the command line, and WindowStyle 7 stays as what Windows
# falls back to if the switch is ever dropped.
#
# The package otherwise forbids this switch, and the test enforcing that now
# makes this one exception by name: a hidden PowerShell window is a shape
# reputation scanners weigh against a file, and it is the shape malware
# uses. It is accepted for the app's own host process and nowhere else -- not
# for the SSH connection, which stays visible because that is where OpenSSH
# asks for a password.
$app = Join-Path $destination 'DELFIN.ps1'
$shell = New-Object -ComObject WScript.Shell
[void][IO.Directory]::CreateDirectory($Desktop)
$shortcut = $shell.CreateShortcut((Join-Path $Desktop 'DELFIN.lnk'))
$shortcut.TargetPath = "$env:SystemRoot\System32\WindowsPowerShell\v1.0\powershell.exe"
$shortcut.Arguments = '-NoProfile -ExecutionPolicy RemoteSigned -WindowStyle Hidden -File "' + $app + '"'
$shortcut.WorkingDirectory = $destination
$shortcut.WindowStyle = 7
$shortcut.IconLocation = (Join-Path $destination 'DELFIN.ico') + ',0'
$shortcut.Description = 'DELFIN SSH dashboard and terminal'
$shortcut.Save()
[void][IO.Directory]::CreateDirectory($Programs)
Copy-Item -LiteralPath (Join-Path $Desktop 'DELFIN.lnk') -Destination (Join-Path $Programs 'DELFIN.lnk') -Force

# Windows Settings > Apps > Installed apps. Uninstall there opens the same
# Uninstall.cmd as a double-click: one question, a result that stays on
# screen, and the choice to keep the saved connections.
[void](New-Item -Path $RegistryKey -Force)
$uninstall = '"' + "$env:SystemRoot\System32\cmd.exe" + '" /c "' + (Join-Path $destination 'Uninstall.cmd') + '"'
$sizeKb = [int](((Get-ChildItem -LiteralPath $destination -File | Measure-Object -Property Length -Sum).Sum) / 1KB)
$registration = @{ DisplayName='DELFIN SSH Launcher'; DisplayVersion='0.7.0'; Publisher='DELFIN'; InstallLocation=$destination; DisplayIcon=(Join-Path $destination 'DELFIN.ico'); UninstallString=$uninstall }
foreach ($entry in $registration.GetEnumerator()) {
    [void](New-ItemProperty -LiteralPath $RegistryKey -Name $entry.Key -Value $entry.Value -PropertyType String -Force)
}
foreach ($entry in @{ NoModify=1; NoRepair=1; EstimatedSize=$sizeKb }.GetEnumerator()) {
    [void](New-ItemProperty -LiteralPath $RegistryKey -Name $entry.Key -Value $entry.Value -PropertyType DWord -Force)
}

# An invitation (delfin-people --invite on the server) travels in the package
# as invite.json and invite_key. It becomes the first saved connection, the
# server's host keys go into known_hosts (so the first login asks no
# fingerprint question), and the invitation key is kept for the app's first
# Start, which registers this PC's own key with it and then deletes it.
$invited = $null
$inviteFile = Join-Path $PSScriptRoot 'invite.json'
$inviteKey = Join-Path $PSScriptRoot 'invite_key'
if ((Test-Path -LiteralPath $inviteFile) -and (Test-Path -LiteralPath $inviteKey)) {
    $invite = Get-Content -LiteralPath $inviteFile -Raw | ConvertFrom-Json
    $target = [string]$invite.target
    if ($target -notmatch '^[a-z_][a-z0-9_-]{0,31}@[A-Za-z0-9][A-Za-z0-9.-]{0,200}$') { throw 'The invitation names no valid server.' }
    $hostKeys = @($invite.host_keys | ForEach-Object { [string]$_ })
    if ($hostKeys.Count -eq 0 -or @($hostKeys | Where-Object { $_ -notmatch '^(ssh-ed25519|ecdsa-sha2-nistp(256|384|521)|ssh-rsa) [A-Za-z0-9+/]+={0,2}$' }).Count -gt 0) { throw 'The invitation carries no valid server key.' }
    $keyText = Get-Content -LiteralPath $inviteKey -Raw
    if (-not $keyText.StartsWith('-----BEGIN OPENSSH ' + 'PRIVATE KEY-----')) { throw 'The invitation key is not readable.' }
    $profilesFile = Join-Path $destination 'profiles.json'
    $profiles = @()
    if (Test-Path -LiteralPath $profilesFile) { $profiles = @(Get-Content -LiteralPath $profilesFile -Raw | ConvertFrom-Json | ForEach-Object { $_ } | Where-Object { $null -ne $_ }) }
    $existing = @($profiles | Where-Object { [string]$_.Target -eq $target })
    if ($existing.Count -gt 0) {
        $profileId = [string]$existing[0].Id
        $existing[0] | Add-Member -NotePropertyName RememberLogin -NotePropertyValue $true -Force
        $profiles = @($existing[0]) + @($profiles | Where-Object { [string]$_.Target -ne $target })
    } else {
        $profileId = [guid]::NewGuid().ToString()
        $profiles = @([pscustomobject]@{ Id=$profileId; Name=$target; Target=$target; Directory=''; Executable=''; Port=0; KeepSession=$true; OpenWorkingTerminal=$false; RememberLogin=$true }) + $profiles
    }
    ConvertTo-Json -InputObject @($profiles) -Depth 3 | Set-Content -LiteralPath $profilesFile -Encoding UTF8
    [void][IO.Directory]::CreateDirectory($SshFolder)
    $knownHosts = Join-Path $SshFolder 'known_hosts'
    $known = if (Test-Path -LiteralPath $knownHosts) { @(Get-Content -LiteralPath $knownHosts) } else { @() }
    $serverName = $target.Split('@')[1]
    foreach ($hostKey in $hostKeys) {
        $entry = $serverName + ' ' + $hostKey
        if ($known -notcontains $entry) { Add-Content -LiteralPath $knownHosts -Value $entry -Encoding ASCII }
    }
    $keptKey = Join-Path $destination ('invite-' + $profileId + '.key')
    if (Test-Path -LiteralPath $keptKey) { Remove-Item -LiteralPath $keptKey -Force }
    Copy-Item -LiteralPath $inviteKey -Destination $keptKey
    # Windows OpenSSH ignores a private key that anyone but its owner (and
    # SYSTEM / Administrators) may read -- "UNPROTECTED PRIVATE KEY FILE" --
    # and a copied file inherits whatever its folder grants. Measured: an
    # inherited extra entry made the invitation unusable. Owner only.
    $me = [Security.Principal.WindowsIdentity]::GetCurrent().User
    $acl = New-Object Security.AccessControl.FileSecurity
    $acl.SetAccessRuleProtection($true, $false)
    $acl.AddAccessRule((New-Object Security.AccessControl.FileSystemAccessRule($me, 'FullControl', 'Allow')))
    # Not Set-Acl: in Windows PowerShell 5.1 it also writes the audit list,
    # which needs administrator rights. This writes the access list only.
    [IO.File]::SetAccessControl($keptKey, $acl)
    $invited = $target
}

Write-Host ''
Write-Host 'DELFIN is installed.' -ForegroundColor Green
if ($invited) { Write-Host "  Your connection to $invited is set up. Start DELFIN and press Start / Reconnect: no password is needed." -ForegroundColor Green }
Write-Host '  Start it: the DELFIN icon on your desktop, or DELFIN in the Start menu.'
Write-Host '  Taskbar:  right-click DELFIN in Start, then Pin to taskbar (under More on some Windows versions).'
Write-Host '  Remove:   Windows Settings > Apps > Installed apps > DELFIN SSH Launcher.'
Write-Host '  Update:   run Install.cmd of a newer version. Saved connections are kept.'
