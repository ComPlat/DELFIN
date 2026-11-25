#requires -Version 5.1
# Install, update, uninstall -- for real, but into a folder of its own: its own
# app folder, desktop, Start menu and registry key, so the person's installed
# DELFIN is not touched. Developer-side (CI: windows-launcher.yml); not part
# of the package.
[CmdletBinding()]
param()
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
function Assert-That([bool]$condition, [string]$message) { if (-not $condition) { throw "FAILED: $message" } }

$powershell = Join-Path $env:SystemRoot 'System32\WindowsPowerShell\v1.0\powershell.exe'
$root = Join-Path ([IO.Path]::GetTempPath()) ('delfin-roundtrip-' + [guid]::NewGuid().ToString('N'))
$dest = Join-Path $root 'app'
$desk = Join-Path $root 'desktop'
$prog = Join-Path $root 'programs'
$reg = 'HKCU:\Software\DELFINLauncherTest\' + [guid]::NewGuid().ToString('N')
$paths = @('-Destination',$dest,'-Desktop',$desk,'-Programs',$prog,'-RegistryKey',$reg)
$fake = $null
function Invoke-Script([string]$script, [string[]]$arguments) {
    & $powershell -NoProfile -NonInteractive -ExecutionPolicy RemoteSigned -File $script @arguments | Out-Host
    return $LASTEXITCODE
}
function Invoke-Uninstall([string]$answer) {
    # The installed Uninstall.cmd, started INSIDE the installed folder (as a
    # double-click there does), answering its question on standard input.
    $answerFile = Join-Path $root ('answer-' + $answer + '.txt')
    Set-Content -LiteralPath $answerFile -Value $answer -Encoding ASCII
    $arguments = '/c ""' + (Join-Path $dest 'Uninstall.cmd') + '" ' + (($paths | ForEach-Object { '"' + $_ + '"' }) -join ' ') + '"'
    $process = Start-Process -FilePath (Join-Path $env:SystemRoot 'System32\cmd.exe') -ArgumentList $arguments -WorkingDirectory $dest -RedirectStandardInput $answerFile -NoNewWindow -Wait -PassThru
    return $process.ExitCode
}
try {
    # 1. Fresh install.
    Assert-That ((Invoke-Script (Join-Path $PSScriptRoot 'Install.ps1') $paths) -eq 0) 'install'
    foreach ($name in @('DELFIN.ps1','Install.cmd','Uninstall.cmd','Uninstall.ps1','remote_bootstrap.sh','remote_launcher.py','DELFIN.ico')) {
        Assert-That (Test-Path -LiteralPath (Join-Path $dest $name)) "installed file $name"
    }
    $shell = New-Object -ComObject WScript.Shell
    $link = $shell.CreateShortcut((Join-Path $desk 'DELFIN.lnk'))
    Assert-That ($link.TargetPath -like '*\WindowsPowerShell\v1.0\powershell.exe' -and $link.Arguments.Contains((Join-Path $dest 'DELFIN.ps1'))) 'desktop shortcut starts the installed app'
    Assert-That (Test-Path -LiteralPath (Join-Path $prog 'DELFIN.lnk')) 'Start menu shortcut'
    $entry = Get-ItemProperty -LiteralPath $reg
    Assert-That ($entry.UninstallString.Contains((Join-Path $dest 'Uninstall.cmd')) -and $entry.UninstallString.Contains('cmd.exe')) 'Installed apps opens Uninstall.cmd'
    Assert-That ($entry.NoModify -eq 1 -and $entry.NoRepair -eq 1 -and $entry.EstimatedSize -gt 0) 'Installed apps entry offers Uninstall only and shows a size'

    # 2. Update over an install with saved connections and leftovers.
    Set-Content -LiteralPath (Join-Path $dest 'profiles.json') -Value '[]' -Encoding UTF8
    Set-Content -LiteralPath (Join-Path $dest 'connection-00000000-0000-0000-0000-000000000000.json') -Value '1' -Encoding ASCII
    Set-Content -LiteralPath (Join-Path $dest 'Starter.cs') -Value '// earlier version' -Encoding ASCII
    Assert-That ((Invoke-Script (Join-Path $PSScriptRoot 'Install.ps1') $paths) -eq 0) 'update'
    Assert-That (Test-Path -LiteralPath (Join-Path $dest 'profiles.json')) 'an update keeps the saved connections'
    Assert-That (-not (Test-Path -LiteralPath (Join-Path $dest 'Starter.cs'))) 'an update removes files of earlier versions'
    Assert-That (@(Get-ChildItem -LiteralPath $dest -Filter 'connection-*').Count -eq 0) 'an update removes leftovers of ended connections'

    # 3. The app is running: no answer means no change; -CloseRunning closes it.
    $fake = Start-Process -FilePath $powershell -ArgumentList @('-NoProfile','-Command',('Start-Sleep 120 # ' + (Join-Path $dest 'DELFIN.ps1'))) -WindowStyle Hidden -PassThru
    Start-Sleep -Seconds 1
    $before = (Get-Item -LiteralPath (Join-Path $dest 'DELFIN.ps1')).LastWriteTimeUtc
    Assert-That ((Invoke-Script (Join-Path $PSScriptRoot 'Install.ps1') $paths) -ne 0) 'an update does not go on while the app runs and nobody answered'
    Assert-That (-not $fake.HasExited) 'the running app was left alone'
    Assert-That ((Invoke-Script (Join-Path $PSScriptRoot 'Install.ps1') ($paths + '-CloseRunning')) -eq 0) 'update after closing the app'
    Assert-That ($fake.WaitForExit(5000)) 'the app was closed for the update'
    $fake = $null

    # 4. Uninstall from the installed folder, keeping the saved connections.
    Assert-That ((Invoke-Uninstall 'K') -eq 0) 'uninstall (keep connections) from inside the installed folder'
    Assert-That (Test-Path -LiteralPath (Join-Path $dest 'profiles.json')) 'K keeps the saved connections'
    Assert-That (-not (Test-Path -LiteralPath (Join-Path $dest 'DELFIN.ps1'))) 'K removes the app'
    Assert-That (-not (Test-Path -LiteralPath (Join-Path $desk 'DELFIN.lnk')) -and -not (Test-Path -LiteralPath (Join-Path $prog 'DELFIN.lnk'))) 'shortcuts removed'
    Assert-That (-not (Test-Path -LiteralPath $reg)) 'Installed apps entry removed'

    # 5. Install again (the connections are still there), then remove everything.
    Assert-That ((Invoke-Script (Join-Path $PSScriptRoot 'Install.ps1') $paths) -eq 0) 'reinstall'
    Assert-That (Test-Path -LiteralPath (Join-Path $dest 'profiles.json')) 'a reinstall finds the kept connections'
    Assert-That ((Invoke-Uninstall 'N') -eq 0 -and (Test-Path -LiteralPath (Join-Path $dest 'DELFIN.ps1'))) 'N changes nothing'
    Assert-That ((Invoke-Uninstall 'Y') -eq 0) 'uninstall (everything) from inside the installed folder'
    Assert-That (-not (Test-Path -LiteralPath $dest)) 'Y removes the app folder'
    # 6. An invitation package (delfin-people --invite): the connection, the
    #    pinned server key and the invitation key are set up by the install.
    $package = Join-Path $root 'invited package'
    [void](New-Item -ItemType Directory -Path $package)
    Get-ChildItem -LiteralPath $PSScriptRoot -File | Copy-Item -Destination $package
    $hostKey = 'ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIGxVZr1q3F0c0o1u2y3m4n5p6q7r8s9t0u1v2w3x4y5z'
    $inviteJson = @{ target='delfin@192.0.2.10'; name='anna'; expires='209901010000Z'; host_keys=@($hostKey) }
    ConvertTo-Json -InputObject $inviteJson | Set-Content -LiteralPath (Join-Path $package 'invite.json') -Encoding ASCII
    # A real key, so OpenSSH's own permission check can be run on the copy.
    $keygen = Join-Path $env:SystemRoot 'System32\OpenSSH\ssh-keygen.exe'
    & $keygen -q -t ed25519 -N '""' -f (Join-Path $package 'invite_key') | Out-Null
    Remove-Item -LiteralPath (Join-Path $package 'invite_key.pub') -Force
    $ssh = Join-Path $root 'ssh'
    Assert-That ((Invoke-Script (Join-Path $package 'Install.ps1') ($paths + @('-SshFolder',$ssh))) -eq 0) 'install with an invitation'
    $saved = @(Get-Content -LiteralPath (Join-Path $dest 'profiles.json') -Raw | ConvertFrom-Json)
    Assert-That ($saved[0].Target -eq 'delfin@192.0.2.10' -and $saved[0].RememberLogin) 'the invitation is the first saved connection, with Remember login'
    Assert-That ((Get-Content -LiteralPath (Join-Path $ssh 'known_hosts')) -contains ('192.0.2.10 ' + $hostKey)) 'the server key is pinned in known_hosts'
    $kept = Join-Path $dest ('invite-' + $saved[0].Id + '.key')
    Assert-That (Test-Path -LiteralPath $kept) 'the invitation key is kept for the first Start'
    # ssh-keygen -y refuses a key file others may read ("UNPROTECTED PRIVATE
    # KEY FILE"), exactly as ssh does when it logs in with it.
    & $keygen -y -f $kept 2>&1 | Out-Null
    Assert-That ($LASTEXITCODE -eq 0) 'OpenSSH accepts the permissions of the kept invitation key'
    Assert-That ((Invoke-Script (Join-Path $package 'Install.ps1') ($paths + @('-SshFolder',$ssh))) -eq 0) 'installing the same invitation again'
    $again = @(Get-Content -LiteralPath (Join-Path $dest 'profiles.json') -Raw | ConvertFrom-Json)
    Assert-That (@($again | Where-Object { $_.Target -eq 'delfin@192.0.2.10' }).Count -eq 1) 'no second connection for the same server'
    Assert-That (@(Get-Content -LiteralPath (Join-Path $ssh 'known_hosts') | Where-Object { $_ -eq ('192.0.2.10 ' + $hostKey) }).Count -eq 1) 'no duplicate known_hosts line'
    Assert-That ((Invoke-Uninstall 'Y') -eq 0 -and -not (Test-Path -LiteralPath $dest)) 'uninstall after an invitation'
    Write-Host 'Install, update, running app, uninstall (keep / cancel / all), invitation checked in an isolated folder.'
} finally {
    if ($fake -and -not $fake.HasExited) { Stop-Process -Id $fake.Id -Force }
    if (Test-Path -LiteralPath $reg) { Remove-Item -LiteralPath $reg -Recurse -Force }
    $parent = 'HKCU:\Software\DELFINLauncherTest'
    if ((Test-Path -LiteralPath $parent) -and @(Get-ChildItem -LiteralPath $parent).Count -eq 0) { Remove-Item -LiteralPath $parent -Force }
    if (Test-Path -LiteralPath $root) { Remove-Item -LiteralPath $root -Recurse -Force -ErrorAction SilentlyContinue }
}
