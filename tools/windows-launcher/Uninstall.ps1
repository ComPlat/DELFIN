#requires -Version 5.1
# Removes the DELFIN app of the current Windows user. Nothing on the server is
# touched. The path parameters exist for the isolated install test
# (Test-InstallRoundTrip.ps1); a normal uninstall passes none of them.
[CmdletBinding(SupportsShouldProcess=$true)]
param(
    [switch]$KeepProfiles,
    [switch]$CloseRunning,
    # Also delete the WSL distribution "DELFIN" that Setup-Wsl.ps1 made,
    # with everything in it. Only when asked for.
    [switch]$RemoveWslDistribution,
    [string]$Destination = (Join-Path $env:LOCALAPPDATA 'DELFIN Launcher'),
    [string]$Desktop = '',
    [string]$Programs = [Environment]::GetFolderPath('Programs'),
    [string]$RegistryKey = 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\DELFINLauncher'
)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$destination = $Destination
$registry = $RegistryKey
# Windows does not remove a folder that is any process's current directory,
# and Remove-Item then removes nothing at all. Uninstall.cmd may have been
# started inside the installed folder.
Set-Location -LiteralPath $env:SystemRoot
[Environment]::CurrentDirectory = $env:SystemRoot

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
    # The app and the ssh/scp it started, as its own Disconnect ends them.
    $ids = @($processes | ForEach-Object { $_.ProcessId })
    $children = @(Get-CimInstance Win32_Process | Where-Object { $ids -contains $_.ParentProcessId -and $_.Name -match '^(ssh|scp|powershell)\.exe$' })
    foreach ($process in @($children) + @($processes)) {
        try { Stop-Process -Id $process.ProcessId -Force -ErrorAction Stop } catch {}
    }
    Start-Sleep -Milliseconds 500
}

if (-not $PSCmdlet.ShouldProcess($destination,'Uninstall the local DELFIN app')) { return }
# A running app holds its files: removing them under it fails halfway.
$running = @(Get-RunningApp $destination)
if ($running.Count -gt 0) {
    Write-Host 'DELFIN is open. It has to be closed to uninstall it.' -ForegroundColor Yellow
    Write-Host 'Open connections end; a session kept on the server stays there.'
    $answer = if ($CloseRunning) { 'y' } else { Read-Host 'Close DELFIN now and continue? (y/N)' }
    if ($answer -ne 'y' -and $answer -ne 'Y') { throw 'DELFIN is still open. Nothing was removed.' }
    Stop-RunningApp $running
    if (@(Get-RunningApp $destination).Count -gt 0) { throw 'DELFIN could not be closed. Close its windows and try again. Nothing was removed.' }
}
if (Test-Path -LiteralPath $destination) {
    if ((Get-Item -LiteralPath $destination).Attributes -band [IO.FileAttributes]::ReparsePoint) { throw 'Installation directory is a link. Refusing recursive removal.' }
    if ($KeepProfiles) {
        foreach ($name in @('Install.cmd','Uninstall.cmd','DELFIN.ps1','Install.ps1','Uninstall.ps1','Test-Launcher.ps1','Setup-Wsl.ps1','wsl_setup.sh','remote_launcher.py','remote_bootstrap.sh','README.md','DELFIN_logo.png','DELFIN.ico','DELFIN.exe','Starter.cs','Build-Starter.ps1','Taskbar.cs')) {
            $path = Join-Path $destination $name
            if (Test-Path -LiteralPath $path) { Remove-Item -LiteralPath $path -Force }
        }
        Get-ChildItem -LiteralPath $destination -Filter 'connection-*' -File | Remove-Item -Force
    } else { Remove-Item -LiteralPath $destination -Recurse -Force }
}
if ($Desktop) {
    $desktopFolders = @($Desktop)
} else {
    $desktopFolders = @([Environment]::GetFolderPath('Desktop'), [Environment]::GetFolderPath('DesktopDirectory'), (Join-Path $env:USERPROFILE 'Desktop'))
    foreach ($cloudRoot in @($env:OneDrive, $env:OneDriveConsumer, $env:OneDriveCommercial)) {
        if ($cloudRoot) { $desktopFolders += Join-Path $cloudRoot 'Desktop' }
    }
}
foreach ($folder in @($desktopFolders | Where-Object { $_ } | Select-Object -Unique)) {
    $launcherShortcut = Join-Path $folder 'DELFIN.lnk'
    if (Test-Path -LiteralPath $launcherShortcut) { Remove-Item -LiteralPath $launcherShortcut -Force }
}
if ($Programs) {
    $startShortcut = Join-Path $Programs 'DELFIN.lnk'
    if (Test-Path -LiteralPath $startShortcut) { Remove-Item -LiteralPath $startShortcut -Force }
}
if (Test-Path -LiteralPath $registry) { Remove-Item -LiteralPath $registry -Recurse -Force }
# A WSL setup still waiting for its restart (Setup-Wsl.ps1) would otherwise run once more.
$resume = 'HKCU:\Software\Microsoft\Windows\CurrentVersion\RunOnce'
if (Get-ItemProperty -LiteralPath $resume -Name 'DELFINWslSetup' -ErrorAction SilentlyContinue) { Remove-ItemProperty -LiteralPath $resume -Name 'DELFINWslSetup' -Force }
$wsl = Join-Path $env:SystemRoot 'System32\wsl.exe'
$linux = $false
if (Test-Path -LiteralPath $wsl) {
    $env:WSL_UTF8 = '1'
    $names = @((& $wsl -l -q 2>$null) -replace "`0", '' | ForEach-Object { $_.Trim() } | Where-Object { $_ })
    $linux = $names -contains 'DELFIN'
}
if ($linux -and -not $RemoveWslDistribution) {
    # Asked only where Setup-Wsl.ps1 made it; the answer defaults to keeping it.
    $answer = Read-Host "Also delete DELFIN's Linux part on this PC (WSL distribution DELFIN, with all files in it)? (y/N)"
    $RemoveWslDistribution = $answer -eq 'y' -or $answer -eq 'Y'
}
if ($RemoveWslDistribution -and $linux) {
    & $wsl --terminate DELFIN | Out-Null
    & $wsl --unregister DELFIN
    if ($LASTEXITCODE -ne 0) { throw 'Could not remove the WSL distribution DELFIN.' }
    $linux = $false
    Write-Host 'The WSL distribution DELFIN and everything in it were removed.'
}
Write-Host ''
Write-Host 'DELFIN is removed from this PC.' -ForegroundColor Green
Write-Host '  Nothing on the server was changed: your sessions and calculations there are untouched.'
if ($KeepProfiles) { Write-Host "  Saved connections are kept in $destination for a later install." }
if ($linux) { Write-Host '  DELFIN on this PC (WSL distribution DELFIN) is kept. "wsl --unregister DELFIN" deletes it.' }
Write-Host '  If you pinned DELFIN to the taskbar: right-click its icon, then Unpin from taskbar.'
Write-Host '  If the desktop icon is still visible, press F5 on the desktop.'
