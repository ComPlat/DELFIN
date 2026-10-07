#requires -Version 5.1
[CmdletBinding(SupportsShouldProcess=$true)]
param([switch]$KeepProfiles)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$destination = Join-Path $env:LOCALAPPDATA 'DELFIN Launcher'
$shortcut = Join-Path ([Environment]::GetFolderPath('Desktop')) 'DELFIN.lnk'
$registry = 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\DELFINLauncher'
Write-Host 'Close the app and DELFIN terminal windows first. The server dashboard will not be stopped.'
if (-not $PSCmdlet.ShouldProcess($destination,'Uninstall the local DELFIN app')) { return }
if (Test-Path -LiteralPath $destination) {
    if ((Get-Item -LiteralPath $destination).Attributes -band [IO.FileAttributes]::ReparsePoint) { throw 'Installation directory is a link. Refusing recursive removal.' }
    if ($KeepProfiles) {
        foreach ($name in @('Install.cmd','Uninstall.cmd','DELFIN.ps1','Install.ps1','Uninstall.ps1','Test-Launcher.ps1','remote_launcher.py','remote_bootstrap.sh','README.md','DELFIN_logo.png','DELFIN.ico','Starter.cs','Build-Starter.ps1','DELFIN.exe')) {
            $path = Join-Path $destination $name
            if (Test-Path -LiteralPath $path) { Remove-Item -LiteralPath $path -Force }
        }
    } else { Remove-Item -LiteralPath $destination -Recurse -Force }
}
$desktopFolders = @([Environment]::GetFolderPath('Desktop'), [Environment]::GetFolderPath('DesktopDirectory'), (Join-Path $env:USERPROFILE 'Desktop'))
foreach ($cloudRoot in @($env:OneDrive, $env:OneDriveConsumer, $env:OneDriveCommercial)) {
    if ($cloudRoot) { $desktopFolders += Join-Path $cloudRoot 'Desktop' }
}
foreach ($desktop in @($desktopFolders | Where-Object { $_ } | Select-Object -Unique)) {
    $launcherShortcut = Join-Path $desktop 'DELFIN.lnk'
    if (Test-Path -LiteralPath $launcherShortcut) { Remove-Item -LiteralPath $launcherShortcut -Force }
}
$programs = [Environment]::GetFolderPath('Programs')
if ($programs) {
    $startShortcut = Join-Path $programs 'DELFIN.lnk'
    if (Test-Path -LiteralPath $startShortcut) { Remove-Item -LiteralPath $startShortcut -Force }
}
Write-Host 'If you pinned DELFIN to the taskbar, right-click its icon and select Unpin from taskbar.'
if (Test-Path -LiteralPath $registry) { Remove-Item -LiteralPath $registry -Recurse -Force }
Write-Host 'If the desktop icon is still visible, refresh the desktop with F5.'
Write-Host 'DELFIN Launcher uninstalled. The server installation and tmux session remain unchanged.'
if ($KeepProfiles) { Write-Host "Saved connections remain in $destination." }
