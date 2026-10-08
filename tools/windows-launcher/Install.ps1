#requires -Version 5.1
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$destination = Join-Path $env:LOCALAPPDATA 'DELFIN Launcher'
[void][IO.Directory]::CreateDirectory($destination)
foreach ($file in @('Install.cmd','Uninstall.cmd','DELFIN.ps1','Install.ps1','Uninstall.ps1','remote_launcher.py','remote_bootstrap.sh','Test-Launcher.ps1','README.md','DELFIN_logo.png','DELFIN.ico','Starter.cs','Build-Starter.ps1','Taskbar.cs')) {
    $source = Join-Path $PSScriptRoot $file
    $target = Join-Path $destination $file
    if ([IO.Path]::GetFullPath($source) -ne [IO.Path]::GetFullPath($target)) { Copy-Item -LiteralPath $source -Destination $target -Force }
}
$starter = Join-Path $destination 'DELFIN.exe'
& (Join-Path $destination 'Build-Starter.ps1') -Destination $starter
$shell = New-Object -ComObject WScript.Shell
$shortcut = $shell.CreateShortcut((Join-Path ([Environment]::GetFolderPath('Desktop')) 'DELFIN.lnk'))
$shortcut.TargetPath = $starter
$shortcut.Arguments = ''
$shortcut.WorkingDirectory = $destination
$shortcut.IconLocation = (Join-Path $destination 'DELFIN.ico') + ',0'
$shortcut.Description = 'DELFIN SSH dashboard and terminal'
$shortcut.Save()
$programs = [Environment]::GetFolderPath('Programs')
if (-not $programs) { throw 'Windows Start Menu Programs folder is unavailable.' }
[void][IO.Directory]::CreateDirectory($programs)
Copy-Item -LiteralPath (Join-Path ([Environment]::GetFolderPath('Desktop')) 'DELFIN.lnk') -Destination (Join-Path $programs 'DELFIN.lnk') -Force
Write-Host "Installed: $destination"
Write-Host 'DELFIN desktop and Start menu shortcuts created.'
Write-Host 'Taskbar: find DELFIN in Start, right-click, then select Pin to taskbar (under More on some Windows versions).'
Write-Host 'Read README for script execution requirements.'

$registry = 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\DELFINLauncher'
[void](New-Item -Path $registry -Force)
$uninstall = '"' + "$env:SystemRoot\System32\WindowsPowerShell\v1.0\powershell.exe" + '" -NoProfile -ExecutionPolicy RemoteSigned -File "' + (Join-Path $destination 'Uninstall.ps1') + '"'
$registration = @{ DisplayName='DELFIN SSH Launcher'; DisplayVersion='0.5.0'; Publisher='DELFIN'; InstallLocation=$destination; DisplayIcon=(Join-Path $destination 'DELFIN.ico'); UninstallString=$uninstall }
foreach ($entry in $registration.GetEnumerator()) {
    [void](New-ItemProperty -LiteralPath $registry -Name $entry.Key -Value $entry.Value -PropertyType String -Force)
}
Write-Host 'Uninstall: Windows Installed apps or Uninstall.ps1.'
