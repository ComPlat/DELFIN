@echo off
setlocal
cd /d "%~dp0"
echo DELFIN - installation
echo.
rem Double-clicking Install.cmd inside the ZIP runs it from a temporary copy
rem without the files beside it: say "extract first" instead of failing on
rem a missing file.
if not exist "%~dp0DELFIN.ps1" goto notextracted
if not exist "%~dp0remote_bootstrap.sh" goto notextracted
rem The app logs in with Windows' own SSH client, an optional Windows feature.
if not exist "%SystemRoot%\System32\OpenSSH\ssh.exe" goto nossh
rem Windows marks a file that came from the internet, and PowerShell's
rem RemoteSigned policy then refuses to run it unless it is signed. That is
rem the first thing most people hit, and what they used to see was the raw
rem security error in their own Windows language followed by "read the
rem error above". The check is inline rather than in a .ps1 of its own: a
rem script that did the checking would be refused for the same reason.
rem Unblocking is offered, never done silently -- the mark is a protection,
rem and removing it is the person's decision, not the installer's.
"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy RemoteSigned -Command "$b = @(Get-ChildItem -LiteralPath . -Filter *.ps1 | Where-Object { Get-Item -LiteralPath $_.FullName -Stream Zone.Identifier -ErrorAction SilentlyContinue }); if ($b.Count -eq 0) { exit 0 }; Write-Host 'Windows marked these files as downloaded, so PowerShell will not run them:' -ForegroundColor Yellow; $b | ForEach-Object { Write-Host ('   ' + $_.Name) }; Write-Host ''; Write-Host 'Unblocking clears that mark on these files and nothing else. No setting'; Write-Host 'and no policy is changed. Do it only if you trust where this came from.'; Write-Host ''; $a = Read-Host 'Unblock them and continue? (y/N)'; if ($a -eq 'y' -or $a -eq 'Y') { $b | Unblock-File; Write-Host 'Unblocked.'; exit 0 }; Write-Host 'Nothing was changed. You can also right-click each file, Properties, Unblock.'; exit 1"
if errorlevel 1 goto failed
echo Checking the package...
"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy RemoteSigned -File "%~dp0Test-Launcher.ps1" >nul
if errorlevel 1 goto failed
echo Package OK.
"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy RemoteSigned -File "%~dp0Install.ps1"
if errorlevel 1 goto failed
echo.
echo DELFIN can also run on this PC itself, in WSL (Linux inside Windows),
echo with no server. That needs about 7 GB and 10-40 minutes. If WSL is not
echo installed yet, Windows asks once for administrator approval and has to
echo restart; setup then continues by itself after you sign in again.
choice /C YN /N /M "Set up DELFIN on this PC too? [Y/N]: "
if errorlevel 2 goto startapp
"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy RemoteSigned -File "%LOCALAPPDATA%\DELFIN Launcher\Setup-Wsl.ps1"
rem 3: WSL was installed and Windows has to restart first; setup continues after it.
if errorlevel 3 (
    pause
    exit /b 0
)
if errorlevel 1 (
    echo.
    echo The DELFIN app is installed and works with a server. Only the setup
    echo on this PC stopped; the reason is above. Run Install.cmd again later.
)
:startapp
echo.
choice /C YN /N /M "Start DELFIN now? [Y/N]: "
if errorlevel 2 exit /b 0
if exist "%APPDATA%\Microsoft\Windows\Start Menu\Programs\DELFIN.lnk" (
    start "" "%APPDATA%\Microsoft\Windows\Start Menu\Programs\DELFIN.lnk"
) else (
    echo Open DELFIN from your desktop.
    pause
)
exit /b 0
:notextracted
echo The files of this package are not next to Install.cmd.
echo.
echo Extract the ZIP first: right-click the ZIP file, "Extract All...",
echo choose a new empty folder, then double-click Install.cmd in that folder.
echo Nothing was installed.
pause
exit /b 1
:nossh
echo Windows' SSH client (OpenSSH Client) is not installed on this PC.
echo DELFIN logs in to your server with it.
echo.
echo Add it in Windows Settings: System, Optional features, "Add a feature"
echo (Windows 10: Apps, Optional features), search for OpenSSH Client, install it.
echo Then double-click Install.cmd again. The settings page opens now.
echo Nothing was installed.
start "" ms-settings:optionalfeatures
pause
exit /b 1
:failed
echo.
echo Installation failed. Nothing was installed. Read the error above,
echo then START-HERE.txt for the usual causes, or README.md for the detail.
pause
exit /b 1
