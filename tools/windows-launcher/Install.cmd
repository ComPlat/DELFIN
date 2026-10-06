@echo off
setlocal
cd /d "%~dp0"
echo DELFIN Launcher - installation
"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy RemoteSigned -File "%~dp0Test-Launcher.ps1"
if errorlevel 1 goto failed
"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy RemoteSigned -File "%~dp0Install.ps1"
if errorlevel 1 goto failed
echo.
echo Installation complete. Open DELFIN from your desktop.
pause
exit /b 0
:failed
echo.
echo Installation failed. Read the error above and README.md.
pause
exit /b 1
