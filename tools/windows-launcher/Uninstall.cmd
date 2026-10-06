@echo off
setlocal
cd /d "%~dp0"
echo DELFIN Launcher - uninstall
echo Close DELFIN app and terminal windows first.
echo The server dashboard will not be stopped.
echo.
choice /C YN /N /M "Uninstall the local app and saved connections? [Y/N]: "
if errorlevel 2 exit /b 0
rem Read the whole block before the uninstaller may remove this installed batch file.
if errorlevel 1 (
    "%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy RemoteSigned -File "%~dp0Uninstall.ps1"
    if errorlevel 1 (
        echo.
        echo Uninstall failed. Read the error above and README.md.
        pause
        exit /b 1
    ) else (
        echo.
        echo Uninstall complete.
        pause
        exit /b 0
    )
)
exit /b 1
