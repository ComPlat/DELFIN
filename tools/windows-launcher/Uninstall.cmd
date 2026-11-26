@echo off
setlocal
rem Not "cd /d %~dp0": Windows does not remove a folder that is a process's
rem current directory, so an uninstall started from the installed folder
rem removed nothing. Every file here is named by its full path instead.
cd /d "%SystemRoot%"
echo DELFIN - uninstall
echo.
echo This removes the DELFIN app from this PC.
echo Nothing on the server is changed: your sessions and calculations stay there.
echo.
echo   Y = remove DELFIN and its saved connections
echo   K = remove DELFIN, keep the saved connections for a later install
echo   N = cancel
echo.
choice /C YKN /N /M "Your choice [Y/K/N]: "
if errorlevel 3 (
    echo Nothing was changed.
    pause
    exit /b 0
)
rem One test: in a .cmd file "set" resets ERRORLEVEL, so the answer is read first.
if errorlevel 2 (set "keep=-KeepProfiles") else (set "keep=")
rem One block, read whole before the uninstaller removes this installed batch
rem file. "(goto) 2>nul" leaves the batch before cmd looks for its next line:
rem without it, a removed Uninstall.cmd ends with "batch file cannot be found"
rem and exit code 1 although everything was removed.
(
    "%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy RemoteSigned -File "%~dp0Uninstall.ps1" %keep% %*
    if errorlevel 1 (
        echo.
        echo Uninstall did not finish. Read the message above.
        pause
        (goto) 2>nul & exit /b 1
    )
    echo.
    pause
    (goto) 2>nul & exit /b 0
)
