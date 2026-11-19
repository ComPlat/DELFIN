@echo off
setlocal
cd /d "%~dp0"
echo DELFIN Launcher - installation
echo.
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
echo Installation failed. Nothing was installed. Read the error above,
echo then START-HERE.txt for the usual causes, or README.md for the detail.
pause
exit /b 1
