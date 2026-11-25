#requires -Version 5.1
# Sets up DELFIN on this PC for the launcher's local mode ("This PC (WSL)"):
#   1. WSL itself, if missing (Windows asks once for administrator approval;
#      a restart follows, after which this script continues once by itself
#      through a RunOnce entry -- no lasting autostart);
#   2. a WSL distribution of its own, "DELFIN" (Ubuntu 24.04), so the
#      person's other distributions stay untouched;
#   3. tmux, git and a compiler, and a Linux user, as root (wsl -u root needs
#      no password, sudo would ask for one);
#   4. DELFIN through its install.sh, as that user;
#   5. the saved connection "This PC (WSL DELFIN)" in the launcher.
# Every step can run again; what is already there is kept.
[CmdletBinding()]
param([switch]$Resume, [switch]$FunctionsOnly)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$script:Distribution = 'DELFIN'
$script:LinuxUser = 'delfin'
$script:Image = 'Ubuntu-24.04'
$script:Repository = 'https://github.com/ComPlat/DELFIN.git'
# The DELFIN commit this launcher installs (git verifies it: the id is the
# checksum of the whole tree). Bumped with each launcher release.
$script:Commit = '6903045ff2cb44530e39c1d67698c77940fa11a7'
$script:ResumeName = 'DELFINWslSetup'
$script:ResumeKey = 'HKCU:\Software\Microsoft\Windows\CurrentVersion\RunOnce'
# WSL 2 and "wsl --install" need Windows 10 version 2004 (build 19041) or newer.
$script:MinimumBuild = 19041

class SetupStop : Exception { SetupStop([string]$message) : base($message) {} }

function Test-WindowsSupported([int]$build) { return $build -ge $script:MinimumBuild }

function Get-VirtualizationVerdict($hypervisorPresent, $firmwareEnabled) {
    # Win32_Processor reports firmware virtualization as off while a
    # hypervisor already runs (Hyper-V, WSL itself), so a running hypervisor
    # wins. $null means Windows could not tell: not a reason to stop.
    if ($hypervisorPresent -eq $true) { return 'on' }
    if ($firmwareEnabled -eq $false) { return 'off' }
    return 'on'
}

function Get-FriendlyWslError([string]$output, [int]$code) {
    # What wsl.exe prints is in the Windows language and full of codes; the
    # codes are the same everywhere, so they are what is matched.
    if ($output -match '0x80370102|HCS_E_HYPERV_NOT_INSTALLED|0x80370114|0x8037010d') {
        return 'Virtualization is switched off in this PC''s firmware (BIOS/UEFI). Switch on "Intel VT-x" / "Intel Virtualization Technology" or "AMD-V" / "SVM Mode" there, start Windows again, and run Install.cmd again. On a work PC, ask your IT.'
    }
    if ($output -match '0x8007019e|0x800701bc|WSL_E_WSL_OPTIONAL_COMPONENT_REQUIRED|0x80370106') {
        return 'WSL is installed but not active yet. Restart Windows, then run Install.cmd again.'
    }
    if ($output -match '0x80072ee7|0x80072efd|0x80072ee2|0x80072f8f|WININET_E') {
        return 'The download failed: no connection to the internet (or it is blocked). Check the connection and run Install.cmd again.'
    }
    if ($output -match '0x80070005|E_ACCESSDENIED|0x800704c7') {
        return 'Windows refused the change: administrator rights are needed. On a work PC, ask your IT to install WSL, then run Install.cmd again.'
    }
    if ($output -match 'WSL_E_DISTRO_NOT_FOUND|0x8007000d.*--name|--name') {
        return 'This WSL version is too old for a distribution of its own. Run "wsl --update" once, then run Install.cmd again.'
    }
    $last = @($output -split "`r?`n" | Where-Object { $_.Trim() }) | Select-Object -Last 1
    return 'WSL reported an error (exit ' + $code + '): ' + $last
}

function ConvertFrom-WslList([string]$text) {
    # wsl -l -q: UTF-16 unless WSL_UTF8=1, a BOM at times, NUL padding.
    return @(($text -replace "`0", '' -replace [char]0xFEFF, '') -split "`r?`n" | ForEach-Object { $_.Trim() } | Where-Object { $_ })
}

function Get-WslProgramPath { return (Join-Path $env:SystemRoot 'System32\wsl.exe') }

function Test-WslInstalled {
    # The WSL from the Store (what "wsl --install" puts there) lives under
    # Program Files\WSL. System32\wsl.exe alone is only Windows' stub, which
    # offers to install WSL.
    return (Test-Path -LiteralPath (Join-Path $env:ProgramFiles 'WSL\wsl.exe')) -and (Test-Path -LiteralPath (Get-WslProgramPath))
}

function Invoke-Wsl([string]$arguments, [byte[]]$stdin = $null, [switch]$Quiet) {
    # wsl.exe reads its raw command line: arguments go in unquoted, built
    # from checked values only (no user text reaches this string).
    if ($arguments -match "[`r`n`0]") { throw 'Invalid WSL arguments.' }
    $info = New-Object Diagnostics.ProcessStartInfo
    $info.FileName = Get-WslProgramPath
    $info.Arguments = $arguments
    $info.UseShellExecute = $false
    $info.RedirectStandardOutput = $true
    $info.RedirectStandardError = $true
    $info.RedirectStandardInput = $true
    $info.StandardOutputEncoding = [Text.Encoding]::UTF8
    $info.StandardErrorEncoding = [Text.Encoding]::UTF8
    $info.EnvironmentVariables['WSL_UTF8'] = '1'
    $process = [Diagnostics.Process]::Start($info)
    $errors = $process.StandardError.ReadToEndAsync()
    if ($stdin) { $process.StandardInput.BaseStream.Write($stdin, 0, $stdin.Length) }
    $process.StandardInput.Close()
    $lines = New-Object Collections.Generic.List[string]
    while ($null -ne ($line = $process.StandardOutput.ReadLine())) {
        $clean = $line -replace "`0", ''
        $lines.Add($clean)
        if (-not $Quiet -and $clean.Trim()) { Write-Host ('   ' + $clean) }
    }
    $process.WaitForExit()
    $all = (($lines -join "`n") + "`n" + $errors.Result) -replace "`0", ''
    return [pscustomobject]@{ Code = $process.ExitCode; Output = $all }
}

function Get-SetupScriptBytes {
    $path = Join-Path $PSScriptRoot 'wsl_setup.sh'
    $bytes = [IO.File]::ReadAllBytes($path)
    # bash must see LF; a checkout or editor may have turned it into CRLF.
    $text = [Text.Encoding]::UTF8.GetString($bytes) -replace "`r`n", "`n"
    return [Text.UTF8Encoding]::new($false).GetBytes($text)
}

function Get-SetupArguments([string]$step) {
    if ($step -eq 'system') { return '-d ' + $script:Distribution + ' -u root --cd / -- bash -s -- system ' + $script:LinuxUser }
    if ($step -eq 'delfin') { return '-d ' + $script:Distribution + ' -u ' + $script:LinuxUser + ' --cd ~ -- bash -s -- delfin ' + $script:Repository + ' ' + $script:Commit }
    throw 'Unknown setup step.'
}

function Get-ResumeCommand([string]$scriptPath) {
    if ($scriptPath.Contains('"')) { throw 'Invalid installation path.' }
    $powershell = Join-Path $env:SystemRoot 'System32\WindowsPowerShell\v1.0\powershell.exe'
    return '"' + $powershell + '" -NoProfile -ExecutionPolicy RemoteSigned -File "' + $scriptPath + '" -Resume'
}

function Register-Resume([string]$scriptPath) {
    # RunOnce: Windows deletes the entry before running it, once, at the next
    # logon of this user. Nothing stays in the autostart.
    [void](New-Item -Path $script:ResumeKey -Force -ErrorAction SilentlyContinue)
    [void](New-ItemProperty -LiteralPath $script:ResumeKey -Name $script:ResumeName -Value (Get-ResumeCommand $scriptPath) -PropertyType String -Force)
}

function Unregister-Resume {
    if (Get-ItemProperty -LiteralPath $script:ResumeKey -Name $script:ResumeName -ErrorAction SilentlyContinue) {
        Remove-ItemProperty -LiteralPath $script:ResumeKey -Name $script:ResumeName -Force
    }
}

function Add-WslProfile([string]$store) {
    # The launcher's saved connections: kept as they are, one WSL connection
    # to the DELFIN distribution added unless there is one already.
    $file = Join-Path $store 'profiles.json'
    $profiles = @()
    if (Test-Path -LiteralPath $file) {
        $loaded = Get-Content -LiteralPath $file -Raw | ConvertFrom-Json
        $profiles = @($loaded | Where-Object { $null -ne $_ })
    }
    foreach ($p in $profiles) {
        $transport = $p.PSObject.Properties['Transport']
        if ($null -ne $transport -and [string]$transport.Value -eq 'WSL' -and [string]$p.Target -eq $script:Distribution) { return $p.Id }
    }
    $id = [guid]::NewGuid().ToString()
    $profiles += [pscustomobject]@{ Id=$id; Name=('This PC (WSL ' + $script:Distribution + ')'); Transport='WSL'; Target=$script:Distribution; Directory=''; Executable=''; Port=0; KeepSession=$false; OpenWorkingTerminal=$false; RememberLogin=$false }
    [void][IO.Directory]::CreateDirectory($store)
    ConvertTo-Json -InputObject @($profiles) | Set-Content -LiteralPath $file -Encoding UTF8
    return $id
}

function Install-WslFeature {
    Write-Host ''
    Write-Host 'WSL (Windows Subsystem for Linux) is not installed yet. DELFIN runs in it.'
    Write-Host 'Windows now asks ONCE for administrator approval to install it.'
    Write-Host 'Afterwards Windows has to restart; setup then continues by itself.'
    Write-Host ''
    try {
        $process = Start-Process -FilePath (Get-WslProgramPath) -ArgumentList '--install --no-distribution' -Verb RunAs -Wait -PassThru
    } catch [ComponentModel.Win32Exception] {
        if ($_.Exception.NativeErrorCode -eq 1223) {
            throw [SetupStop]::new('Administrator approval was declined, so WSL was not installed. Run Install.cmd again and approve, or on a work PC ask your IT to install WSL.')
        }
        throw
    }
    if ($process.ExitCode -ne 0) {
        throw [SetupStop]::new((Get-FriendlyWslError ('exit 0x{0:x8}' -f $process.ExitCode) $process.ExitCode))
    }
}

function Invoke-Setup([string]$store) {
    $build = [Environment]::OSVersion.Version.Build
    if (-not (Test-WindowsSupported $build)) {
        throw [SetupStop]::new('This Windows (build ' + $build + ') is too old for WSL 2. It needs Windows 10 version 2004 or newer, or Windows 11. Use DELFIN on a server (SSH) instead.')
    }
    if (-not (Test-WslInstalled)) {
        $system = Get-CimInstance Win32_ComputerSystem
        $processor = @(Get-CimInstance Win32_Processor)[0]
        $firmware = $processor.PSObject.Properties['VirtualizationFirmwareEnabled']
        if ((Get-VirtualizationVerdict $system.HypervisorPresent $(if ($firmware) { $firmware.Value } else { $null })) -eq 'off') {
            throw [SetupStop]::new((Get-FriendlyWslError 'HCS_E_HYPERV_NOT_INSTALLED' 0))
        }
        Install-WslFeature
        Register-Resume (Join-Path $PSScriptRoot 'Setup-Wsl.ps1')
        return 'restart'
    }
    # With no distribution at all, "wsl -l" fails; the install below then
    # reports any real WSL problem in its own words.
    $listing = Invoke-Wsl '-l -q' -Quiet
    $names = if ($listing.Code -eq 0) { ConvertFrom-WslList $listing.Output } else { @() }
    if ($names -notcontains $script:Distribution) {
        Write-Host ('Downloading and setting up Linux (' + $script:Image + ') as "' + $script:Distribution + '". This takes a few minutes.')
        $created = Invoke-Wsl ('--install ' + $script:Image + ' --name ' + $script:Distribution + ' --no-launch')
        if ($created.Code -ne 0) { throw [SetupStop]::new((Get-FriendlyWslError $created.Output $created.Code)) }
    }
    Write-Host 'Installing tmux, git and a compiler in Linux...'
    $system = Invoke-Wsl (Get-SetupArguments 'system') (Get-SetupScriptBytes)
    if ($system.Code -ne 0) { throw [SetupStop]::new((Get-FriendlyWslError $system.Output $system.Code)) }
    # The default user from /etc/wsl.conf takes effect after a restart of the distribution.
    [void](Invoke-Wsl ('--terminate ' + $script:Distribution) -Quiet)
    Write-Host 'Installing DELFIN in Linux. This takes a while (about 10-20 minutes).'
    $delfin = Invoke-Wsl (Get-SetupArguments 'delfin') (Get-SetupScriptBytes)
    if ($delfin.Code -ne 0 -or $delfin.Output -notmatch 'DELFIN_WSL_READY') {
        throw [SetupStop]::new('Installing DELFIN in Linux failed (exit ' + $delfin.Code + '). The lines above say why. Run Install.cmd again to continue where it stopped.')
    }
    [void](Add-WslProfile $store)
    return 'ready'
}

if ($FunctionsOnly) { return }

Write-Host 'DELFIN on this PC (WSL) - setup'
if ($Resume) { Write-Host 'Continuing after the restart.' }
$exitCode = 0
try {
    $result = Invoke-Setup (Join-Path $env:LOCALAPPDATA 'DELFIN Launcher')
    Write-Host ''
    if ($result -eq 'restart') {
        Write-Host 'WSL is installed. Restart Windows now.' -ForegroundColor Yellow
        Write-Host 'After the restart and your next sign-in, setup continues by itself (once).'
        $exitCode = 3
    } else {
        Unregister-Resume
        Write-Host 'Done. Open DELFIN, choose "This PC (WSL DELFIN)" and click Start.' -ForegroundColor Green
    }
} catch [SetupStop] {
    Write-Host ''
    Write-Host $_.Exception.Message -ForegroundColor Red
    Write-Host 'The server mode (SSH) works without WSL.'
    $exitCode = 2
}
if ($Resume) { Read-Host 'Press Enter to close' | Out-Null }
exit $exitCode
