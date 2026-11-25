#requires -Version 5.1
[CmdletBinding()]
param([switch]$SkipOpenSSHCheck)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
foreach ($name in @('DELFIN.ps1','Install.ps1','Uninstall.ps1','Test-Launcher.ps1','Setup-Wsl.ps1')) {
    $tokens = $null
    $parseErrors = $null
    [void][Management.Automation.Language.Parser]::ParseFile((Join-Path $PSScriptRoot $name),[ref]$tokens,[ref]$parseErrors)
    if ($parseErrors.Count -gt 0) { throw ($parseErrors | Out-String) }
}
foreach ($name in @('Install.cmd','Uninstall.cmd','remote_launcher.py','remote_bootstrap.sh','wsl_setup.sh','README.md','DELFIN_logo.png','DELFIN.ico')) {
    if (-not (Test-Path -LiteralPath (Join-Path $PSScriptRoot $name))) { throw "Missing file: $name" }
}
if (-not $SkipOpenSSHCheck -and -not (Test-Path "$env:SystemRoot\System32\OpenSSH\ssh.exe")) { throw 'Windows OpenSSH Client is missing.' }
# Load only pure validation functions: no GUI, network, subprocess, or profile writes.
. (Join-Path $PSScriptRoot 'DELFIN.ps1') -Mode Validate
function Assert-That([bool]$condition, [string]$message) { if (-not $condition) { throw $message } }
function Assert-Rejected([scriptblock]$action) {
    $rejected = $false
    try { & $action | Out-Null } catch { $rejected = $true }
    Assert-That $rejected 'Unsafe input was accepted.'
}
$legacy = Normalize-Profile ([pscustomobject]@{ Id=[guid]::NewGuid().ToString(); Target='alice@login.example.org' })
Assert-That (-not $legacy.KeepSession -and -not $legacy.OpenWorkingTerminal) 'Legacy connections must get explicit session defaults.'
Assert-That ($legacy.Directory -eq '' -and $legacy.Executable -eq '' -and $legacy.Port -eq 0) 'Older connections must receive optional defaults.'
Check-Profile $legacy
# Verify PowerShell 5.1 JSON arrays are enumerated into profiles, not normalized as one object.
$originalProfilesFile = $profilesFile
$roundTripFile = [IO.Path]::GetTempFileName()
try {
    $profilesFile = $roundTripFile
    ConvertTo-Json -InputObject @($legacy) | Set-Content -LiteralPath $profilesFile -Encoding UTF8
    $firstLoad = @(Read-Profiles)
    $secondLoad = @(Read-Profiles)
    Assert-That ($firstLoad.Count -eq 1 -and $firstLoad[0].Id -eq $legacy.Id -and $secondLoad[0].Id -eq $legacy.Id) 'Saved connection IDs must survive JSON array loading.'
    Assert-That ($firstLoad[0].Target -eq $legacy.Target) 'Saved SSH destination must survive JSON loading.'
} finally {
    $profilesFile = $originalProfilesFile
    Remove-Item -LiteralPath $roundTripFile -Force
}
$testProfile = [pscustomobject]@{ Id=[guid]::NewGuid().ToString(); Name='Cluster'; Target='alice@login.example.org'; Directory='/home/alice/project space'; Executable='/home/alice/env/bin/delfin-voila'; Port=8866 }
Check-Profile $testProfile
foreach ($bad in @('-oProxyCommand=bad','alice@host;bad','host with space',"host`ncommand")) {
    $testProfile.Target=$bad
    Assert-Rejected { Check-Profile $testProfile }
}
$testProfile.Target='cluster'
foreach ($bad in @('/home/a;command',"/home/a'command",'C:\Users\alice','/home/a|command')) {
    $testProfile.Directory=$bad
    Assert-Rejected { Check-Profile $testProfile }
}
foreach ($valid in @('relative/path','../repo','~/repo','/home/alice/projects','"/home/alice/projects"')) { $testProfile.Directory=$valid; Check-Profile $testProfile }
$testProfile.Directory='/home/alice/projects'
$testProfile.Port=80
Assert-Rejected { Check-Profile $testProfile }
$testProfile.Port=0
$testProfile.Directory=''
$testProfile.Executable=''
Check-Profile $testProfile
Assert-That ((Shell-Quote '/home/alice/project space') -eq "'/home/alice/project space'") 'Shell path was not quoted correctly.'
foreach ($kind in @('Dashboard','Terminal')) {
    $options=@(Get-SSHOptions $kind 8866)
    foreach ($required in @('StrictHostKeyChecking=ask','ForwardAgent=no','ForwardX11=no','PermitLocalCommand=no','ControlMaster=no','ControlPath=none')) {
        Assert-That ($options -contains $required) "Missing SSH protection: $required"
    }
    for ($i=0; $i -lt $options.Count; $i++) {
        if ($options[$i] -eq '-o') { Assert-That ($i+1 -lt $options.Count -and $options[$i+1].Contains('=')) 'Missing SSH option value.' }
    }
    Assert-That ($options -contains 'ClearAllForwardings=yes' -and $options -notcontains '-L') 'Connections without a tunnel must not inherit forwarding settings.'
    Assert-That ($options[0] -eq '-tt') 'SSH sessions need a terminal.'
}
foreach ($kind in @('Copy')) {
    $options=@(Get-SSHOptions $kind)
    foreach ($required in @('StrictHostKeyChecking=ask','ForwardAgent=no','ForwardX11=no','ClearAllForwardings=yes','PermitLocalCommand=no','ControlMaster=no','ControlPath=none')) {
        Assert-That ($options -contains $required) "Missing SSH protection for ${kind}: $required"
    }
    Assert-That ($options -notcontains '-tt' -and $options -notcontains '-L') "$kind must not ask for a terminal or a tunnel."
}
Assert-That ((@(Get-SSHOptions 'Copy')) -notcontains '-T') 'scp reads -T differently; the copy must not carry it.'
# The server half travels as two plain files into a content-addressed folder.
$version = Get-ScriptVersion
Assert-That ($version.Folder -match '^\.cache/delfin-launcher/[a-f0-9]{32}$') 'Script folder must be named by the scripts digests.'
$expectedBootstrap = (Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $PSScriptRoot 'remote_bootstrap.sh')).Hash.ToLowerInvariant()
Assert-That ($version.Bootstrap -eq $expectedBootstrap -and $version.Folder.EndsWith($expectedBootstrap.Substring(0,16) + $version.Launcher.Substring(0,16))) 'Script version must come from both files.'
$check = Get-ScriptCheckCommand $version
Assert-That ($check.StartsWith("bash -c 'umask 077") -and $check.EndsWith(' delfin ' + $version.Folder + ' ' + $version.Bootstrap + ' ' + $version.Launcher)) 'Script check must pass the folder and digests as arguments.'
Assert-That (-not $check.Contains('"')) 'Script check must not contain double quotes.'
$commandProfile = Normalize-Profile $testProfile
$remote = Get-RemoteCommand $commandProfile 'dashboard' '' '' $version
Assert-That ($remote.StartsWith($check + " && bash -ilc '. ~/" + $version.Folder + "/remote_bootstrap.sh' delfin 0 ")) 'One connection: the script check, then the login shell on the copied bootstrap.'
Assert-That ($check.Contains('echo DELFIN_SCRIPTS_MISSING; exit 10')) 'A missing copy must be reported by marker and exit code 10.'
Assert-That ($remote.Contains("'~/" + $version.Folder + "/remote_launcher.py' 'dashboard'")) 'The bootstrap must get the path of the copied launcher.'
Assert-That (-not $remote.Contains('"') -and -not $remote.Contains('base64') -and -not $remote.Contains('eval')) 'The remote command must carry paths and data only.'
Assert-Rejected { Get-RemoteCommand $commandProfile 'dashboard' '' '' ([pscustomobject]@{ Folder='.cache/delfin-launcher/../../x' }) }

$tunnelOptions = @(Get-SSHOptions 'Terminal' 8870 9000)
Assert-That ($tunnelOptions -contains '127.0.0.1:8870:127.0.0.1:9000' -and $tunnelOptions -contains 'ExitOnForwardFailure=yes') 'Local and remote tunnel ports must be independent.'
$occupied = [Net.Sockets.TcpListener]::new([Net.IPAddress]::Loopback,0)
try {
    $occupied.Start()
    $usedPort=([Net.IPEndPoint]$occupied.LocalEndpoint).Port
    $freePort=Get-FreeLocalPort $usedPort
    Assert-That ($freePort -ne $usedPort -and $freePort -ge 1024) 'Local port discovery must skip occupied ports.'
} finally { $occupied.Stop() }

$socketOptions = @(Get-SSHOptions 'Dashboard' 8870 0 '/tmp/delfin-win-tunnel-00000000-0000-0000-0000-000000000000/dashboard.sock')
Assert-That ($socketOptions -contains '127.0.0.1:8870:/tmp/delfin-win-tunnel-00000000-0000-0000-0000-000000000000/dashboard.sock') 'The dashboard must carry its own private tunnel without a second login.'
function Encode-Url([string]$url) { return [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($url)) }
$valid='http://127.0.0.1:8866/voila/render/dashboard.ipynb?token=' + ('a' * 43)
Assert-That ((Get-BrowserUrl (Encode-Url $valid) 8866) -eq $valid) 'Local token URL must be accepted.'
foreach ($bad in @($valid.Replace('127.0.0.1','evil.example'), $valid.Replace('8866','9000'), $valid.Replace('http:','file:'), $valid.Replace('/voila/render/','/api/'), $valid.Replace('127.0.0.1','user@127.0.0.1'), ($valid+'#fragment'), ($valid+'&extra=bad'))) {
    Assert-Rejected { Get-BrowserUrl (Encode-Url $bad) 8866 }
}
Add-Type -AssemblyName System.Drawing
$icon=New-Object Drawing.Icon((Join-Path $PSScriptRoot 'DELFIN.ico'))
$icon.Dispose()
Assert-That ((Quote-NativeArgument 'a "b" c') -eq '"a \"b\" c"') 'Native Windows arguments must escape quotation marks.'

# Local mode: DELFIN in a WSL distribution on this PC. Pure checks only; no
# WSL is needed here (Install.cmd runs this on PCs without it).
Assert-That ((Get-Transport $legacy) -eq 'SSH') 'A connection saved before WSL must stay an SSH connection.'
Assert-That ((Normalize-Profile ([pscustomobject]@{ Target='x' })).Transport -eq 'SSH') 'New connections default to SSH.'
$wslProfile = Normalize-Profile ([pscustomobject]@{ Id=[guid]::NewGuid().ToString(); Transport='WSL'; Target='Ubuntu-24.04' })
Check-Profile $wslProfile
foreach ($good in @('Ubuntu','DELFIN','Debian','openSUSE-Tumbleweed','Ubuntu_22.04')) { $wslProfile.Target=$good; Check-Profile $wslProfile }
foreach ($bad in @('','-d','Ubuntu;bad','Ubuntu bad','alice@Ubuntu','--exec',"Ubuntu`nbad",'Ubuntu"')) {
    $wslProfile.Target=$bad
    Assert-Rejected { Check-Profile $wslProfile }
}
$wslProfile.Target='Ubuntu'
foreach ($bad in @('wsl','ssh','Local','')) {
    $odd = [pscustomobject]@{ Id=$wslProfile.Id; Transport=$bad; Target='Ubuntu'; Directory=''; Executable=''; Port=0 }
    if ($bad) { Assert-Rejected { Check-Profile $odd } } else { Check-Profile $odd }
}
# wsl.exe reads its raw command line: plain options, then the command verbatim.
Assert-That ((Get-WslArguments 'Ubuntu' "sh -c 'a b' x") -ceq "-d Ubuntu --cd ~ -- sh -c 'a b' x") 'WSL starts the command in the Linux home, verbatim after --.'
Assert-Rejected { Get-WslArguments '-e' 'x' }
Assert-Rejected { Get-WslArguments 'Ubuntu' "x`nrm -rf ~" }
$localId = 'delfin-win-local-' + [guid]::NewGuid().ToString()
$wslRemote = Get-RemoteCommand $commandProfile 'dashboard' '' '' $version $localId
Assert-That ($wslRemote -eq ($remote + ' ' + (Shell-Quote $localId))) 'A WSL start is the SSH command plus its local connection as the last argument.'
foreach ($bad in @('delfin-win-local-x',"delfin-win-local-$([guid]::Empty)'; rm -rf ~",'delfin-win-tunnel-' + [guid]::Empty)) {
    Assert-Rejected { Get-RemoteCommand $commandProfile 'dashboard' '' '' $version $bad }
}
# wsl.exe -l -q: UTF-16 (NULs once read as UTF-8), a BOM, CRLF, Docker's own distributions.
$listing = [string][char]0xFEFF + "Ubuntu`r`ndocker-desktop`r`nDELFIN`r`n`r`nbad name`r`n"
$listing = -join ($listing.ToCharArray() | ForEach-Object { [string]$_ + [char]0 })
$names = @(ConvertFrom-WslListing $listing)
Assert-That (($names -join '|') -eq 'Ubuntu|DELFIN') "WSL distributions must be read from wsl.exe's listing: $($names -join '|')"
# The browser's token goes only to a listener the launcher trusts.
Assert-That (Test-TrustedListener 'SSH' @(4242) 4242) 'SSH: the connection''s own ssh.exe holds the tunnel.'
Assert-That (-not (Test-TrustedListener 'SSH' @(4243) 4242)) 'SSH: another program on the tunnel port is refused.'
Assert-That (-not (Test-TrustedListener 'SSH' @() 4242)) 'SSH: no tunnel yet is not ready.'
Assert-That (Test-TrustedListener 'WSL' @() 0) 'WSL mirrored networking: no Windows listener is the Linux side.'
Assert-That (-not (Test-TrustedListener 'WSL' @($PID) 0)) 'WSL: a Windows program other than the WSL relay is refused.'
$relay = @(Get-Process -Name wslrelay -ErrorAction SilentlyContinue | Where-Object { $_.Path })
if ($relay.Count) { Assert-That (Test-TrustedListener 'WSL' @($relay[0].Id) 0) 'WSL: the WSL relay from its install folder is trusted.' }
# The browser job gets the check as a function of its own (another runspace).
$trustCheck = [scriptblock]::Create('function Test-TrustedListener {' + ${function:Test-TrustedListener}.ToString() + '}')
$probe = Start-Job -InitializationScript $trustCheck -ScriptBlock { @((Test-TrustedListener 'SSH' @(7) 7), (Test-TrustedListener 'WSL' @($PID) 0)) }
try {
    [void](Wait-Job $probe -Timeout 60)
    $answers = @(Receive-Job $probe)
    Assert-That ($answers.Count -eq 2 -and $answers[0] -eq $true -and $answers[1] -eq $false) 'The browser job must run the same listener check.'
} finally { Remove-Job $probe -Force }
Write-Host 'Local mode (WSL): connection type, distribution names, command, listing and listener trust checked.'

# Setup of DELFIN on this PC (Setup-Wsl.ps1). Pure checks: nothing is
# installed, no registry entry is written, wsl.exe is not started.
. (Join-Path $PSScriptRoot 'Setup-Wsl.ps1') -FunctionsOnly
Assert-That ((Get-VirtualizationVerdict $false $false) -eq 'off') 'Setup: firmware virtualization off must stop the setup.'
Assert-That ((Get-VirtualizationVerdict $true $false) -eq 'on') 'Setup: a running hypervisor means virtualization is on.'
Assert-That ((Get-VirtualizationVerdict $false $null) -eq 'on') 'Setup: an unknown answer must not stop the setup.'
Assert-That ((Test-WindowsSupported 19041) -and -not (Test-WindowsSupported 18363)) 'Setup: WSL 2 needs build 19041.'
foreach ($case in @(@('Fehler 0x80370102','BIOS/UEFI'), @('Error code: Wsl/0x8007019e','Restart Windows'), @('0x80072ee7','internet'), @('E_ACCESSDENIED','administrator rights'))) {
    $text = Get-FriendlyWslError $case[0] 1
    Assert-That ($text.Contains($case[1]) -and -not $text.Contains('0x')) ('Setup: ' + $case[0] + ' must be explained in words.')
}
Assert-That ((Get-FriendlyWslError "line one`nsomething unknown" 5) -eq 'WSL reported an error (exit 5): something unknown') 'Setup: an unknown error shows its last line.'
$listed = ConvertFrom-WslList ([string][char]0xFEFF + "U`0b`0u`0n`0t`0u`0`r`0`n`0D`0E`0L`0F`0I`0N`0`r`0`n`0")
Assert-That ($listed.Count -eq 2 -and $listed[1] -eq 'DELFIN') 'Setup: wsl -l output must be read in every encoding.'
Assert-That ((Get-SetupArguments 'system') -ceq '-d DELFIN -u root --cd / -- bash -s -- system delfin') 'Setup: system step arguments.'
Assert-That ((Get-SetupArguments 'delfin') -ceq ('-d DELFIN -u delfin --cd ~ -- bash -s -- delfin https://github.com/ComPlat/DELFIN.git ' + $script:Commit) -and $script:Commit -match '^[0-9a-f]{40}$') 'Setup: DELFIN step arguments.'
Assert-Rejected { Get-SetupArguments 'other' }
Assert-Rejected { Invoke-Wsl ("-l`n-q") }
Assert-That (-not ((Get-SetupScriptBytes) -contains 13)) 'Setup: the Linux half must reach bash with LF line ends.'
Assert-That ((Get-ResumeCommand 'C:\x y\Setup-Wsl.ps1') -ceq ('"' + $env:SystemRoot + '\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy RemoteSigned -File "C:\x y\Setup-Wsl.ps1" -Resume')) 'Setup: resume command.'
Assert-Rejected { Get-ResumeCommand 'C:\x"y\Setup-Wsl.ps1' }
$setupStore = Join-Path ([IO.Path]::GetTempPath()) ('delfin-setup-test-' + [guid]::NewGuid().ToString())
try {
    [void][IO.Directory]::CreateDirectory($setupStore)
    $server = [pscustomobject]@{ Id=[guid]::NewGuid().ToString(); Name='cluster'; Target='alice@login.example.org'; Directory=''; Executable=''; Port=0; KeepSession=$true; OpenWorkingTerminal=$false; RememberLogin=$true }
    ConvertTo-Json -InputObject @($server) | Set-Content -LiteralPath (Join-Path $setupStore 'profiles.json') -Encoding UTF8
    $first = Add-WslProfile $setupStore
    $second = Add-WslProfile $setupStore
    # Windows PowerShell 5.1 emits a JSON array as one object: unwrap it.
    $savedJson = Get-Content -LiteralPath (Join-Path $setupStore 'profiles.json') -Raw | ConvertFrom-Json
    $saved = @($savedJson | ForEach-Object { $_ })
    Assert-That ($first -eq $second -and $saved.Count -eq 2) 'Setup: the WSL connection is added once.'
    Assert-That ($saved[0].Id -eq $server.Id -and $saved[0].KeepSession) 'Setup: existing connections stay as they are.'
    $local = Normalize-Profile $saved[1]
    Check-Profile $local
    Assert-That ((Get-Transport $local) -eq 'WSL' -and $local.Target -eq 'DELFIN' -and $local.Name -eq 'This PC (WSL DELFIN)' -and -not $local.Executable) 'Setup: the saved connection is the app''s WSL connection.'
} finally { Remove-Item -LiteralPath $setupStore -Recurse -Force }
Write-Host 'Setup on this PC (WSL): checks, error texts, arguments and saved connection checked.'
Write-Host 'Syntax, connections, SSH options, browser destination and icon checked. Test a real OTP connection separately.'

# Start the actual DELFIN window the way its shortcut does and inspect it.
$guiTestFolder = Join-Path ([IO.Path]::GetTempPath()) ('delfin-gui-test-' + [guid]::NewGuid().ToString())
[void][IO.Directory]::CreateDirectory($guiTestFolder)
try {
    foreach ($name in @('DELFIN.ps1','DELFIN_logo.png','DELFIN.ico')) {
        Copy-Item -LiteralPath (Join-Path $PSScriptRoot $name) -Destination (Join-Path $guiTestFolder $name) -Force
    }
    $powershell = Join-Path $env:SystemRoot 'System32\WindowsPowerShell\v1.0\powershell.exe'
    $arguments = '-NoProfile -ExecutionPolicy RemoteSigned -File ' + (Quote-NativeArgument (Join-Path $guiTestFolder 'DELFIN.ps1')) + ' -GuiSmokeTest'
    $guiProbe = Start-Process -FilePath $powershell -ArgumentList $arguments -PassThru
    try {
        if (-not $guiProbe.WaitForExit(30000)) { throw 'Actual DELFIN GUI did not complete its smoke test.' }
        Assert-That ($guiProbe.ExitCode -eq 0) 'Actual DELFIN GUI startup failed.'
        $guiReport = Get-Content -LiteralPath (Join-Path $guiTestFolder 'gui-smoke.json') -Raw | ConvertFrom-Json
        Assert-That ($guiReport.Visible -and $guiReport.ShowInTaskbar -and $guiReport.HasIcon) 'DELFIN GUI must be visible with its icon and taskbar entry enabled.'
        Assert-That ((-not $guiReport.HasOwner) -and $guiReport.TopLevel -and $guiReport.Title -eq 'DELFIN - SSH Dashboard') 'DELFIN GUI must be an independent top-level window.'
    } finally { if (-not $guiProbe.HasExited) { $guiProbe.Kill() }; $guiProbe.Dispose() }
} finally { Remove-Item -LiteralPath $guiTestFolder -Recurse -Force }
