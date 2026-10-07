#requires -Version 5.1
[CmdletBinding()]
param([switch]$SkipOpenSSHCheck)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
foreach ($name in @('DELFIN.ps1','Install.ps1','Uninstall.ps1','Test-Launcher.ps1','Build-Starter.ps1')) {
    $tokens = $null
    $parseErrors = $null
    [void][Management.Automation.Language.Parser]::ParseFile((Join-Path $PSScriptRoot $name),[ref]$tokens,[ref]$parseErrors)
    if ($parseErrors.Count -gt 0) { throw ($parseErrors | Out-String) }
}
foreach ($name in @('Install.cmd','Uninstall.cmd','remote_launcher.py','remote_bootstrap.sh','README.md','DELFIN_logo.png','DELFIN.ico','Starter.cs','Build-Starter.ps1')) {
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

}

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
Write-Host 'Syntax, connections, SSH options, browser destination and icon checked. Test a real OTP connection separately.'

# Compile the windowless starter and check its Windows GUI subsystem and entry point.
$starterTestFolder = Join-Path ([IO.Path]::GetTempPath()) ('delfin-starter-test-' + [guid]::NewGuid().ToString())
[void][IO.Directory]::CreateDirectory($starterTestFolder)
try {
    $starterTest = Join-Path $starterTestFolder 'DELFIN.exe'
    & (Join-Path $PSScriptRoot 'Build-Starter.ps1') -Destination $starterTest
    $pe = [IO.File]::ReadAllBytes($starterTest)
    $peOffset = [BitConverter]::ToInt32($pe,0x3c)
    Assert-That ([BitConverter]::ToUInt16($pe,$peOffset+24+68) -eq 2) 'Starter must use the Windows GUI subsystem.'
    $selfTest = Start-Process -FilePath $starterTest -ArgumentList '--self-test' -Wait -PassThru
    Assert-That ($selfTest.ExitCode -eq 0) 'Starter self-test failed.'
    # Exercise the actual app-start path with a harmless script instead of the GUI.
    $probeScript = @'
Add-Type -TypeDefinition 'using System; using System.Runtime.InteropServices; public static class ConsoleProbe { [DllImport("kernel32.dll")] public static extern IntPtr GetConsoleWindow(); }'
[IO.File]::WriteAllText((Join-Path $PSScriptRoot 'console.txt'),[ConsoleProbe]::GetConsoleWindow().ToInt64().ToString())
'@
    $probeScript | Set-Content -LiteralPath (Join-Path $starterTestFolder 'DELFIN.ps1') -Encoding UTF8
    $appProbe = Start-Process -FilePath $starterTest -Wait -PassThru
    Assert-That ($appProbe.ExitCode -eq 0) 'Starter could not run its app script.'
    Assert-That ((Get-Content -LiteralPath (Join-Path $starterTestFolder 'console.txt') -Raw) -eq '0') 'App script must run without an allocated console.'
    # Launch the actual DELFIN GUI through its real starter, then inspect its HWND.
    foreach ($name in @('DELFIN.ps1','DELFIN_logo.png','DELFIN.ico')) {
        Copy-Item -LiteralPath (Join-Path $PSScriptRoot $name) -Destination (Join-Path $starterTestFolder $name) -Force
    }
    $guiProbe = Start-Process -FilePath $starterTest -ArgumentList '--gui-smoke-test' -PassThru
    try {
        if (-not $guiProbe.WaitForExit(30000)) { throw 'Actual DELFIN GUI did not complete its smoke test.' }
        Assert-That ($guiProbe.ExitCode -eq 0) 'Actual DELFIN GUI startup failed.'
        $guiReport = Get-Content -LiteralPath (Join-Path $starterTestFolder 'gui-smoke.json') -Raw | ConvertFrom-Json
        Assert-That ($guiReport.Visible -and $guiReport.ShowInTaskbar -and $guiReport.HasIcon) 'DELFIN GUI must be visible with its icon and taskbar entry enabled.'
        Assert-That ($guiReport.Owner -eq 0 -and ($guiReport.ExtendedStyle -band 0x80) -eq 0) 'DELFIN GUI must be an independent window, not an owned tool window.'
    } finally { if (-not $guiProbe.HasExited) { $guiProbe.Kill() }; $guiProbe.Dispose() }
} finally { Remove-Item -LiteralPath $starterTestFolder -Recurse -Force }
