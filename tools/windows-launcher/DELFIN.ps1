#requires -Version 5.1
[CmdletBinding()]
param([ValidateSet('App','Dashboard','Terminal','Validate')][string]$Mode = 'App', [string]$ProfileId, [int]$LocalPort = 0, [int]$RemotePort = 0, [string]$ConnectionId, [switch]$HideLauncherConsole, [switch]$GuiSmokeTest)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$store = Join-Path $env:LOCALAPPDATA 'DELFIN Launcher'
$profilesFile = Join-Path $store 'profiles.json'

function Normalize-Profile($value) {
    $defaults = [ordered]@{ Id=[guid]::NewGuid().ToString(); Name=''; Target=''; Directory=''; Executable=''; Port=0; KeepSession=$false; OpenWorkingTerminal=$false }
    foreach ($key in @($defaults.Keys)) {
        $property = $value.PSObject.Properties[$key]
        if ($null -ne $property -and $null -ne $property.Value) { $defaults[$key] = $property.Value }
    }
    if (-not $defaults.Name) { $defaults.Name = $defaults.Target }
    return [pscustomobject]$defaults
}
function Read-Profiles {
    if (Test-Path $profilesFile) {
        $loaded = Get-Content -LiteralPath $profilesFile -Raw | ConvertFrom-Json
        foreach ($value in $loaded) { if ($null -ne $value) { Normalize-Profile $value } }
    }
}
function Check-Profile($p) {
    if ($p.Id -notmatch '^[a-fA-F0-9]{8}-[a-fA-F0-9]{4}-[a-fA-F0-9]{4}-[a-fA-F0-9]{4}-[a-fA-F0-9]{12}$') { throw 'Invalid connection ID.' }
    if ($p.Target -notmatch '^[a-zA-Z0-9][a-zA-Z0-9_.@-]{0,200}$') { throw 'Enter an SSH alias or user@hostname (no options).' }
    foreach ($field in @('Directory','Executable')) {
        $path = ([string]$p.$field).Trim()
        if ($path.Length -ge 2 -and (($path.StartsWith('"') -and $path.EndsWith('"')) -or ($path.StartsWith("'") -and $path.EndsWith("'")))) { $path = $path.Substring(1,$path.Length-2).Trim() }
        $label = if ($field -eq 'Directory') { 'Working directory' } else { 'DELFIN location' }
        if ($path -and $path -notmatch '^(~(/.*)?|[\p{L}\p{N}_./ -]+)$') { throw "$label must be a Linux path, not a Windows path or command. Leave it empty for automatic discovery, or use /path, ~/path, or a relative path." }
        $p.$field = $path
    }
    if ([int]$p.Port -ne 0 -and ([int]$p.Port -lt 1024 -or [int]$p.Port -gt 65535)) { throw 'Starting port must be 0 (automatic) or between 1024 and 65535.' }
}
function Shell-Quote([string]$value) {
    # POSIX single-quote escaping, without double quotes in the native Windows argument.
    return "'" + $value.Replace("'", "'\''") + "'"
}
function Quote-NativeArgument([string]$value) {
    $escaped = [regex]::Replace($value,'(\\*)"','$1$1\"')
    $escaped = [regex]::Replace($escaped,'(\\+)$','$1$1')
    return '"' + $escaped + '"'
}
function Get-RemoteCommand($p, [string]$role, [string]$tunnelPath = '') {
    $bootstrap = [Convert]::ToBase64String([IO.File]::ReadAllBytes((Join-Path $PSScriptRoot 'remote_bootstrap.sh')))
    $payload = [Convert]::ToBase64String([IO.File]::ReadAllBytes((Join-Path $PSScriptRoot 'remote_launcher.py')))
    $script = 'eval "$(printf %s ' + $bootstrap + ' | base64 -d)"'
    return 'bash -ilc ' + (Shell-Quote $script) + ' delfin ' + [int]$p.Port + ' ' + (Shell-Quote $p.Directory) + ' ' + (Shell-Quote $p.Executable) + ' ' + (Shell-Quote $payload) + ' ' + (Shell-Quote $role) + ' ' + (Shell-Quote $tunnelPath) + ' ' + $(if ($p.KeepSession) { '1' } else { '0' }) + ' ' + $(if ($p.OpenWorkingTerminal) { '1' } else { '0' })
}
function Get-FreeLocalPort([int]$preferred) {
    $candidates = if ($preferred -eq 0) { @(0) } else { @($preferred..([Math]::Min(65535,$preferred+100))) + @(0) }
    foreach ($candidate in $candidates) {
        $listener = [Net.Sockets.TcpListener]::new([Net.IPAddress]::Loopback,[int]$candidate)
        try { $listener.Start(); return ([Net.IPEndPoint]$listener.LocalEndpoint).Port }
        catch { continue }
        finally { $listener.Stop() }
    }
    throw 'No local port available.'
}
function Get-SSHOptions([string]$kind, [int]$local = 0, [int]$remote = 0, [string]$tunnelPath = '') {
    $hasTunnel = $kind -eq 'Terminal' -and $local -ge 1024 -and $local -le 65535 -and $remote -ge 1024 -and $remote -le 65535
    $privateTunnel = $kind -eq 'Dashboard' -and $local -ge 1024 -and $local -le 65535 -and $tunnelPath -match '^/tmp/delfin-win-tunnel-[a-f0-9-]{36}/dashboard.sock$'
    $forwardingPolicy = if ($hasTunnel -or $privateTunnel) { 'ClearAllForwardings=no' } else { 'ClearAllForwardings=yes' }
    $options = @('-tt','-o','StrictHostKeyChecking=ask','-o','ForwardAgent=no','-o','ForwardX11=no','-o',$forwardingPolicy,'-o','PermitLocalCommand=no','-o','ServerAliveInterval=30','-o','ServerAliveCountMax=3','-o','ControlMaster=no','-o','ControlPath=none')
    if ($hasTunnel) {
        $options += @('-o','GatewayPorts=no','-o','ExitOnForwardFailure=yes','-L',"127.0.0.1:${local}:127.0.0.1:${remote}")
    }
    if ($privateTunnel) {
        $options += @('-o','GatewayPorts=no','-o','ExitOnForwardFailure=yes','-L',"127.0.0.1:${local}:$tunnelPath")
    }
    return $options
}
function Get-BrowserUrl([string]$encoded, [int]$port) {
    $url = [Text.Encoding]::UTF8.GetString([Convert]::FromBase64String($encoded))
    $uri = [Uri]$url
    if (-not $uri.IsAbsoluteUri -or $uri.Scheme -ne 'http' -or $uri.Host -ne '127.0.0.1' -or ($port -ne 0 -and $uri.Port -ne $port) -or $uri.Port -lt 1024 -or $uri.Port -gt 65535 -or $uri.UserInfo -or $uri.Fragment -or ($uri.AbsolutePath -ne '/' -and -not $uri.AbsolutePath.StartsWith('/voila/render/')) -or $uri.Query -notmatch '^\?token=[A-Za-z0-9_-]{32,128}(&session=[A-Za-z0-9%_.~+-]+)?$') {
        throw 'Invalid browser URL from the server.'
    }
    return $url
}
function Start-Window([string]$kind, [string]$id, [int]$local = 0, [int]$remote = 0, [string]$connection = '') {
    if ($PSCommandPath.Contains('"')) { throw 'Invalid installation path.' }
    if (-not $connection) { $connection = [guid]::NewGuid().ToString() }
    $starter = Join-Path $PSScriptRoot 'DELFIN.exe'
    if ($kind -ne 'Dashboard' -or -not (Test-Path -LiteralPath $starter)) { throw 'Reinstall DELFIN using Install.cmd to enable the browser-only starter.' }
    $worker = Start-Process -FilePath $starter -ArgumentList @('--dashboard-worker',$id,$connection) -PassThru
    [void]$script:ownedConnections.Add([pscustomobject]@{ Profile=$id; Connection=$connection; Worker=$worker })
}

if ($Mode -eq 'Validate') { return }

# Only the installed GUI shortcut opts in; never hide an interactive SSH worker
# or a user's existing PowerShell window when they run this script manually.
if ($Mode -eq 'App' -and $HideLauncherConsole) {
    Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
public static class DelfinLauncherWindow {
    [DllImport("kernel32.dll")] public static extern IntPtr GetConsoleWindow();
    [DllImport("user32.dll")] public static extern bool ShowWindow(IntPtr window, int command);
    [DllImport("user32.dll")] public static extern bool SetForegroundWindow(IntPtr window);
    [DllImport("user32.dll")] public static extern bool IsWindowVisible(IntPtr window);
}
'@
    $launcherConsole = [DelfinLauncherWindow]::GetConsoleWindow()
    if ($launcherConsole -ne [IntPtr]::Zero) {
        [void][DelfinLauncherWindow]::ShowWindow($launcherConsole, 0)
    }
}

if ($Mode -ne 'App') {
    $env:TERM = 'xterm-256color'
    # tmux sends ANSI terminal control sequences; enable Windows VT rendering.
    Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
public static class DelfinConsole {
    [DllImport("kernel32.dll")] public static extern IntPtr GetConsoleWindow();
    [DllImport("user32.dll")] public static extern bool ShowWindow(IntPtr window, int command);
    [DllImport("user32.dll")] public static extern IntPtr SendMessage(IntPtr window, uint message, IntPtr wParam, IntPtr lParam);
    [DllImport("kernel32.dll")] public static extern IntPtr GetStdHandle(int handle);
    [DllImport("kernel32.dll")] public static extern bool GetConsoleMode(IntPtr handle, out uint mode);
    [DllImport("kernel32.dll")] public static extern bool SetConsoleMode(IntPtr handle, uint mode);
}
'@
    $consoleHandle = [DelfinConsole]::GetStdHandle(-11)
    [uint32]$consoleMode = 0
    if ([DelfinConsole]::GetConsoleMode($consoleHandle,[ref]$consoleMode)) {
        [void][DelfinConsole]::SetConsoleMode($consoleHandle,($consoleMode -bor 4))
    }
    function Update-WorkerControls {
        if ($stopFile -and (Test-Path -LiteralPath $stopFile)) {
            $script:disconnectRequested = $true
            if ($null -ne $sshProcess -and -not $sshProcess.HasExited) { $sshProcess.Kill() }
        }
        if ($showFile -and (Test-Path -LiteralPath $showFile)) {
            [void][DelfinConsole]::ShowWindow($loginWindow,9)
            [void][DelfinConsole]::SetForegroundWindow($loginWindow)
            $result = if ($loginWindow -ne [IntPtr]::Zero -and [DelfinConsole]::IsWindowVisible($loginWindow)) { 'Terminal restored.' } else { 'Windows could not restore this terminal window.' }
            [IO.File]::WriteAllText((Join-Path $store ('connection-' + $connection + '.show-status')),$result)
            Remove-Item -LiteralPath $showFile -Force
        }
    }
    $sshProcess = $null
    $browserJob = $null
    $readyFile = $null
    $stopFile = $null
    $browserReadyFile = $null
    $showFile = $null
    $disconnectRequested = $false
    $remoteClosed = $false
    $loginWindow = [DelfinConsole]::GetConsoleWindow()
    Add-Type -AssemblyName System.Drawing
    $terminalIcon = [Drawing.Icon]::ExtractAssociatedIcon("$env:SystemRoot\System32\WindowsPowerShell\v1.0\powershell.exe")
    if ($loginWindow -ne [IntPtr]::Zero -and $terminalIcon) {
        [void][DelfinConsole]::SendMessage($loginWindow,0x80,[IntPtr]::Zero,$terminalIcon.Handle)
        [void][DelfinConsole]::SendMessage($loginWindow,0x80,[IntPtr]1,$terminalIcon.Handle)
    }
    try {
        $matchingProfiles = @(Read-Profiles | Where-Object { $_.Id -eq $ProfileId })
        if ($matchingProfiles.Count -ne 1) { throw 'Connection not found.' }
        $p = $matchingProfiles[0]
        Check-Profile $p
        $ssh = Join-Path $env:SystemRoot 'System32\OpenSSH\ssh.exe'
        if (-not (Test-Path $ssh)) { throw 'Windows OpenSSH Client is missing. See README.' }
        $Host.UI.RawUI.WindowTitle = "SSH $Mode - $($p.Name)"
        $tunnelPath = ''
        if ($Mode -eq 'Dashboard') {
            $LocalPort = Get-FreeLocalPort ([int]$p.Port)
            $tunnelPath = '/tmp/delfin-win-tunnel-' + [guid]::NewGuid().ToString() + '/dashboard.sock'
        }
        $options = @(Get-SSHOptions $Mode $LocalPort $RemotePort $tunnelPath)
        if ($ConnectionId) {
            if ($ConnectionId -notmatch '^[a-f0-9-]{36}$' -or ($Mode -eq 'Terminal' -and ($LocalPort -lt 1024 -or $RemotePort -lt 1024))) { throw 'Invalid tunnel request.' }
            $readyFile = Join-Path $store ('connection-' + $ConnectionId + '.json')
        }
        if ($Mode -eq 'Terminal') {
            Write-Host 'Working terminal. Password and OTP prompts are handled exclusively by OpenSSH.'
            $terminalCommand = Get-RemoteCommand $p 'terminal'
            $terminalInfo = New-Object Diagnostics.ProcessStartInfo
            $terminalInfo.FileName = $ssh
            $terminalArgs = @($options) + @([string]$p.Target,$terminalCommand)
            $terminalInfo.Arguments = ($terminalArgs | ForEach-Object { Quote-NativeArgument $_ }) -join ' '
            $terminalInfo.UseShellExecute = $false
            $sshProcess = [Diagnostics.Process]::Start($terminalInfo)
            if ($readyFile) { [string]$sshProcess.Id | Set-Content -LiteralPath $readyFile -Encoding ASCII }
            $sshProcess.WaitForExit()
            Write-Host "SSH ended (exit code $($sshProcess.ExitCode))."
        } else {
            $port = [int]$p.Port
            $remote = Get-RemoteCommand $p 'dashboard' $tunnelPath
            Write-Host 'Dashboard terminal: enter your password and OTP directly in OpenSSH.'
            Write-Host 'On first login, compare the host fingerprint with the official server fingerprint.'
            Write-Host 'The browser signs in automatically.'
            if ($p.KeepSession) {
                Write-Host 'Keep session ON: close this SSH window or detach with Ctrl+B, then D. Reconnect: start the app again.'
            } else {
                Write-Host 'Keep session OFF: closing this SSH window stops a dashboard newly started by this connection.'
                Write-Host 'A dashboard that was already running keeps its existing lifetime settings.'
            }
            Write-Host 'Ctrl+C in the dashboard explicitly stops it, even with Keep session ON.'
            # Password/OTP use the inherited console; stdout is pumped without logging.
            $info = New-Object Diagnostics.ProcessStartInfo
            $info.FileName = $ssh
            $allArgs = @($options) + @([string]$p.Target, $remote)
            $info.Arguments = ($allArgs | ForEach-Object { Quote-NativeArgument $_ }) -join ' '
            $info.UseShellExecute = $false
            $info.RedirectStandardOutput = $true
            $info.StandardOutputEncoding = [Text.UTF8Encoding]::new($false)
            $sshProcess = [Diagnostics.Process]::Start($info)
            $connection = if ($ConnectionId) { $ConnectionId } else { [guid]::NewGuid().ToString() }
            $readyFile = Join-Path $store ('connection-' + $connection + '.json')
            [string]$sshProcess.Id | Set-Content -LiteralPath $readyFile -Encoding ASCII
            $stopFile = Join-Path $store ('connection-' + $connection + '.stop')
            $browserReadyFile = Join-Path $store ('connection-' + $connection + '.browser-ready')
            $showFile = Join-Path $store ('connection-' + $connection + '.show')
            $opened = $false
            while ($true) {
                Update-WorkerControls
                $readTask = $sshProcess.StandardOutput.ReadLineAsync()
                while (-not $readTask.Wait(100)) {
                    Update-WorkerControls
                }
                $line = $readTask.Result
                if ($null -eq $line) { break }
                if ($line -match '^DELFIN_BROWSER:([A-Za-z0-9+/=]+)$') {
                    $url = Get-BrowserUrl $Matches[1] 0
                    $remoteUri = [Uri]$url
                    $actualLocal = $LocalPort
                    $builder = New-Object UriBuilder($url)
                    $builder.Port = $actualLocal
                    $url = $builder.Uri.AbsoluteUri
                    $handoffFile = $readyFile
                    Write-Host "Dashboard server port: $($remoteUri.Port); local tunnel port: $actualLocal"
                    $browserJob = Start-Job -ArgumentList $url,$actualLocal,$handoffFile,$browserReadyFile -ScriptBlock {
                        param($url,$local,$readyFile,$browserReadyFile)
                        $deadline = [DateTime]::UtcNow.AddMinutes(5)
                        while ([DateTime]::UtcNow -lt $deadline) {
                            try {
                                if (Test-Path -LiteralPath $readyFile) {
                                    $owner = [int](Get-Content -LiteralPath $readyFile -Raw)
                                    $connection = Get-NetTCPConnection -LocalAddress '127.0.0.1' -LocalPort $local -State Listen -ErrorAction Stop
                                    if ($connection.OwningProcess -contains $owner) {
                                        $probe = Invoke-WebRequest -Uri "http://127.0.0.1:$local/login" -UseBasicParsing -TimeoutSec 2
                                        if ($probe.StatusCode -eq 200) { Start-Process $url; [IO.File]::WriteAllText($browserReadyFile,'ready'); return }
                                    }
                                    if (-not (Get-Process -Id $owner -ErrorAction SilentlyContinue)) { throw 'SSH tunnel ended.' }
                                }
                            } catch {}
                            Start-Sleep -Seconds 2
                        }
                        Write-Output 'Browser startup timed out after 5 minutes. Check the dashboard terminal and SSH login.'
                    }
                    $url = $null
                    $opened = $true
                    break
                }
                [Console]::WriteLine($line)
            }
            if ($opened) {
                # Read one decoded character so a short final marker cannot wait for a full buffer.
                $closedMarker = 'DELFIN_CONNECTION_CLOSED:' + ([IO.Path]::GetFileName([IO.Path]::GetDirectoryName($tunnelPath)))
                $tail = ''
                $oneChar = New-Object char[] 1
                while ($true) {
                    Update-WorkerControls
                    $readTask = $sshProcess.StandardOutput.ReadAsync($oneChar,0,1)
                    while (-not $readTask.Wait(100)) {
                        Update-WorkerControls
                    }
                    if ($readTask.Result -eq 0) { break }
                    $chunk = [string]$oneChar[0]
                    [Console]::Write($chunk)
                    $tail += $chunk
                    if ($tail.Contains($closedMarker)) {
                        $remoteClosed = $true
                        # Remote cleanup has finished. Close lingering SSH channels/PTY holders.
                        if (-not $sshProcess.HasExited) { $sshProcess.Kill() }
                        break
                    }
                    if ($tail.Length -gt 512) { $tail = $tail.Substring($tail.Length-512) }
                }
            }
            $sshProcess.WaitForExit()
            Write-Host "SSH ended (exit code $($sshProcess.ExitCode))."

        }

    } catch {
        [void][DelfinConsole]::ShowWindow($loginWindow,5)
        Write-Host $_.Exception.Message -ForegroundColor Red
        if ($ConnectionId -and -not $disconnectRequested) { [void](Read-Host 'Press Enter to close') }
    }
    finally {
        if ($browserJob) { Stop-Job $browserJob; Receive-Job $browserJob; Remove-Job $browserJob }
        if ($readyFile -and (Test-Path -LiteralPath $readyFile)) { Remove-Item -LiteralPath $readyFile -Force }
        if ($null -ne $sshProcess) {
            if (-not $sshProcess.HasExited) { $sshProcess.Kill(); $sshProcess.WaitForExit() }
            if ($sshProcess.ExitCode -ne 0 -and -not $disconnectRequested -and -not $remoteClosed -and $ConnectionId) {
                [void][DelfinConsole]::ShowWindow($loginWindow,5)
                [void](Read-Host 'SSH ended. Press Enter to close')
            }
            $sshProcess.Dispose()
        }
        if ($terminalIcon) { $terminalIcon.Dispose() }
        if ($showFile) {
            $ack = [IO.Path]::ChangeExtension($showFile,'show-status')
            if (Test-Path -LiteralPath $ack) { Remove-Item -LiteralPath $ack -Force }
        }
        foreach ($file in @($stopFile,$browserReadyFile,$showFile)) {
            if ($file -and (Test-Path -LiteralPath $file)) { Remove-Item -LiteralPath $file -Force }
        }
    }
    return
}

$script:ownedConnections = New-Object Collections.ArrayList
function Disconnect-Owned([string]$profile = '') {
    foreach ($item in $script:ownedConnections) {
        if (($profile -eq '' -or $item.Profile -eq $profile) -and -not $item.Worker.HasExited) {
            [IO.File]::WriteAllText((Join-Path $store ('connection-' + $item.Connection + '.stop')),'disconnect')
        }
    }
}
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing
Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
public static class DelfinAppWindow {
    [DllImport("shell32.dll", CharSet=CharSet.Unicode)] public static extern int SetCurrentProcessExplicitAppUserModelID(string id);
    [DllImport("user32.dll")] public static extern bool IsWindowVisible(IntPtr window);
    [DllImport("user32.dll")] public static extern IntPtr GetWindow(IntPtr window, uint command);
    [DllImport("user32.dll", EntryPoint="GetWindowLongW")] public static extern int GetWindowLong(IntPtr window, int index);
}
'@
[void][DelfinAppWindow]::SetCurrentProcessExplicitAppUserModelID('ComPlat.DELFIN.Launcher')
if ($GuiSmokeTest) {
    $store = Join-Path $PSScriptRoot 'smoke-state'
    $profilesFile = Join-Path $store 'profiles.json'
}
[Windows.Forms.Application]::EnableVisualStyles()
$form = New-Object Windows.Forms.Form
$form.Text = 'DELFIN - SSH Dashboard'
$form.ShowInTaskbar = $true
$form.ShowIcon = $true
$form.ClientSize = New-Object Drawing.Size(590,480)
$form.BackColor = [Drawing.Color]::White
$form.Font = New-Object Drawing.Font('Segoe UI',9)
$form.StartPosition = 'CenterScreen'
$form.FormBorderStyle = 'FixedDialog'
$form.MaximizeBox = $false
$logo = New-Object Windows.Forms.PictureBox
$logo.Location = New-Object Drawing.Point(20,10)
$logo.Size = New-Object Drawing.Size(85,85)
$logo.SizeMode = 'Zoom'
$logo.Image = [Drawing.Image]::FromFile((Join-Path $PSScriptRoot 'DELFIN_logo.png'))
$form.Controls.Add($logo)
$form.Icon = New-Object Drawing.Icon((Join-Path $PSScriptRoot 'DELFIN.ico'))
Add-Type -Path (Join-Path $PSScriptRoot 'Taskbar.cs')
$form.add_HandleCreated({ [DelfinTaskbar]::Configure($form.Handle,(Join-Path $PSScriptRoot 'DELFIN.exe'),(Join-Path $PSScriptRoot 'DELFIN.ico')) })
$heading = New-Object Windows.Forms.Label
$heading.Text = 'DELFIN'
$heading.Font = New-Object Drawing.Font('Segoe UI',22,[Drawing.FontStyle]::Bold)
$heading.ForeColor = [Drawing.Color]::FromArgb(0,80,110)
$heading.Location = New-Object Drawing.Point(125,15)
$heading.Size = New-Object Drawing.Size(400,42)
$form.Controls.Add($heading)
$subtitle = New-Object Windows.Forms.Label
$subtitle.Text = 'SSH access, dashboard and working terminal'
$subtitle.Location = New-Object Drawing.Point(128,65)
$subtitle.Size = New-Object Drawing.Size(420,25)
$form.Controls.Add($subtitle)
$combo = New-Object Windows.Forms.ComboBox
$combo.Location = New-Object Drawing.Point(180,120)
$combo.Size = New-Object Drawing.Size(385,25)
$combo.DropDownStyle = 'DropDownList'
$form.Controls.Add($combo)
$fields = @{}
$labels = @('SSH server / user@host','Working directory (optional)','DELFIN location (optional)','Starting port (0=Auto)')
$keys = @('Target','Directory','Executable','Port')
$advancedControls = @()
for ($i=0; $i -lt $keys.Count; $i++) {
    $label = New-Object Windows.Forms.Label
    $label.Text = $labels[$i]
    $row = if ($i -eq 0) { 160 } else { 285 + $i*40 }
    $label.Location = New-Object Drawing.Point(15,$row)
    $label.Size = New-Object Drawing.Size(165,25)
    $form.Controls.Add($label)
    $text = New-Object Windows.Forms.TextBox
    $text.Location = New-Object Drawing.Point(180,$row)
    $text.Size = New-Object Drawing.Size(385,25)
    $fields[$keys[$i]] = $text
    $form.Controls.Add($text)
    if ($i -gt 0) { $advancedControls += @($label,$text) }
}
$fields.Port.Text = '0'
try { $script:profiles = @(Read-Profiles) } catch { [void][Windows.Forms.MessageBox]::Show('Could not load connections: ' + $_.Exception.Message,'DELFIN'); return }
$workingToggle = New-Object Windows.Forms.CheckBox
$workingToggle.Text = 'Working shell (Ctrl+B, then N)'
$workingToggle.Checked = $false
$workingToggle.Location = New-Object Drawing.Point(180,205)
$workingToggle.Size = New-Object Drawing.Size(380,25)
$form.Controls.Add($workingToggle)
$keepToggle = New-Object Windows.Forms.CheckBox
$keepToggle.Text = 'Keep session after disconnect'
$keepToggle.Location = New-Object Drawing.Point(180,245)
$keepToggle.Size = New-Object Drawing.Size(380,25)
$form.Controls.Add($keepToggle)
$advanced = New-Object Windows.Forms.CheckBox
$advanced.Text = 'Advanced settings'
$advanced.Location = New-Object Drawing.Point(180,285)
$advanced.Size = New-Object Drawing.Size(360,25)
$form.Controls.Add($advanced)
$script:selectedId = $null
foreach ($p in $script:profiles) { [void]$combo.Items.Add($p.Name) }
$combo.add_SelectedIndexChanged({
    if ($combo.SelectedIndex -ge 0) {
        $p = $script:profiles[$combo.SelectedIndex]
        $script:selectedId = $p.Id
        foreach ($key in $keys) { $fields[$key].Text = [string]$p.$key }
        $workingToggle.Checked = [bool]$p.OpenWorkingTerminal
        $keepToggle.Checked = [bool]$p.KeepSession
        $advanced.Checked = [bool]($p.Directory -or $p.Executable -or [int]$p.Port -ne 0)
    }
})
$note = New-Object Windows.Forms.Label
$note.Text = 'Dashboard terminal stays visible. Ctrl+C stops the dashboard. Disconnect respects Keep session.'
$note.Location = New-Object Drawing.Point(15,370)
$note.Size = New-Object Drawing.Size(555,35)
$form.Controls.Add($note)
function Save-Current {
    $id = $script:selectedId
    if (-not $id) { $id = [guid]::NewGuid().ToString() }
    $p = [pscustomobject]@{ Id=$id; Name=$fields.Target.Text.Trim(); Target=$fields.Target.Text.Trim(); Directory=$fields.Directory.Text.Trim(); Executable=$fields.Executable.Text.Trim(); Port=[int]$fields.Port.Text; KeepSession=$keepToggle.Checked; OpenWorkingTerminal=$workingToggle.Checked }
    Check-Profile $p
    $script:profiles = @($script:profiles | Where-Object { $_.Id -ne $id }) + @($p)
    [void][IO.Directory]::CreateDirectory($store)
    ConvertTo-Json -InputObject @($script:profiles) | Set-Content -LiteralPath $profilesFile -Encoding UTF8
    $script:selectedId = $id
    $combo.Items.Clear()
    foreach ($item in $script:profiles) { [void]$combo.Items.Add($item.Name) }
    $combo.SelectedIndex = $script:profiles.Count - 1
    return $id
}
$actionButtons = @()
foreach ($spec in @(@('New',15),@('Save',130),@('Start / Reconnect',280))) {
    $button = New-Object Windows.Forms.Button
    $button.Text = $spec[0]
    $button.Location = New-Object Drawing.Point($spec[1],420)
    $button.Size = New-Object Drawing.Size($(if ($spec[1] -eq 280) { 285 } else { 105 }),35)
    if ($spec[0] -eq 'New') {
        $button.add_Click({ $script:selectedId=$null; $combo.SelectedIndex=-1; foreach ($key in $keys) { $fields[$key].Text='' }; $fields.Port.Text='0'; $advanced.Checked=$false; $workingToggle.Checked=$false; $keepToggle.Checked=$false })
    } elseif ($spec[0] -eq 'Save') {
        $button.add_Click({ try { [void](Save-Current) } catch { [void][Windows.Forms.MessageBox]::Show($_.Exception.Message,'DELFIN') } })
    } else {
        $button.add_Click({
            try {
                $active = @($script:ownedConnections | Where-Object { -not $_.Worker.HasExited })
                if ($active.Count) {
                    $script:connectionStatus.Text = 'A connection is already active. Use the dashboard browser, or Disconnect first.'
                    return
                }
                $id=Save-Current
                Start-Window 'Dashboard' $id
            } catch { [void][Windows.Forms.MessageBox]::Show($_.Exception.Message,'DELFIN') }
        })
    }
    $form.Controls.Add($button)
    $actionButtons += $button
}
$connectionButtons = @()
foreach ($label in @('Show terminal','Disconnect')) {
    $control = New-Object Windows.Forms.Button
    $control.Text = $label
    $control.Size = New-Object Drawing.Size(160,30)
    $control.Left = if ($label -eq 'Show terminal') { 15 } else { 190 }
    if ($label -eq 'Disconnect') {
        $control.add_Click({ Disconnect-Owned; $script:connectionStatus.Text = 'Disconnect requested...' })
    } else {
        $control.add_Click({
            $active = @($script:ownedConnections | Where-Object { -not $_.Worker.HasExited })
            if (-not $active.Count) { $script:connectionStatus.Text = 'No active SSH terminal to show.'; return }
            foreach ($item in $active) {
                [IO.File]::WriteAllText((Join-Path $store ('connection-' + $item.Connection + '.show')),'show')
            }
            $script:connectionStatus.Text = 'Terminal show request sent...'
        })
    }
    $form.Controls.Add($control)
    $connectionButtons += $control
}
$script:connectionStatus = New-Object Windows.Forms.Label
$script:connectionStatus.Left = 15
$script:connectionStatus.Size = New-Object Drawing.Size(555,30)
$script:connectionStatus.Text = 'Not connected.'
$form.Controls.Add($script:connectionStatus)
$statusTimer = New-Object Windows.Forms.Timer
$statusTimer.Interval = 500
$statusTimer.add_Tick({
    $active = @($script:ownedConnections | Where-Object { -not $_.Worker.HasExited })
    if ($active.Count -eq 0) {
        $script:connectionStatus.Text = if ($script:ownedConnections.Count) { 'Disconnected. Keep ON dashboards can remain on the server.' } else { 'Not connected.' }
    } else {
        $stopping = @($active | Where-Object { Test-Path -LiteralPath (Join-Path $store ('connection-' + $_.Connection + '.stop')) })
        $ready = @($active | Where-Object { Test-Path -LiteralPath (Join-Path $store ('connection-' + $_.Connection + '.browser-ready')) })
        if ($stopping.Count) { $script:connectionStatus.Text = 'Disconnecting SSH...' }
        elseif ($ready.Count) { $script:connectionStatus.Text = 'Connected. Dashboard SSH terminal stays open for Ctrl+C.' }
        else { $script:connectionStatus.Text = 'Connecting. Complete password/OTP in the SSH window.' }
        foreach ($item in $active) {
            $ack = Join-Path $store ('connection-' + $item.Connection + '.show-status')
            try {
                if ((Test-Path -LiteralPath $ack) -and ([DateTime]::UtcNow - (Get-Item -LiteralPath $ack).LastWriteTimeUtc).TotalSeconds -lt 8) {
                    $script:connectionStatus.Text = Get-Content -LiteralPath $ack -Raw
                }
            } catch { } # The worker can remove the acknowledgement while disconnecting.
        }
    }
})
$statusTimer.Start()
$form.add_FormClosing({ Disconnect-Owned; $statusTimer.Stop() })
function Show-Advanced {
    foreach ($control in $advancedControls) { $control.Visible = $advanced.Checked }
    $note.Location = New-Object Drawing.Point(15,$(if ($advanced.Checked) { 450 } else { 325 }))
    foreach ($button in $actionButtons) { $button.Top = $(if ($advanced.Checked) { 500 } else { 375 }) }
    foreach ($control in $connectionButtons) { $control.Top = $(if ($advanced.Checked) { 545 } else { 420 }) }
    $script:connectionStatus.Top = $(if ($advanced.Checked) { 580 } else { 455 })
    $form.ClientSize = New-Object Drawing.Size(590,$(if ($advanced.Checked) { 625 } else { 500 }))
}
$advanced.add_CheckedChanged({ Show-Advanced })
Show-Advanced
if ($combo.Items.Count -gt 0) { $combo.SelectedIndex=0 }
if ($GuiSmokeTest) {
    $smokeTimer = New-Object Windows.Forms.Timer
    $smokeTimer.Interval = 500
    $smokeTimer.add_Tick({
        $smokeTimer.Stop()
        $report = @{ Visible=[DelfinAppWindow]::IsWindowVisible($form.Handle); ShowInTaskbar=$form.ShowInTaskbar; Owner=[DelfinAppWindow]::GetWindow($form.Handle,4).ToInt64(); ExtendedStyle=[DelfinAppWindow]::GetWindowLong($form.Handle,-20); HasIcon=($null -ne $form.Icon); RelaunchCommand=[DelfinTaskbar]::Read($form.Handle,2); RelaunchIcon=[DelfinTaskbar]::Read($form.Handle,3); RelaunchName=[DelfinTaskbar]::Read($form.Handle,4); AppId=[DelfinTaskbar]::Read($form.Handle,5) }
        $report | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $PSScriptRoot 'gui-smoke.json') -Encoding UTF8
        $form.Close()
    })
    $form.add_Shown({ $smokeTimer.Start() })
}
try { [Windows.Forms.Application]::Run($form) } finally { $logo.Image.Dispose(); $form.Icon.Dispose(); $form.Dispose() }
