#requires -Version 5.1
[CmdletBinding()]
param([ValidateSet('App','Dashboard','Terminal','Validate')][string]$Mode = 'App', [string]$ProfileId, [int]$LocalPort = 0, [int]$RemotePort = 0, [string]$ConnectionId, [switch]$GuiSmokeTest)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$store = Join-Path $env:LOCALAPPDATA 'DELFIN Launcher'
$profilesFile = Join-Path $store 'profiles.json'

function Normalize-Profile($value) {
    $defaults = [ordered]@{ Id=[guid]::NewGuid().ToString(); Name=''; Target=''; Directory=''; Executable=''; Port=0; KeepSession=$false; OpenWorkingTerminal=$false; RememberLogin=$true }
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
function Get-PersonName {
    # Who sits at this Windows account. A server that runs several people on one
    # account (remote_bootstrap.sh, ~/.config/delfin/launcher-people.conf) gives
    # each a home of their own; every other server ignores it.
    $name = [regex]::Replace(([string]$env:USERNAME).ToLowerInvariant(), '[^a-z0-9_-]', '-').Trim('-')
    if ($name.Length -gt 32) { $name = $name.Substring(0, 32).TrimEnd('-') }
    return $name
}
function Get-LoginKey {
    # Remember login: one ed25519 key per Windows user, without a passphrase,
    # protected by the Windows account. remote_bootstrap.sh adds its public
    # half to the server on the first (password) login. $null when OpenSSH's
    # ssh-keygen is missing or the key cannot be made: the login stays as is.
    $folder = Join-Path $env:USERPROFILE '.ssh'
    $key = Join-Path $folder 'delfin_launcher_ed25519'
    if (-not (Test-Path -LiteralPath $key)) {
        $keygen = Join-Path $env:SystemRoot 'System32\OpenSSH\ssh-keygen.exe'
        if (-not (Test-Path -LiteralPath $keygen)) { return $null }
        [void][IO.Directory]::CreateDirectory($folder)
        $info = New-Object Diagnostics.ProcessStartInfo
        $info.FileName = $keygen
        $info.Arguments = (@('-q','-t','ed25519','-N','','-C',('delfin-launcher-' + (Get-PersonName)),'-f',$key) | ForEach-Object { Quote-NativeArgument $_ }) -join ' '
        $info.UseShellExecute = $false
        $made = [Diagnostics.Process]::Start($info)
        $made.WaitForExit()
        if ($made.ExitCode -ne 0) { return $null }
    }
    $public = $key + '.pub'
    if (-not (Test-Path -LiteralPath $public)) { return $null }
    $line = (Get-Content -LiteralPath $public -Raw).Trim()
    if ($line -notmatch '^ssh-ed25519 [A-Za-z0-9+/]+={0,2}( [A-Za-z0-9@._-]{1,64})?$') { return $null }
    return [pscustomobject]@{ Key=$key; Public=$line }
}
function Get-ScriptVersion {
    # The server half (remote_bootstrap.sh, remote_launcher.py) is copied as
    # plain files into ~/.cache/delfin-launcher/<version>/ on the server and
    # run from there. <version> is the first 16 hex digits of each file's
    # SHA-256, so a changed script gets a folder of its own; the full digests
    # are checked on the server before each start.
    $digests = @{}
    foreach ($name in @('remote_bootstrap.sh','remote_launcher.py')) {
        $digests[$name] = (Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $PSScriptRoot $name)).Hash.ToLowerInvariant()
    }
    $folder = '.cache/delfin-launcher/' + $digests['remote_bootstrap.sh'].Substring(0,16) + $digests['remote_launcher.py'].Substring(0,16)
    return [pscustomobject]@{ Folder=$folder; Bootstrap=$digests['remote_bootstrap.sh']; Launcher=$digests['remote_launcher.py'] }
}
function Get-ScriptCheckCommand($version) {
    # Exit 0: both files are present with these digests. Exit 10: missing or
    # different (the folder now exists, owner-only). Exit 11: no folder.
    # $1 is the folder below HOME, $2 and $3 the digests (hex only).
    $checkScript = 'umask 077 && cd && mkdir -p $1 && cd $1 || exit 11; printf ''%s  %s\n'' $2 remote_bootstrap.sh $3 remote_launcher.py | sha256sum -c --status 2>/dev/null || exit 10'
    if ($version.Folder -notmatch '^\.cache/delfin-launcher/[a-f0-9]{32}$' -or $version.Bootstrap -notmatch '^[a-f0-9]{64}$' -or $version.Launcher -notmatch '^[a-f0-9]{64}$') { throw 'Invalid launcher script version.' }
    return 'bash -c ' + (Shell-Quote $checkScript) + ' delfin ' + $version.Folder + ' ' + $version.Bootstrap + ' ' + $version.Launcher
}
function Get-RemoteCommand($p, [string]$role, [string]$tunnelPath = '', [string]$publicKey = '', $version) {
    # The login shell sources the copied bootstrap; its positional parameters
    # are the arguments after 'delfin'. Only file paths and data travel.
    if ($version.Folder -notmatch '^\.cache/delfin-launcher/[a-f0-9]{32}$') { throw 'Invalid launcher script version.' }
    $bootstrap = '. ~/' + $version.Folder + '/remote_bootstrap.sh'
    $launcher = '~/' + $version.Folder + '/remote_launcher.py'
    return 'bash -ilc ' + (Shell-Quote $bootstrap) + ' delfin ' + [int]$p.Port + ' ' + (Shell-Quote $p.Directory) + ' ' + (Shell-Quote $p.Executable) + ' ' + (Shell-Quote $launcher) + ' ' + (Shell-Quote $role) + ' ' + (Shell-Quote $tunnelPath) + ' ' + $(if ($p.KeepSession) { '1' } else { '0' }) + ' ' + $(if ($p.OpenWorkingTerminal) { '1' } else { '0' }) + ' ' + (Shell-Quote (Get-PersonName)) + ' ' + (Shell-Quote $publicKey)
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
    $options = @('-o','StrictHostKeyChecking=ask','-o','ForwardAgent=no','-o','ForwardX11=no','-o',$forwardingPolicy,'-o','PermitLocalCommand=no','-o','ServerAliveInterval=30','-o','ServerAliveCountMax=3','-o','ControlMaster=no','-o','ControlPath=none')
    # ssh sessions get a terminal (-tt); the short check runs without one (-T);
    # scp ('Copy') takes neither.
    if ($kind -eq 'Dashboard' -or $kind -eq 'Terminal') { $options = @('-tt') + $options } elseif ($kind -eq 'Check') { $options = @('-T') + $options }
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
function Start-Window([string]$kind, [string]$id, [string]$title, [string]$connection = '') {
    # The dashboard connection runs in a console window of its own: a plain
    # Windows PowerShell started on this script, where OpenSSH asks for the
    # password and OTP.
    if ($PSCommandPath.Contains('"')) { throw 'Invalid installation path.' }
    if ($kind -ne 'Dashboard') { throw 'Only the dashboard connection is started from the app.' }
    if (-not $connection) { $connection = [guid]::NewGuid().ToString() }
    $powershell = Join-Path $env:SystemRoot 'System32\WindowsPowerShell\v1.0\powershell.exe'
    $arguments = (@('-NoProfile','-ExecutionPolicy','RemoteSigned','-File',$PSCommandPath,'-Mode','Dashboard','-ProfileId',$id,'-ConnectionId',$connection) | ForEach-Object { Quote-NativeArgument $_ }) -join ' '
    $worker = Start-Process -FilePath $powershell -ArgumentList $arguments -WorkingDirectory $PSScriptRoot -PassThru
    [void]$script:ownedConnections.Add([pscustomobject]@{ Profile=$id; Connection=$connection; Worker=$worker; Title=$title })
}
function Invoke-OpenSSH([string]$program, [string[]]$arguments) {
    # One OpenSSH client run in this console (password and host-key prompts
    # included); a Disconnect from the app ends it. Returns the exit code.
    $info = New-Object Diagnostics.ProcessStartInfo
    $info.FileName = $program
    $info.Arguments = ($arguments | ForEach-Object { Quote-NativeArgument $_ }) -join ' '
    $info.WorkingDirectory = $PSScriptRoot
    $info.UseShellExecute = $false
    $script:sshProcess = [Diagnostics.Process]::Start($info)
    try {
        while (-not $script:sshProcess.WaitForExit(100)) { Update-WorkerControls }
        if ($script:disconnectRequested) { throw 'Disconnected.' }
        return $script:sshProcess.ExitCode
    } finally { $script:sshProcess.Dispose(); $script:sshProcess = $null }
}
function Sync-ServerScripts([string]$target, [string[]]$identity) {
    # A short SSH call checks the server's copy; only a missing or different
    # copy is sent, with Windows OpenSSH's scp, as two plain files.
    $version = Get-ScriptVersion
    $ssh = Join-Path $env:SystemRoot 'System32\OpenSSH\ssh.exe'
    $scp = Join-Path $env:SystemRoot 'System32\OpenSSH\scp.exe'
    Write-Host 'Checking the DELFIN launcher scripts on the server...'
    $code = Invoke-OpenSSH $ssh (@(Get-SSHOptions 'Check') + $identity + @($target, (Get-ScriptCheckCommand $version)))
    if ($code -eq 0) { return $version }
    if ($code -ne 10) { throw "Could not check the launcher scripts on the server (exit code $code)." }
    if (-not (Test-Path -LiteralPath $scp)) { throw 'Windows OpenSSH scp.exe is missing. See README.' }
    Write-Host 'Copying the DELFIN launcher scripts to the server (once per launcher version).'
    Write-Host 'Without a remembered login, OpenSSH asks for the password and OTP again for this copy.'
    $code = Invoke-OpenSSH $scp (@('-q') + @(Get-SSHOptions 'Copy') + $identity + @('remote_bootstrap.sh','remote_launcher.py',($target + ':' + $version.Folder + '/')))
    if ($code -ne 0) { throw "Could not copy the launcher scripts to the server (exit code $code)." }
    return $version
}

if ($Mode -eq 'Validate') { return }

if ($Mode -ne 'App') {
    # tmux draws with ANSI control sequences; Windows PowerShell 5.1 and
    # Windows Terminal render them in their console.
    $env:TERM = 'xterm-256color'
    function Update-WorkerControls {
        if ($stopFile -and (Test-Path -LiteralPath $stopFile)) {
            $script:disconnectRequested = $true
            if ($null -ne $script:sshProcess -and -not $script:sshProcess.HasExited) { $script:sshProcess.Kill() }
        }
    }
    $sshProcess = $null
    $browserJob = $null
    $readyFile = $null
    $stopFile = $null
    $browserReadyFile = $null
    $disconnectRequested = $false
    $remoteClosed = $false
    try {
        $matchingProfiles = @(Read-Profiles | Where-Object { $_.Id -eq $ProfileId })
        if ($matchingProfiles.Count -ne 1) { throw 'Connection not found.' }
        $p = $matchingProfiles[0]
        Check-Profile $p
        $ssh = Join-Path $env:SystemRoot 'System32\OpenSSH\ssh.exe'
        if (-not (Test-Path $ssh)) { throw 'Windows OpenSSH Client is missing. See README.' }
        $Host.UI.RawUI.WindowTitle = "SSH $Mode - $($p.Name)"
        if ($ConnectionId -and ($ConnectionId -notmatch '^[a-f0-9-]{36}$' -or ($Mode -eq 'Terminal' -and ($LocalPort -lt 1024 -or $RemotePort -lt 1024)))) { throw 'Invalid tunnel request.' }
        $connection = if ($ConnectionId) { $ConnectionId } else { [guid]::NewGuid().ToString() }
        $readyFile = Join-Path $store ('connection-' + $connection + '.json')
        $stopFile = Join-Path $store ('connection-' + $connection + '.stop')
        $tunnelPath = ''
        if ($Mode -eq 'Dashboard') {
            $LocalPort = Get-FreeLocalPort ([int]$p.Port)
            $tunnelPath = '/tmp/delfin-win-tunnel-' + [guid]::NewGuid().ToString() + '/dashboard.sock'
        }
        $options = @(Get-SSHOptions $Mode $LocalPort $RemotePort $tunnelPath)
        $loginKey = $null
        if ($p.RememberLogin) { try { $loginKey = Get-LoginKey } catch { $loginKey = $null } }
        $publicKey = ''
        $identity = @()
        if ($loginKey) { $options += @('-i',$loginKey.Key); $identity = @('-i',$loginKey.Key); $publicKey = $loginKey.Public }
        if ($Mode -eq 'Terminal') {
            Write-Host 'Working terminal. Password and OTP prompts are handled exclusively by OpenSSH.'
            Write-Host 'To paste a password here: right-click in this window. Ctrl+V does not paste into this console; it sends a control character, and the server sees a wrong password.'
            $version = Sync-ServerScripts ([string]$p.Target) $identity
            $terminalCommand = Get-RemoteCommand $p 'terminal' '' $publicKey $version
            $terminalInfo = New-Object Diagnostics.ProcessStartInfo
            $terminalInfo.FileName = $ssh
            $terminalArgs = @($options) + @([string]$p.Target,$terminalCommand)
            $terminalInfo.Arguments = ($terminalArgs | ForEach-Object { Quote-NativeArgument $_ }) -join ' '
            $terminalInfo.UseShellExecute = $false
            $sshProcess = [Diagnostics.Process]::Start($terminalInfo)
            [string]$sshProcess.Id | Set-Content -LiteralPath $readyFile -Encoding ASCII
            $sshProcess.WaitForExit()
            Write-Host "SSH ended (exit code $($sshProcess.ExitCode))."
        } else {
            Write-Host 'Dashboard terminal: enter your password and OTP directly in OpenSSH.'
            Write-Host 'To paste a password here: right-click in this window. Ctrl+V does not paste into this console; it sends a control character, and the server sees a wrong password.'
            Write-Host 'On first login, compare the host fingerprint with the official server fingerprint.'
            if ($loginKey) { Write-Host 'Remember login ON: after this login the server knows this PC; the next connection needs no password (where the server allows keys).' }
            $version = Sync-ServerScripts ([string]$p.Target) $identity
            $remote = Get-RemoteCommand $p 'dashboard' $tunnelPath $publicKey $version
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
            [string]$sshProcess.Id | Set-Content -LiteralPath $readyFile -Encoding ASCII
            $browserReadyFile = Join-Path $store ('connection-' + $connection + '.browser-ready')
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
        Write-Host $_.Exception.Message -ForegroundColor Red
        if ($ConnectionId -and -not $disconnectRequested) { [void](Read-Host 'Press Enter to close') }
    }
    finally {
        if ($browserJob) { Stop-Job $browserJob; Receive-Job $browserJob; Remove-Job $browserJob }
        if ($readyFile -and (Test-Path -LiteralPath $readyFile)) { Remove-Item -LiteralPath $readyFile -Force }
        if ($null -ne $sshProcess) {
            if (-not $sshProcess.HasExited) { $sshProcess.Kill(); $sshProcess.WaitForExit() }
            if ($sshProcess.ExitCode -ne 0 -and -not $disconnectRequested -and -not $remoteClosed -and $ConnectionId) {
                [void](Read-Host 'SSH ended. Press Enter to close')
            }
            $sshProcess.Dispose()
        }
        foreach ($file in @($stopFile,$browserReadyFile)) {
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
if ($GuiSmokeTest) {
    $store = Join-Path $PSScriptRoot 'smoke-state'
    $profilesFile = Join-Path $store 'profiles.json'
}
# One settings app per Windows login: a second start brings the first window
# to the front (WScript.Shell.AppActivate on its title) and ends.
$appMutex = $null
if (-not $GuiSmokeTest) {
    $created = $false
    $appMutex = [Threading.Mutex]::new($true, ('Local\DELFIN_Launcher_' + [Security.Principal.WindowsIdentity]::GetCurrent().User.Value), [ref]$created)
    if (-not $created) {
        $appMutex.Dispose()
        if (-not (New-Object -ComObject WScript.Shell).AppActivate('DELFIN - SSH Dashboard')) {
            [void][Windows.Forms.MessageBox]::Show('DELFIN is already running or starting. Use its existing settings window.','DELFIN')
        }
        return
    }
}
# The app's console window starts minimized (the shortcut's window style);
# an error here would close it unseen, so it is shown in a message box.
try {
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
    $rememberToggle = New-Object Windows.Forms.CheckBox
    $rememberToggle.Text = 'Remember login on this PC (SSH key)'
    $rememberToggle.Checked = $true
    $rememberToggle.Location = New-Object Drawing.Point(180,445)
    $rememberToggle.Size = New-Object Drawing.Size(380,25)
    $form.Controls.Add($rememberToggle)
    $advancedControls += @($rememberToggle)
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
            $rememberToggle.Checked = [bool]$p.RememberLogin
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
        $p = [pscustomobject]@{ Id=$id; Name=$fields.Target.Text.Trim(); Target=$fields.Target.Text.Trim(); Directory=$fields.Directory.Text.Trim(); Executable=$fields.Executable.Text.Trim(); Port=[int]$fields.Port.Text; KeepSession=$keepToggle.Checked; OpenWorkingTerminal=$workingToggle.Checked; RememberLogin=$rememberToggle.Checked }
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
            $button.add_Click({ $script:selectedId=$null; $combo.SelectedIndex=-1; foreach ($key in $keys) { $fields[$key].Text='' }; $fields.Port.Text='0'; $advanced.Checked=$false; $workingToggle.Checked=$false; $keepToggle.Checked=$false; $rememberToggle.Checked=$true })
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
                    Start-Window 'Dashboard' $id ('SSH Dashboard - ' + $fields.Target.Text.Trim())
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
                # WScript.Shell.AppActivate: the worker's console by process ID,
                # else by its window title (Windows Terminal hosts it elsewhere).
                $shell = New-Object -ComObject WScript.Shell
                $shown = $false
                foreach ($item in $active) {
                    if ($shell.AppActivate($item.Worker.Id) -or $shell.AppActivate($item.Title)) { $shown = $true }
                }
                $script:notice = if ($shown) { 'Terminal brought to the front. If it is minimized, click it in the taskbar.' } else { 'Windows could not bring the SSH terminal to the front; find it in the taskbar.' }
                $script:noticeUntil = [DateTime]::UtcNow.AddSeconds(8)
            })
        }
        $form.Controls.Add($control)
        $connectionButtons += $control
    }
    $script:connectionStatus = New-Object Windows.Forms.Label
    $script:connectionStatus.Left = 15
    $script:connectionStatus.Size = New-Object Drawing.Size(555,30)
    $script:connectionStatus.Text = 'Not connected.'
    $script:notice = ''
    $script:noticeUntil = [DateTime]::MinValue
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
            if ([DateTime]::UtcNow -lt $script:noticeUntil) { $script:connectionStatus.Text = $script:notice }
        }
    })
    $statusTimer.Start()
    $form.add_FormClosing({ Disconnect-Owned; $statusTimer.Stop() })
    function Show-Advanced {
        foreach ($control in $advancedControls) { $control.Visible = $advanced.Checked }
        $note.Location = New-Object Drawing.Point(15,$(if ($advanced.Checked) { 490 } else { 325 }))
        foreach ($button in $actionButtons) { $button.Top = $(if ($advanced.Checked) { 540 } else { 375 }) }
        foreach ($control in $connectionButtons) { $control.Top = $(if ($advanced.Checked) { 585 } else { 420 }) }
        $script:connectionStatus.Top = $(if ($advanced.Checked) { 620 } else { 455 })
        $form.ClientSize = New-Object Drawing.Size(590,$(if ($advanced.Checked) { 665 } else { 500 }))
    }
    $advanced.add_CheckedChanged({ Show-Advanced })
    Show-Advanced
    if ($combo.Items.Count -gt 0) { $combo.SelectedIndex=0 }
    if ($GuiSmokeTest) {
        $smokeTimer = New-Object Windows.Forms.Timer
        $smokeTimer.Interval = 500
        $smokeTimer.add_Tick({
            $smokeTimer.Stop()
            $report = @{ Visible=$form.Visible; ShowInTaskbar=$form.ShowInTaskbar; HasOwner=($null -ne $form.Owner); TopLevel=$form.TopLevel; HasIcon=($null -ne $form.Icon); Title=$form.Text }
            $report | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $PSScriptRoot 'gui-smoke.json') -Encoding UTF8
            $form.Close()
        })
        $form.add_Shown({ $smokeTimer.Start() })
    }
    try { [Windows.Forms.Application]::Run($form) } finally { $logo.Image.Dispose(); $form.Icon.Dispose(); $form.Dispose() }
} catch {
    [void][Windows.Forms.MessageBox]::Show($_.Exception.Message,'DELFIN startup failed',[Windows.Forms.MessageBoxButtons]::OK,[Windows.Forms.MessageBoxIcon]::Error)
    exit 1
} finally {
    if ($appMutex) { $appMutex.ReleaseMutex(); $appMutex.Dispose() }
}
