<#
.SYNOPSIS
  Build a ready-to-boot DELFIN agent card (SD card or SSD) for a Raspberry Pi 5.

.DESCRIPTION
  Runs on Windows. One command does everything:
    1. download the latest Raspberry Pi OS (64-bit, desktop) + verify SHA-256
    2. write it raw to the disk via write_image.ps1 (one UAC prompt)
    3. copy the current DELFIN working tree (LF endings) + setup scripts onto
       the boot partition and write a cloud-init user-data
  On first boot (LAN needed) the Pi provisions itself: user 'pi' (SSH key
  only, desktop autologin), agent user 'delfin' with passwordless sudo,
  Python 3.11 venv with all DELFIN extras, Claude/Codex CLIs, dashboard
  service on 127.0.0.1:8866. Progress: /var/log/delfin-firstboot.log
  Later updates: scripts/pi/deploy_pi.sh delfin@<host>

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File scripts\pi\make_pi_card.ps1 -DiskNumber 1
.EXAMPLE
  # card already flashed, only refresh DELFIN payload + cloud-init
  powershell -ExecutionPolicy Bypass -File scripts\pi\make_pi_card.ps1 -DiskNumber 1 -SkipFlash
#>
param(
    [Parameter(Mandatory = $true)][int]$DiskNumber,
    [string]$Hostname = "delfin-pi",
    [string]$LoginUser = "pi",
    [string]$Timezone = "Europe/Berlin",
    [string]$KeyboardLayout = "de",
    [string]$PubKeyFile = "$env:USERPROFILE\.ssh\id_ed25519.pub",
    [string]$ImageUrl = "https://downloads.raspberrypi.com/raspios_arm64_latest",
    [string]$CacheDir = "$env:LOCALAPPDATA\delfin-pi",
    [switch]$SkipFlash,
    [switch]$Yes
)
$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$repo = (git -C $here rev-parse --show-toplevel).Trim()

function Write-LfFile([string]$path, [string]$text) {
    [IO.File]::WriteAllText($path, ($text -replace "`r`n", "`n"), (New-Object Text.UTF8Encoding $false))
}

# ---- 0. target disk sanity ---------------------------------------------------
$disk = Get-Disk -Number $DiskNumber
if ($disk.IsSystem -or $disk.IsBoot -or $disk.BusType -notin @("USB", "SD", "MMC")) {
    throw "Disk $DiskNumber ($($disk.FriendlyName), $($disk.BusType)) is not a removable USB/SD disk - refusing."
}
if (-not (Test-Path $PubKeyFile)) { throw "SSH public key not found: $PubKeyFile" }
$pubKey = (Get-Content $PubKeyFile -Raw).Trim()

if (-not $SkipFlash) {
    $gb = [math]::Round($disk.Size / 1GB, 1)
    Write-Host "Disk $DiskNumber : $($disk.FriendlyName) $gb GB ($($disk.BusType)) will be ERASED." -ForegroundColor Yellow
    if (-not $Yes -and (Read-Host "Type YES to continue") -ne "YES") { throw "aborted" }

    # ---- 1. image ------------------------------------------------------------
    New-Item -ItemType Directory -Force $CacheDir | Out-Null
    $url = (& curl.exe -sIL -o NUL -w "%{url_effective}" $ImageUrl).Trim()
    $img = Join-Path $CacheDir ([IO.Path]::GetFileName($url))
    $want = ((& curl.exe -sL "$url.sha256") -split "\s+")[0].ToLower()
    if (-not (Test-Path $img) -or (Get-FileHash $img -Algorithm SHA256).Hash.ToLower() -ne $want) {
        Write-Host "==> downloading $url"
        & curl.exe -L --fail -o $img $url
        if ($LASTEXITCODE -ne 0) { throw "download failed" }
    }
    if ((Get-FileHash $img -Algorithm SHA256).Hash.ToLower() -ne $want) { throw "SHA-256 mismatch for $img" }
    Write-Host "==> image ok: $img"

    # ---- 2. flash (elevated raw write; rpi-imager --cli hangs on Windows) -----
    $xz = @("$env:ProgramFiles\Git\mingw64\bin\xz.exe", "$env:ProgramFiles\Git\usr\bin\xz.exe") |
        Where-Object { Test-Path $_ } | Select-Object -First 1
    if (-not $xz) { throw "xz.exe not found - install Git for Windows" }
    $log = Join-Path $CacheDir "write_image.log"
    Remove-Item $log -ErrorAction SilentlyContinue
    Write-Host "==> flashing (confirm the UAC prompt, takes a few minutes) ..."
    $argLine = "-NoProfile -ExecutionPolicy Bypass -File `"$here\write_image.ps1`" -DiskNumber $DiskNumber " +
               "-Image `"$img`" -ExpectedSize $($disk.Size) -Xz `"$xz`" -Log `"$log`""
    $p = Start-Process powershell -ArgumentList $argLine -Verb RunAs -Wait -PassThru -WindowStyle Hidden
    Get-Content $log -ErrorAction SilentlyContinue | Select-Object -Last 3 | ForEach-Object { Write-Host "    $_" }
    if ($p.ExitCode -ne 0) { throw "flashing failed, see $log" }
}

# ---- 3. boot partition -------------------------------------------------------
$boot = $null
for ($i = 0; $i -lt 30 -and -not $boot; $i++) {
    try { Update-HostStorageCache -ErrorAction Stop } catch { }   # needs admin; best-effort
    $boot = Get-Partition -DiskNumber $DiskNumber |
        Where-Object { $_.DriveLetter -and (Test-Path "$($_.DriveLetter):\config.txt") } |
        Select-Object -First 1
    if (-not $boot) { Start-Sleep -Seconds 2 }
}
if (-not $boot) { throw "boot partition (config.txt) of disk $DiskNumber not mounted - replug the card and rerun with -SkipFlash" }
$bootDir = "$($boot.DriveLetter):\"
Write-Host "==> boot partition: $bootDir"

# ---- 4. DELFIN payload ---------------------------------------------------------
$dest = Join-Path $bootDir "delfin"
New-Item -ItemType Directory -Force $dest | Out-Null
Push-Location $repo
try {
    $tmpIndex = [IO.Path]::GetTempFileName()
    Copy-Item (git rev-parse --git-path index).Trim() $tmpIndex -Force
    $env:GIT_INDEX_FILE = $tmpIndex
    git -c core.safecrlf=false add -A
    $tree = (git write-tree).Trim()
    Remove-Item Env:GIT_INDEX_FILE
    Remove-Item $tmpIndex
    $version = (git rev-parse --short HEAD).Trim()
    if ($tree -ne (git rev-parse "HEAD^{tree}").Trim()) { $version += "-dirty" }
    git -c core.autocrlf=false -c core.eol=lf archive --format=tar.gz -o (Join-Path $dest "delfin-src.tgz") $tree
    if ($LASTEXITCODE -ne 0) { throw "git archive failed" }
} finally {
    Remove-Item Env:GIT_INDEX_FILE -ErrorAction SilentlyContinue
    Pop-Location
}
Write-LfFile (Join-Path $dest "VERSION") "$version`n"
foreach ($f in "pi_bootstrap.sh", "pi_install.sh", "pi_firstboot.sh") {
    Write-LfFile (Join-Path $dest $f) (Get-Content (Join-Path $here $f) -Raw)
}

# ---- 5. cloud-init -------------------------------------------------------------
$userData = @"
#cloud-config
hostname: $Hostname
manage_etc_hosts: true
timezone: $Timezone
keyboard:
  model: pc105
  layout: $KeyboardLayout
users:
  - name: $LoginUser
    groups: users,adm,dialout,audio,netdev,video,plugdev,cdrom,games,input,gpio,spi,i2c,render,sudo
    shell: /bin/bash
    lock_passwd: true
    sudo: ALL=(ALL) NOPASSWD:ALL
    ssh_authorized_keys:
      - $pubKey
ssh_pwauth: false
runcmd:
  - [systemctl, enable, --now, ssh]
  - [raspi-config, nonint, do_boot_behaviour, B4]
  - [systemd-run, --unit, delfin-firstboot, --setenv=LOGIN_USER=$LoginUser, bash, /boot/firmware/delfin/pi_firstboot.sh]

"@
Write-LfFile (Join-Path $bootDir "user-data") $userData

Write-Host ""
Write-Host "Card ready: DELFIN $version, host $Hostname, login $LoginUser (SSH key only)." -ForegroundColor Green
Write-Host "Insert into the Pi, connect LAN, power on. First boot installs everything (~20-40 min)."
Write-Host "  progress:  ssh $LoginUser@$Hostname.local tail -f /var/log/delfin-firstboot.log"
Write-Host "  password:  sudo passwd $LoginUser   (on the Pi, for local logins)"
Write-Host "  updates:   scripts/pi/deploy_pi.sh delfin@$Hostname.local"
