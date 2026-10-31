<#
.SYNOPSIS
  Write a .img.xz image raw onto a removable disk (elevated helper of make_pi_card.ps1).

.DESCRIPTION
  Streams `xz -dc <image>` straight onto \\.\PhysicalDrive<N>: no temporary
  .img, no Raspberry Pi Imager. Needs admin rights. Refuses system/boot disks,
  non-USB/SD disks and disks whose size differs from -ExpectedSize, so the
  caller's confirmation cannot hit a different disk after a replug.
#>
param(
    [Parameter(Mandatory = $true)][int]$DiskNumber,
    [Parameter(Mandatory = $true)][string]$Image,
    [Parameter(Mandatory = $true)][long]$ExpectedSize,
    [Parameter(Mandatory = $true)][string]$Xz,
    [Parameter(Mandatory = $true)][string]$Log
)
$ErrorActionPreference = "Stop"
function Log([string]$m) { Add-Content -Path $Log -Value "$(Get-Date -Format HH:mm:ss) $m" }

try {
    $disk = Get-Disk -Number $DiskNumber
    if ($disk.IsSystem -or $disk.IsBoot -or $disk.BusType -notin @("USB", "SD", "MMC") -or $disk.Size -ne $ExpectedSize) {
        throw "disk $DiskNumber failed safety checks ($($disk.FriendlyName), $($disk.BusType), $($disk.Size))"
    }
    # A hung Raspberry Pi Imager keeps the disk locked.
    Get-Process rpi-imager -ErrorAction SilentlyContinue | Stop-Process -Force
    Start-Sleep -Seconds 1

    # Clear-Disk silently leaves removable media alone; diskpart clean does not.
    Log "clearing disk $DiskNumber"
    $dp = Join-Path $env:TEMP "delfin-diskpart.txt"
    Set-Content -Path $dp -Value "select disk $DiskNumber", "clean", "rescan" -Encoding ASCII
    $dpOut = & diskpart.exe /s $dp 2>&1
    Remove-Item $dp
    if ($LASTEXITCODE -ne 0) { throw "diskpart clean failed: $($dpOut -join ' ')" }
    if (Get-Partition -DiskNumber $DiskNumber -ErrorAction SilentlyContinue) { throw "disk $DiskNumber still has partitions after clean" }

    Add-Type -Namespace Native -Name Disk -MemberDefinition @'
[DllImport("kernel32.dll", SetLastError = true, CharSet = CharSet.Unicode)]
public static extern Microsoft.Win32.SafeHandles.SafeFileHandle CreateFile(
    string name, uint access, uint share, IntPtr sec, uint disposition, uint flags, IntPtr template);
'@
    # GENERIC_READ|GENERIC_WRITE, FILE_SHARE_READ|WRITE, OPEN_EXISTING, FILE_FLAG_WRITE_THROUGH
    # Decimal on purpose: PowerShell parses 0xC0000000 as a negative Int32.
    $h = [Native.Disk]::CreateFile("\\.\PhysicalDrive$DiskNumber", [uint32]3221225472, 3, [IntPtr]::Zero, 3, [uint32]2147483648, [IntPtr]::Zero)
    if ($h.IsInvalid) { throw "CreateFile failed: Win32 error $([Runtime.InteropServices.Marshal]::GetLastWin32Error())" }
    $out = New-Object IO.FileStream($h, [IO.FileAccess]::ReadWrite)

    $psi = New-Object Diagnostics.ProcessStartInfo($Xz, "-dc `"$Image`"")
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.CreateNoWindow = $true
    $proc = [Diagnostics.Process]::Start($psi)
    $in = $proc.StandardOutput.BaseStream

    $chunk = 4MB
    $buf = New-Object byte[] $chunk
    $total = 0L
    $nextLog = 0L
    # The first chunk holds the partition table. Write it last: once Windows
    # sees partitions it mounts them and blocks raw writes into them.
    $first = $null
    Log "writing $Image"
    while ($true) {
        $fill = 0
        while ($fill -lt $chunk) {
            $n = $in.Read($buf, $fill, $chunk - $fill)
            if ($n -le 0) { break }
            $fill += $n
        }
        if ($fill -eq 0) { break }
        $len = [int]([math]::Ceiling($fill / 512) * 512)   # raw disk writes must be sector-aligned
        if ($len -gt $fill) { [Array]::Clear($buf, $fill, $len - $fill) }
        if ($null -eq $first) {
            $first = New-Object byte[] $len
            [Array]::Copy($buf, $first, $len)
            [void]$out.Seek($len, [IO.SeekOrigin]::Begin)
        } else {
            $out.Write($buf, 0, $len)
        }
        $total += $fill
        if ($total -ge $nextLog) { Log ("{0:N0} MB" -f ($total / 1MB)); $nextLog += 256MB }
        if ($fill -lt $chunk) { break }
    }
    $out.Flush()
    [void]$out.Seek(0, [IO.SeekOrigin]::Begin)
    $out.Write($first, 0, $first.Length)
    $out.Flush()
    $out.Close()
    $proc.WaitForExit()
    if ($proc.ExitCode -ne 0) { throw "xz exited with $($proc.ExitCode)" }
    Update-Disk -Number $DiskNumber
    Log ("DONE {0:N0} MB" -f ($total / 1MB))
} catch {
    Log "ERROR $($_.Exception.Message)"
    exit 1
}
