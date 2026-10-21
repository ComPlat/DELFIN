#requires -Version 5.1
[CmdletBinding()]
param([Parameter(Mandatory=$true)][string]$Destination)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$compiler = @(
    (Join-Path $env:SystemRoot 'Microsoft.NET\Framework64\v4.0.30319\csc.exe'),
    (Join-Path $env:SystemRoot 'Microsoft.NET\Framework\v4.0.30319\csc.exe')
) | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
if (-not $compiler) { throw 'Windows .NET Framework C# compiler is missing.' }
$arguments = @('/nologo','/target:winexe','/reference:System.Windows.Forms.dll',('/win32icon:' + (Join-Path $PSScriptRoot 'DELFIN.ico')),('/out:' + $Destination),(Join-Path $PSScriptRoot 'Starter.cs'))
& $compiler @arguments
if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $Destination)) { throw 'Could not build the DELFIN Windows starter.' }
