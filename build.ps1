$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$PythonBin = Get-Command python3 -ErrorAction SilentlyContinue

if (-not $PythonBin) {
    $PythonBin = Get-Command python -ErrorAction SilentlyContinue
}

if (-not $PythonBin) {
    Write-Error "[build] Python 3 is required but was not found on PATH."
    exit 1
}

& $PythonBin.Source (Join-Path $ScriptDir "scripts/build.py") @args
exit $LASTEXITCODE
