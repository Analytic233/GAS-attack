param(
    [string]$BuildDir = "build-artifact",
    [string]$LogsDir = "artifact-logs"
)

$ErrorActionPreference = "Stop"

if (-not (Test-Path $LogsDir)) {
    New-Item -ItemType Directory -Path $LogsDir | Out-Null
}

Write-Host "[1/4] Configure and build..."
cmake -G "MinGW Makefiles" -S . -B $BuildDir | Out-File "$LogsDir\cmake_configure.log" -Encoding utf8
cmake --build $BuildDir | Out-File "$LogsDir\cmake_build.log" -Encoding utf8

Write-Host "[2/4] Run FiLIP verifier..."
& ".\$BuildDir\FiLIP_verifier.exe" | Out-File "$LogsDir\FiLIP_verifier.log" -Encoding utf8

Write-Host "[3/4] Run toy attack..."
& ".\$BuildDir\toy_attack.exe" | Out-File "$LogsDir\toy_attack.log" -Encoding utf8

Write-Host "[4/4] Summarize logs..."
& ".\scripts\summarize_harness_logs.ps1" -LogsDir $LogsDir | Tee-Object "$LogsDir\summary.txt"

Write-Host "Done. Logs stored in $LogsDir"
