param(
    [string]$LogsDir = "artifact-logs"
)

$ErrorActionPreference = "Stop"

function Parse-ToyLog {
    param([string]$Path)
    if (-not (Test-Path $Path)) { return "toy_attack: log not found" }
    $content = Get-Content $Path -Raw
    if ($content -match "\[SUCCESS\]") {
        return "toy_attack: SUCCESS"
    }
    if ($content -match "\[FAIL\]") {
        return "toy_attack: FAIL"
    }
    return "toy_attack: UNKNOWN"
}

function Parse-FiLIPLog {
    param([string]$Path)
    if (-not (Test-Path $Path)) { return "FiLIP_verifier: log not found" }
    $lines = Get-Content $Path
    $good = ($lines | Select-String -Pattern "Good Group" | Select-Object -Last 1).Line
    $bad = ($lines | Select-String -Pattern "Bad Group" | Select-Object -Last 1).Line
    if (-not $good -and -not $bad) {
        return "FiLIP_verifier: summary lines not found"
    }
    return "FiLIP_verifier:`n  $good`n  $bad"
}

Write-Output "=== Artifact Harness Summary ==="
Write-Output (Parse-FiLIPLog -Path (Join-Path $LogsDir "FiLIP_verifier.log"))
Write-Output (Parse-ToyLog -Path (Join-Path $LogsDir "toy_attack.log"))
