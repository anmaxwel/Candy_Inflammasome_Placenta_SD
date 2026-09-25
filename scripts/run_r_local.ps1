param(
    [Parameter(Mandatory = $true, Position = 0)]
    [string]$Script,

    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$ScriptArgs
)

$ErrorActionPreference = 'Stop'
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$allowedRoot = 'C:\Users\antho\Desktop\git_repos\Candy_Inflammasome_Placenta_SD'

if ($projectRoot -ne $allowedRoot) {
    throw "R execution is restricted to $allowedRoot; resolved project root was $projectRoot"
}

$localDirs = @(
    (Join-Path $projectRoot '.r-local\library'),
    (Join-Path $projectRoot '.r-user'),
    (Join-Path $projectRoot '.r-cache'),
    (Join-Path $projectRoot '.r-tmp'),
    (Join-Path $projectRoot 'renv\cache')
)
foreach ($dir in $localDirs) {
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
}

$env:R_USER = Join-Path $projectRoot '.r-user'
$env:R_LIBS_USER = Join-Path $projectRoot '.r-local\library'
$env:R_ENVIRON_USER = Join-Path $projectRoot '.Renviron'
$env:R_PROFILE_USER = Join-Path $projectRoot '.Rprofile'
$env:R_CACHE_DIR = Join-Path $projectRoot '.r-cache'
$env:XDG_CACHE_HOME = Join-Path $projectRoot '.r-cache'
$env:TEMP = Join-Path $projectRoot '.r-tmp'
$env:TMP = Join-Path $projectRoot '.r-tmp'
$env:TMPDIR = Join-Path $projectRoot '.r-tmp'
$env:RENV_PATHS_CACHE = Join-Path $projectRoot 'renv\cache'
$env:RENV_PATHS_LIBRARY = Join-Path $projectRoot 'renv\library'
$env:RENV_CONFIG_CACHE_ENABLED = 'FALSE'
$env:LANG = 'English_United States.utf8'
$env:LC_ALL = 'English_United States.utf8'
$env:LC_CTYPE = 'English_United States.utf8'

$rscript = 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe'
if (-not (Test-Path -LiteralPath $rscript)) {
    throw "Expected R runtime not found at $rscript"
}

$scriptPath = if ([IO.Path]::IsPathRooted($Script)) {
    $Script
} else {
    Join-Path $projectRoot $Script
}
$scriptPath = (Resolve-Path -LiteralPath $scriptPath).Path
if (-not $scriptPath.StartsWith($projectRoot, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to run an R script outside the project root: $scriptPath"
}

Push-Location $projectRoot
try {
    & $rscript $scriptPath @ScriptArgs
    exit $LASTEXITCODE
} finally {
    Pop-Location
}
