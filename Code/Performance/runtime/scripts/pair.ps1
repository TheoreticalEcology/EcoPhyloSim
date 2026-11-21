# Paired CPU-time benchmark (Windows): two builds run at the same time, pinned to the two logical CPUs of one
# physical core, at priority BelowNormal; the sides are swapped in every repetition. Both runs then share the
# same core and see the same background load, so the ratio of their CPU times is stable even on a busy machine
# (A/A test: pair ratios 0.999-1.023, see results/pairs-aa.csv). Each repetition uses seed seed0 + i for both runs.
# Example:
#   .\pair.ps1 -tagA master -libA C:\libs\master -tagB proto -libB C:\libs\proto `
#              -setting "dispersal=1,runs=2000" -reps 4 -out ..\results\new.csv
param([string]$tagA, [string]$libA, [string]$tagB, [string]$libB, [string]$setting, [int]$reps, [string]$out,
      [int]$seed0 = 1, [int[]]$cpus = @(14, 15), [string]$Rscript = "C:\Program Files\R\R-4.6.2\bin\x64\Rscript.exe")
$script = Join-Path $PSScriptRoot "run-one.R"
for ($i = 0; $i -lt $reps; $i++) {
  $seed = $seed0 + $i
  if ($i % 2 -eq 0) { $c1 = $cpus[0]; $c2 = $cpus[1] } else { $c1 = $cpus[1]; $c2 = $cpus[0] }
  $pa = Start-Process $Rscript -ArgumentList "`"$script`"", $tagA, "`"$($libA.Replace('\','/'))`"", $setting, $seed, "`"$out`"", $c1 -PassThru -NoNewWindow
  $pb = Start-Process $Rscript -ArgumentList "`"$script`"", $tagB, "`"$($libB.Replace('\','/'))`"", $setting, $seed, "`"$out`"", $c2 -PassThru -NoNewWindow
  foreach ($x in @(@($pa, $c1), @($pb, $c2))) { $x[0].ProcessorAffinity = [IntPtr](1 -shl $x[1]); $x[0].PriorityClass = 'BelowNormal' }
  if (-not $pa.WaitForExit(900000)) { $pa.Kill(); "A killed" }
  if (-not $pb.WaitForExit(900000)) { $pb.Kill(); "B killed" }
}
