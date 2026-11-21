# Sampling profile (Windows, Very Sleepy 0.91) of one setting of profile-run.R. The library should be built with
# R CMD INSTALL --debug (keeps -O2, adds DWARF symbols and does not strip the DLL). Sampling starts about 1 s after the
# start of the simulation (profile-run.R waits 4 s after writing its pid) and lasts $seconds; the R process is then stopped. Unzip the .sleepy file and analyse it with
#   Rscript profile-analyse.R <unzipped folder> <library>/EcoPhyloSim/libs/x64/EcoPhyloSim.dll
param([string]$setting, [string]$lib, [string]$out, [int]$gens = 1000000, [int]$seconds = 20,
      [string]$Rscript = "C:\Program Files\R\R-4.6.2\bin\x64\Rscript.exe", [string]$sleepy = "C:\Program Files\Very Sleepy\sleepy.exe")
New-Item -ItemType Directory -Force $out | Out-Null
$pidfile = "$out\$setting.pid"; Remove-Item $pidfile -ErrorAction SilentlyContinue
$r = Start-Process $Rscript -ArgumentList "`"$(Join-Path $PSScriptRoot 'profile-run.R')`"", $setting, "`"$pidfile`"", "`"$($lib.Replace('\','/'))`"", $gens `
     -PassThru -NoNewWindow -RedirectStandardOutput "$out\$setting.out" -RedirectStandardError "$out\$setting.err"
while (-not (Test-Path $pidfile) -and -not $r.HasExited) { Start-Sleep -Milliseconds 200 }
Start-Sleep -Milliseconds 300
$rpid = [int](Get-Content $pidfile).Trim()   # the R process itself (Rscript.exe starts R as a child process)
(Get-Process -Id $rpid).PriorityClass = 'BelowNormal'
Start-Sleep -Seconds 5
$sl = Start-Process $sleepy -ArgumentList "/a:$rpid", "/t:$seconds", "/o:$out\$setting.sleepy", "/q" -PassThru
$sl.WaitForExit(($seconds + 300) * 1000) | Out-Null
Stop-Process -Id $rpid -ErrorAction SilentlyContinue
