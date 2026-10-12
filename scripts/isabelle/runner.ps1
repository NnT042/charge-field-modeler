# CFM job runner for ISABELLE. Runs C:\cfm\queue\*.cmd in name order, one at
# a time, and stops when the queue is empty. Each job is a .cmd that sets its
# settings and runs C:\cfm\cfm.exe; a line "rem minutes=NN" gives its length.
# Before each job it asks cfm.exe (guard_check) whether the GPU may be used
# for that long: pause file, school hours (Mon-Fri 07:00-17:00), Genshin
# Impact running. A job that gives way (exit 75) stays queued and runs again.
$root = 'C:\cfm'
$mutex = $null
try { $mutex = New-Object System.Threading.Mutex($false, 'Global\cfm_runner') } catch { $mutex = New-Object System.Threading.Mutex($false, 'cfm_runner') }
if (-not $mutex.WaitOne(0)) { exit 0 }   # another runner is already going
function Log($m) { Add-Content "$root\runner.log" ("{0}  {1}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $m) }
$env:CFM_PAUSE_FILE = "$root\pause\cfm_pause"
New-Item -ItemType Directory -Force "$root\queue", "$root\done", "$root\out" | Out-Null
Log "runner started"
$waiting = $null
while ($true) {
    $job = Get-ChildItem "$root\queue\*.cmd" | Sort-Object Name | Select-Object -First 1
    if (-not $job) { break }
    $need = 110
    $m = Select-String -Path $job.FullName -Pattern 'rem minutes=(\d+)' | Select-Object -First 1
    if ($m) { $need = [int]$m.Matches[0].Groups[1].Value }
    $env:CFM_GUARD_NEED_MIN = "$need"
    $chk = & "$root\cfm.exe" guard_check --ignored --nocapture --exact gpu::dsmc::guard::tests::guard_check 2>$null
    $ok = ($LASTEXITCODE -eq 0)
    $env:CFM_GUARD_NEED_MIN = $null
    if (-not $ok) {
        $why = "$(($chk | Select-String '^blocked' | Select-Object -First 1).Line)"
        $kind = ($why -replace '\d+', '#')
        if ($kind -ne $waiting) { Log "waiting ($($job.BaseName)): $why"; $waiting = $kind }
        Start-Sleep -Seconds 60
        continue
    }
    $waiting = $null
    Log "start $($job.BaseName) (about $need min)"
    cmd /c "`"$($job.FullName)`""
    $code = $LASTEXITCODE
    if ($code -eq 75) {
        Log "gave way: $($job.BaseName), will run again"
        Start-Sleep -Seconds 60
        continue
    }
    Move-Item -Force $job.FullName "$root\done\"
    Log "done $($job.BaseName) exit $code"
}
Log "queue empty, runner stopped"
$mutex.ReleaseMutex()
