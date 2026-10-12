# Start the CFM runner detached (it outlives the SSH session). With -Restart,
# stop a runner that is already going first (only when it is waiting: a job
# in progress is left alone unless -Force is also given).
param([switch]$Restart, [switch]$Force)
$cmd = 'powershell -NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File C:\cfm\runner.ps1'
$old = Get-CimInstance Win32_Process | Where-Object { $_.Name -eq 'powershell.exe' -and $_.CommandLine -like '*-File C:\cfm\runner.ps1*' }
if ($old -and $Restart) {
    $busy = Get-Process cfm -ErrorAction SilentlyContinue
    if ($busy -and -not $Force) { "a job is running; runner left alone"; exit 1 }
    $old | ForEach-Object { Stop-Process -Id $_.ProcessId -Force }
    Start-Sleep -Seconds 1
    $old = $null
}
if ($old) { "runner already going (pid $($old.ProcessId))"; exit 0 }
$r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd }
"runner started (pid $($r.ProcessId), result $($r.ReturnValue))"
