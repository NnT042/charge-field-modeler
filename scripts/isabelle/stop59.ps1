Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -like '*-File C:\cfm\runner.ps1*' } | ForEach-Object { "stop runner $($_.ProcessId)"; Stop-Process -Id $_.ProcessId -Force }
Get-Process cfm* -ErrorAction SilentlyContinue | ForEach-Object { "stop $($_.Name) $($_.Id)"; Stop-Process -Id $_.Id -Force }
Start-Sleep 2
Get-Process cfm* -ErrorAction SilentlyContinue
