Get-Content C:\cfm\runner.log -Tail 4
Get-ChildItem C:\cfm\queue | ForEach-Object Name
Get-CimInstance Win32_Process | Where-Object { $_.Name -like 'cfm*' -or $_.CommandLine -like '*runner.ps1*' } | ForEach-Object { "$($_.ProcessId) $($_.Name) $($_.CommandLine)" }
