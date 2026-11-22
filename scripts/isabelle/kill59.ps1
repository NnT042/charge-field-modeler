Get-CimInstance Win32_Process | Where-Object { $_.Name -like 'cfm59*.exe' } | ForEach-Object { "$($_.ProcessId) $($_.Name) started $($_.CreationDate)" }
Get-Content C:\cfm\runner.log -Tail 3
