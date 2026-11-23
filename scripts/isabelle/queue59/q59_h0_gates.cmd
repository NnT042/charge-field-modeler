@echo off
rem minutes=20
set CFM_PAUSE_FILE=C:\cfm\pause\cfm_pause
C:\cfm\cfm59h.exe gpu::dsmc --test-threads=1 > C:\cfm\out\q59_h0_gates.txt 2>&1
exit /b %ERRORLEVEL%
