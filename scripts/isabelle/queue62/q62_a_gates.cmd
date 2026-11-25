@echo off
rem minutes=20
set CFM_PAUSE_FILE=C:\cfm\pause\cfm_pause
C:\cfm\cfm62.exe gpu::dsmc --test-threads=1 > C:\cfm\out\q62_a_gates.txt 2>&1
exit /b %ERRORLEVEL%
