@echo off
rem minutes=5
C:\cfm\cfm67.exe gpu::dsmc --test-threads=1 > C:\cfm\out\q67_a_gates.txt 2>&1
exit /b %ERRORLEVEL%
