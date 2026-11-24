@echo off
rem minutes=30
set CFM_PAUSE_FILE=C:\cfm\pause\cfm_pause
set CFM_C_SEED=49344
set CFM_D_INTAKE=ball
set CFM_D_COPIES=64
set CFM_D_POLAR=through
set CFM_D_BEAM=0.222
set CFM_C_DT=0.06
set CFM_F_RD=6
set CFM_D_CALWARMT=150
set CFM_D_CALT=100
set CFM_C_LAMBDA=0.3
set CFM_C_H=0.25
set CFM_P_MIX=1,1,1,1
set CFM_P_INNER=random
set CFM_P_SIZE=1
set CFM_P_MASS=1
set CFM_P_PREFILL=1
set CFM_P_NORM=masspath
set CFM_P_CELLS=1
set CFM_D_LONEONLY=1
set CFM_F_PAIR=p
C:\cfm\cfm60.exe --exact gpu::dsmc::shell_gates::report_deuteron --ignored --nocapture > C:\cfm\out\q60_b_p_up.txt 2>&1
exit /b %ERRORLEVEL%
