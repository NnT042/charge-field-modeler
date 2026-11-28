#!/bin/bash
# Doc 15.73 (ESTELLE half): long free n-p-n, aim 1 -> edge -1, seed 49344, 1200 tu. ISABELLE runs seed 1 (q73_a).
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=C:/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm72.exe
O=../reference/tc73
P=$O/progress_estelle.txt
echo "START $(date)" > $P
export CFM_C_SEED=49344 CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_C_DT=0.06 CFM_F_RD=6 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_LOG=167 CFM_C_LAMBDA=0.3 CFM_C_H=0.25 CFM_P_MIX=1,1,1,1 CFM_P_INNER=random CFM_P_SIZE=1 CFM_P_MASS=1 CFM_P_PREFILL=1 CFM_P_NORM=masspath CFM_P_CELLS=1
go() { n=$1; shift; ( export "$@"; $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$n.txt 2>&1; echo "$n done $(date)" >> $P ); }
FR="CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_PIPE=pn CFM_D_PIPE_DEG=90 CFM_D_PIPE_LANE=1"
AIM="CFM_D_PIPE_FOCUS=1 CFM_D_PIPE_FOCUS_EDGE=-1"
go q73_b_npn_free_f1e-1_s49344_long CFM_F_PAIR=npn CFM_D_LAYOUT=custom "CFM_D_POS=0,0,-2.05;0,0,0;0,0,2.05" CFM_D_MASS=7e6 CFM_D_STEPS=20000 CFM_D_PIPE=pn CFM_D_PIPE_DEG=90 CFM_D_PIPE_LANE=1 CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=8.75e5 $AIM
echo "ALL DONE $(date)" >> $P
