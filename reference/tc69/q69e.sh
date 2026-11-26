#!/bin/bash
# Doc 15.69 (ESTELLE half): free p-n seed 1, p-n-p frozen + lane, p-n-p free + lane, p-n frozen + lane seed 1. ISABELLE runs q69_a..d.
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm68.exe
O=../reference/tc69
P=$O/progress_estelle.txt
echo "START $(date)" > $P
export CFM_C_SEED=49344 CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_C_DT=0.06 CFM_F_RD=6 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_LOG=167 CFM_C_LAMBDA=0.3 CFM_C_H=0.25 CFM_P_MIX=1,1,1,1 CFM_P_INNER=random CFM_P_SIZE=1 CFM_P_MASS=1 CFM_P_PREFILL=1 CFM_P_NORM=masspath CFM_P_CELLS=1
go() { n=$1; shift; ( export "$@"; $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$n.txt 2>&1; echo "$n done $(date)" >> $P ); }
go q69_e_pn_free_s1 CFM_C_SEED=1 CFM_F_PAIR=pn CFM_D_LAYOUT=custom "CFM_D_POS=0,0,-1.025;0,0,1.025" CFM_D_MASS=7e6 CFM_D_STEPS=10000 CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=8.75e5
go q69_f_pnp_frozen_lane CFM_F_PAIR=pnp CFM_D_LAYOUT=custom "CFM_D_POS=0,0,-2.05;0,0,0;0,0,2.05" CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_PIPE=pn CFM_D_PIPE_DEG=90 CFM_D_PIPE_LANE=1
go q69_g_pnp_free_lane CFM_F_PAIR=pnp CFM_D_LAYOUT=custom "CFM_D_POS=0,0,-2.05;0,0,0;0,0,2.05" CFM_D_MASS=7e6 CFM_D_STEPS=10000 CFM_D_PIPE=pn CFM_D_PIPE_DEG=90 CFM_D_PIPE_LANE=1 CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=8.75e5
go q69_h_pn_frozen_lane_s1 CFM_C_SEED=1 CFM_F_PAIR=pn CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_SEP0=2.05 CFM_D_PIPE=pn CFM_D_PIPE_DEG=90 CFM_D_PIPE_LANE=1
echo "ALL DONE $(date)" >> $P
