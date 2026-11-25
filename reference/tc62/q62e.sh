#!/bin/bash
# Doc 15.62 (ESTELLE half): free lone proton in its own frame (field drifts past at -V), masspath mix.
# ISABELLE runs q62_a..c (gates; tilt 30 at V0 0.017 with spin L 8.75e5 and with no spin).
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm62.exe
O=../reference/tc62
P=$O/progress_estelle.txt
echo "START $(date)" > $P
export CFM_C_SEED=49344 CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_C_DT=0.06 CFM_F_RD=6
export CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_C_LAMBDA=0.3 CFM_C_H=0.25
export CFM_P_MIX=1,1,1,1 CFM_P_INNER=random CFM_P_SIZE=1 CFM_P_MASS=1 CFM_P_PREFILL=1 CFM_P_NORM=masspath CFM_D_LONEONLY=1 CFM_F_PAIR=p
export CFM_D_LONEFREE_T=1200 CFM_D_LOG=500
go() { n=$1; shift; ( export "$@"; $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$n.txt 2>&1; echo "$n done $(date)" >> $P ); }
go q62_d_free_tilt30_spin_x10 CFM_D_MASS=7e5 CFM_D_INERTIA=2.8e5 CFM_D_LSPIN=8.75e4 CFM_D_WIND=0,0,-0.017 CFM_D_LONEAXIS=0.5,0,0.8660254
go q62_e_free_rest_up CFM_D_MASS=7e6 CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=8.75e5
echo "ALL DONE $(date)" >> $P
