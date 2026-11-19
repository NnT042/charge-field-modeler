#!/bin/bash
# Doc 15.53 (ESTELLE half): lane distance sweep. Proton and neutron both piped at 90 deg with the
# protected lane (the +55/+43 attractive setting of 15.52), lambda 0.3, h 0.25, seed 49344.
# Finds where the contact pull turns into the d-3 push (-393 p-only). ISABELLE: 2.35, 3, and 2.05 seed 2.
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm_q53.exe
O=../reference/tc53
P=$O/progress.txt
echo "START $(date)" > $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_F_RD=6 CFM_F_PAIR=pn
export CFM_D_FREEZE=1 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_STEPS=1667 CFM_D_LOG=167 CFM_C_DT=0.06
export CFM_C_LAMBDA=0.3 CFM_C_H=0.25
export CFM_D_PIPE=pn CFM_D_PIPE_DEG=90 CFM_D_PIPE_LANE=1
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $P; }
grab() { grep "^self-$2 $3" $O/$1.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g'; }
fself() { echo "proton:$(grab $1 force proton),$(grab $1 torque proton);neutron:$(grab $1 force neutron),$(grab $1 torque neutron)"; }
unset CFM_D_FSELF; CFM_D_SEP0=2.2 run pn_2.2_lane90pn
FS=$(fself pn_2.2_lane90pn); echo "FSELF $FS" >> $P
for d in 2.5 2.75; do CFM_D_FSELF="$FS" CFM_D_SEP0=$d run pn_${d}_lane90pn; done
echo "ALL DONE $(date)" >> $P
