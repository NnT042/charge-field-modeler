#!/bin/bash
# Doc 15.54 (ESTELLE): quantum channel exchange (photon5.pdf), straight pipe, p and n piped at 90 deg, lambda 0.3, h 0.25,
# seed 49344. ISABELLE runs the held helix (levels 8 and 10).
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm54.exe
O=../reference/tc54
P=$O/progress.txt
echo "START $(date)" > $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_F_RD=6 CFM_F_PAIR=pn
export CFM_D_FREEZE=1 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_STEPS=1667 CFM_D_LOG=167 CFM_C_DT=0.06
export CFM_C_LAMBDA=0.3 CFM_C_H=0.25
export CFM_D_PIPE=pn CFM_D_PIPE_DEG=90 CFM_D_PIPE_SPIN=1 CFM_D_PIPE_ETA=0 CFM_D_PIPE_QUANTUM=1
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $P; }
grab() { grep "^self-$2 $3" $O/$1.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g'; }
fself() { echo "proton:$(grab $1 force proton),$(grab $1 torque proton);neutron:$(grab $1 force neutron),$(grab $1 torque neutron)"; }
unset CFM_D_FSELF; CFM_D_SEP0=2.05 run pn_2.05_quantum
FS=$(fself pn_2.05_quantum); echo "FSELF $FS" >> $P
CFM_D_FSELF="$FS" CFM_D_SEP0=3 run pn_3_quantum
echo "ALL DONE $(date)" >> $P
