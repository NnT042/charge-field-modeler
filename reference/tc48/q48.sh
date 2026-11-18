#!/bin/bash
# Doc 15.48 (ESTELLE half): beam-fate tally on the He3 at contact, then the p-n frozen distance
# sweep with trace bodies on seed 49344 (ISABELLE runs the same sweep on seed 777: q48_pn_*_s777).
# Trace bodies: 64 copies, surface emission, through release, beam 0.222, lambda 0.3.
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=../reference/tc48
P=$O/progress.txt
echo "START $(date)" > $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_C_LAMBDA=0.3 CFM_C_H=0.25 CFM_C_DT=0.06 CFM_F_RD=6
export CFM_D_FREEZE=1 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_STEPS=1667 CFM_D_LOG=167
export CFM_D_FSELF="proton:-106.9,-12.0,22.6,1.5,1.8,-4.2;neutron:1.5,5.8,-5.0,-3.3,1.5,-5.8"
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $P; }
CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 run he3_205_b222_tally
for d in 2.05 2.3 2.6 3 4 5; do CFM_F_PAIR=pn CFM_D_SEP0=$d run pn_${d}_s49344; done
echo "ALL DONE $(date)" >> $P
