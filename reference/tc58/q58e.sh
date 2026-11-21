#!/bin/bash
# Doc 15.58 (ESTELLE): flip tests (free + frozen) and the photon/antiphoton balance sweep (phi).
# Standard tc56 config. ISABELLE (C:\cfm\queue, q58_*) runs the sideways-stiffness and tilt-curve screens.
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm56.exe
O=../reference/tc58
P=$O/progress_estelle.txt
echo "START $(date)" > $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_F_RD=6 CFM_F_PAIR=pn CFM_C_SEED=49344
export CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_LOG=167 CFM_C_DT=0.06 CFM_C_LAMBDA=0.3 CFM_C_H=0.25
go() { n=$1; shift; ( export "$@"; run $n ); }
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $P; }
LANE="CFM_D_PIPE=pn CFM_D_PIPE_DEG=90 CFM_D_PIPE_LANE=1"
ROT="CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=8.75e5"
TRI="CFM_D_LAYOUT=custom CFM_D_POS=0,0,-2.05;0,0,0;0,0,2.05"
# 1. Which end loses its proton? Free p-n-p + lane with the middle neutron flipped, and the whole chain flipped.
go pnp_free_lane_flipmid $LANE $ROT $TRI CFM_F_PAIR=pnp CFM_D_AXES="0,0,1;0,0,-1;0,0,1" CFM_D_MASS=7e6 CFM_D_STEPS=10000
go pnp_free_lane_flipall $LANE $ROT $TRI CFM_F_PAIR=pnp CFM_D_AXES="0,0,-1;0,0,-1;0,0,-1" CFM_D_MASS=7e6 CFM_D_STEPS=10000
# 2. Frozen versions of the same (the end-to-end force difference).
go pnp_frozen_lane_flipmid $LANE $TRI CFM_F_PAIR=pnp CFM_D_AXES="0,0,1;0,0,-1;0,0,1" CFM_D_FREEZE=1 CFM_D_STEPS=1667
go pnp_frozen_lane_flipall $LANE $TRI CFM_F_PAIR=pnp CFM_D_AXES="0,0,-1;0,0,-1;0,0,-1" CFM_D_FREEZE=1 CFM_D_STEPS=1667
# 3. Photon/antiphoton balance (CFM_I_PHI = photon fraction of the ambient; 0.667 = Earth, the tc56 baseline).
#    phi 0.5 = balanced (Mathis: proton and neutron indistinguishable), 0.8 and 1.0 = stronger imbalance.
for phi in 0.5 0.8 1.0; do
  go pn_2.05_phi$phi CFM_I_PHI=$phi CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_SEP0=2.05
  go pn_2.05_lane_phi$phi $LANE CFM_I_PHI=$phi CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_SEP0=2.05
done
echo "ALL DONE $(date)" >> $P
