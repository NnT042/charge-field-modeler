#!/bin/bash
# Doc 15.51 (ESTELLE half): pipe test. Frozen coaxial p-n, trace bodies, through release, beam 0.222,
# lambda 0.3, h 0.25, dt 0.06, random pairing, seed 49344 (same as tc48 baselines -833 at d 3, -1758 at 2.05).
# CFM_D_PIPE turns photons entering a body's polar mouth (0.222, the clear hole) within the pipe angle
# onto its axis, so they run down the hole instead of hitting the walls and being captured.
#  a) pipe on protons, 20 deg: d 3, 2.05, 4
#  b) pipe on both, 20 deg: d 3, 2.05
#  c) pipe on protons, 45 deg: d 3, 2.05
# ISABELLE: pipe on protons 90 deg (whole mouth), and 20 deg with local pairing (q51_*).
# The first run of each setting measures the lone self-forces (with the pipe on); later runs reuse them.
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=../reference/tc51
P=$O/progress.txt
echo "START $(date)" > $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_F_RD=6 CFM_F_PAIR=pn
export CFM_D_FREEZE=1 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_STEPS=1667 CFM_D_LOG=167 CFM_C_DT=0.06
export CFM_C_LAMBDA=0.3 CFM_C_H=0.25
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $P; }
grab() { grep "^self-$2 $3" $O/$1.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g'; }
fself() { echo "proton:$(grab $1 force proton),$(grab $1 torque proton);neutron:$(grab $1 force neutron),$(grab $1 torque neutron)"; }
sweep() { # tag pipe deg dists...
  local tag=$1; export CFM_D_PIPE=$2 CFM_D_PIPE_DEG=$3; shift 3
  local first=$1; shift
  unset CFM_D_FSELF; CFM_D_SEP0=$first run pn_${first}_${tag}
  local FS=$(fself pn_${first}_${tag}); echo "FSELF $tag $FS" >> $P
  for d in "$@"; do CFM_D_FSELF="$FS" CFM_D_SEP0=$d run pn_${d}_${tag}; done
}
sweep p20 p 20 3 2.05 4
sweep pn20 pn 20 3 2.05
sweep p45 p 45 3 2.05
echo "ALL DONE $(date)" >> $P
