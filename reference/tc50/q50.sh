#!/bin/bash
# Doc 15.50 (ESTELLE half): is the scattering between the bodies grid or physics, and how does it
# change with field density? Frozen coaxial p-n, trace bodies, through release, beam 0.222, seed 49344.
#  a) grid check: lambda 0.3 with h 0.125 (vs h 0.25 in tc48/49)
#  b) thinner field: lambda 0.6 and 1.0 (h 0.25, dt 0.06)
# ISABELLE runs the denser field, lambda 0.15 (h 0.15, dt 0.03), q50_*.
# The first run of each setting measures the lone self-forces; later runs reuse them.
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=../reference/tc50
P=$O/progress.txt
echo "START $(date)" > $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_F_RD=6 CFM_F_PAIR=pn
export CFM_D_FREEZE=1 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_STEPS=1667 CFM_D_LOG=167 CFM_C_DT=0.06
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $P; }
grab() { grep "^self-$2 $3" $O/$1.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g'; }
fself() { echo "proton:$(grab $1 force proton),$(grab $1 torque proton);neutron:$(grab $1 force neutron),$(grab $1 torque neutron)"; }
sweep() { # tag lambda h dists...
  local tag=$1 lam=$2 hh=$3; shift 3
  export CFM_C_LAMBDA=$lam CFM_C_H=$hh
  local first=$1; shift
  unset CFM_D_FSELF; run pn_${first}_${tag}
  local FS=$(fself pn_${first}_${tag}); echo "FSELF $tag $FS" >> $P
  for d in "$@"; do CFM_D_FSELF="$FS" CFM_D_SEP0=$d run pn_${d}_${tag}; done
}
export CFM_D_SEP0=3; sweep l03_h0125 0.3 0.125 3 2.05
export CFM_D_SEP0=3; sweep l06 0.6 0.25 3 2.05 4
export CFM_D_SEP0=3; sweep l10 1.0 0.25 3 2.05 4
echo "ALL DONE $(date)" >> $P
