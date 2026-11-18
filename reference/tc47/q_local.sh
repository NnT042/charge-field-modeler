#!/bin/bash
# Doc 15.47: runs after tc46/q_beam.sh. Local (nearest-neighbour) photon pairing so dense directed
# emission can shadow what is behind it at sub-cell scale (user, 2026-11-17). Same settings as tc46,
# beam 0.222. Step 0 builds and runs the default-off regression gates first.
export PATH=/usr/bin:/bin:$PATH
ROOT=/c/Users/nnt04/Project/charge-field-modeler
O=$ROOT/reference/tc47
P=$O/progress.txt
echo "WAITING $(date)" > $P
until grep -q "ALL DONE" $ROOT/reference/tc46/progress.txt 2>/dev/null; do sleep 60; done
cd $ROOT
echo "BUILD $(date)" >> $P
cargo build --manifest-path rust/Cargo.toml --profile fast --features gpu --tests > $O/build.txt 2>&1 || { echo "BUILD FAILED $(date)" >> $P; exit 1; }
cd rust
B=target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
# 0. Regression: the non-ignored GPU DSMC gates with local pairing OFF (default).
$B gpu::dsmc > $O/gates_default.txt 2>&1; echo "gates done (exit $?) $(date)" >> $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_C_LAMBDA=0.3 CFM_C_H=0.25 CFM_C_DT=0.06 CFM_F_RD=6 CFM_D_BEAM=0.222
export CFM_D_FREEZE=1 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_STEPS=1667 CFM_D_LOG=167
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $P; }
fself() {
  F=$(grep "^self-force neutron" $O/$1.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g')
  T=$(grep "^self-torque neutron" $O/$1.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g')
  echo "proton:-106.9,-12.0,22.6,1.5,1.8,-4.2;neutron:$F,$T"
}
# 1. Deuteron at 3 with local pairing (compare tc46 pn_axis_3_b222: by-source rows + collisions).
CFM_D_LOCALPAIR=1 CFM_F_PAIR=pn CFM_D_SEP0=3 run pn_axis_3_b222_local
FS=$(fself pn_axis_3_b222_local); echo "FSELF local $FS" >> $P
# 2. He3 at contact with local pairing, two seeds (compare tc46 he3_205_b222_s0 / _s777).
CFM_D_LOCALPAIR=1 CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 CFM_D_FSELF="$FS" run he3_205_b222_local_s0
CFM_C_SEED=777 CFM_D_LOCALPAIR=1 CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 CFM_D_FSELF="$FS" run he3_205_b222_local_s777
echo "ALL DONE $(date)" >> $P
