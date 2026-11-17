#!/bin/bash
# Doc 15.45: same trace bodies (64 copies, lambda 0.3) but emission leaves from the ball's own
# surface at c relative to it (MazeSpeed::EjectC) instead of being placed on the envelope
# (EjectSurface), so emitted charge must also get out through the trace.
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=../reference/tc45
export CFM_I_RUNS="ball:swept:0" CFM_C_LAMBDA=0.3 CFM_C_NA=8000 CFM_C_H=0.25 CFM_C_DT=0.06 CFM_C_WARMUP=1667 CFM_C_AVG=1667 CFM_I_SPEED=ejectc CFM_C_COPIES=64
for cl in proton neutron; do
  CFM_I_CLASS=$cl $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_intake > $O/lone_${cl}_c64_ejectc.txt 2>&1
  echo "$cl ejectc done $(date)" >> $O/progress_eject.txt
done
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_SPEED=ejectc CFM_F_RD=6
export CFM_D_FREEZE=1 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_STEPS=1667 CFM_D_LOG=167
CFM_F_PAIR=pn CFM_D_SEP0=3 $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/tr_pn_axis_3_ejectc.txt 2>&1
echo "ALL DONE $(date)" >> $O/progress_eject.txt
