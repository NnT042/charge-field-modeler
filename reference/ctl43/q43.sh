#!/bin/bash
# Doc 15.43 controls: isotropic-emission proton pair; inert-sphere alpha.
cd /c/Users/nnt04/Project/charge-field-modeler/rust
BIN=/c/Users/nnt04/Project/charge-field-modeler/rust/target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O="/c/Users/nnt04/AppData/Local/Temp/claude/1e8dfebd78308e59/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/ctl43"
export CFM_F_RD=6 CFM_D_FREEZE=1 CFM_D_STEPS=7500 CFM_D_LOG=500
run() { "$BIN" report_deuteron --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > "$O/$1.txt" 2>&1; echo "$1 done $(date)" >> "$O/progress.txt"; }
# A: isotropic pair, pole to pole, 3 apart (swept maze gives +591 / +641)
CFM_F_PAIR=pp CFM_D_SEP0=3 CFM_D_MAZE=iso CFM_D_SRCCAL=1 CFM_D_SRCDIR="$O/src_iso" run iso_pp_axial
# D: alpha with inert spheres for neutrons (narrow neutrons gave p +1015 / +1090)
CFM_F_PAIR=pxxp CFM_D_LAYOUT=custom CFM_D_POS="0,0,-1.9;1,0,0;-1,0,0;0,0,1.9" CFM_D_FSELF="proton:-10.4,-14.1,63.8,9.5,14.0,-7.6" run dead_alpha
# C: isotropic pair side by side (swept gives -256); reuse A's lone-body numbers
F=$(grep "^self-force proton" "$O/iso_pp_axial.txt" | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g')
T=$(grep "^self-torque proton" "$O/iso_pp_axial.txt" | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g')
CFM_F_PAIR=pp CFM_F_GEOM=side CFM_D_SEP0=3 CFM_D_MAZE=iso CFM_D_FSELF="proton:$F,$T" CFM_D_SRCDIR="$O/src_iso" run iso_pp_side
echo "ALL DONE $(date)" >> "$O/progress.txt"
