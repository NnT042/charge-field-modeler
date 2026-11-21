#!/bin/bash
# Vast.ai box setup for CFM (Ubuntu CUDA template, NVIDIA_DRIVER_CAPABILITIES=all).
# Run once on the box:  bash bootstrap.sh
# Installs Vulkan + Rust, checks the GPU is seen by Vulkan (not llvmpipe), clones, builds,
# runs the 23 gpu::dsmc gates and prints the storage-binding limit (sets the photon cap).
set -e
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq libvulkan1 vulkan-tools git build-essential pkg-config curl tmux > /dev/null

# The NVIDIA container toolkit mounts the driver libraries but some images lack the ICD file.
if ! vulkaninfo --summary 2>/dev/null | grep -qi nvidia; then
  echo "no NVIDIA Vulkan device yet; writing the ICD file"
  mkdir -p /etc/vulkan/icd.d
  cat > /etc/vulkan/icd.d/nvidia_icd.json <<'J'
{ "file_format_version" : "1.0.0",
  "ICD": { "library_path": "libGLX_nvidia.so.0", "api_version" : "1.3" } }
J
fi
vulkaninfo --summary | grep -iE "deviceName|driverName|apiVersion" || true
if ! vulkaninfo --summary 2>/dev/null | grep -qi nvidia; then
  echo "STOP: Vulkan still does not see the NVIDIA GPU (check NVIDIA_DRIVER_CAPABILITIES=all in the template)"; exit 1
fi
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv

if ! command -v cargo > /dev/null; then
  curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal > /dev/null
fi
source "$HOME/.cargo/env"
cd "$HOME"
[ -d cfm ] || git clone -q https://github.com/NnT042/charge-field-modeler.git cfm
cd cfm && git pull -q
cargo build --manifest-path rust/Cargo.toml --profile fast --features gpu --tests 2>&1 | tail -3
# New cargo puts the test binary under target/fast/build/.../out/; take the path printed for src/lib.rs.
BIN=$(cargo test --manifest-path rust/Cargo.toml --profile fast --features gpu --no-run 2>&1 | grep 'src/lib.rs' | sed 's/.*(\(.*\))/\1/')
[ -n "$BIN" ] || { echo "STOP: test binary not found"; exit 1; }
echo "test binary: $BIN"
rm -f "$HOME/cfm_test"; ln -s "$HOME/cfm/$BIN" "$HOME/cfm_test"
cd rust
# One test per process: in a Vast container repeated GpuContext::new in one process fails.
for t in $("$HOME/cfm_test" gpu::dsmc --list 2>/dev/null | grep ': test' | grep -v report_ | sed 's/: test//'); do
  if "$HOME/cfm_test" --exact "$t" 2>&1 | grep -q 'test result: ok'; then echo "ok"; else echo "FAIL $t"; fi
done | sort | uniq -c
"$HOME/cfm_test" --ignored --nocapture --exact gpu::tests::bench_particle_stream 2>&1 | grep -E "adapter:|limits:" || true
echo "BOOTSTRAP OK"
