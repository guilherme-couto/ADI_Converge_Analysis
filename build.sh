#!/bin/bash

# =============================================================================
# Build script for Monodomain Simulation
# =============================================================================

USE_CUDA=ON
USE_OPENMP=ON
USE_FLOAT=OFF
BUILD_DIR="build"
CMAKE_ARGS=""

# Parse command-line arguments
for arg in "$@"; do
  case $arg in
    --no-cuda)       USE_CUDA=OFF ;;
    --no-openmp)     USE_OPENMP=OFF ;;
    --float)         USE_FLOAT=ON ;;
    --clean)         CLEAN=true ;;
  esac
done

# CUDA architecture detection
detect_cuda_arch() {
  if ! command -v nvidia-smi &> /dev/null; then
    echo "75"
    return
  fi

  compute_cap=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits 2>/dev/null | head -n1)

  if [[ -z "$compute_cap" ]]; then
    echo "75"
    return
  fi

  major=$(echo "$compute_cap" | cut -d'.' -f1)
  minor=$(echo "$compute_cap" | cut -d'.' -f2)
  arch="${major}${minor}"

  echo "$arch"
}

# Clean build directory
if [ "$CLEAN" = true ]; then
  echo "🧹 Cleaning build directory..."
  rm -rf "$BUILD_DIR"
fi

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR" || exit

# Add precision flag
if [ "$USE_FLOAT" = "ON" ]; then
  CMAKE_ARGS+=" -DUSE_FLOAT=ON"
else
  CMAKE_ARGS+=" -DUSE_FLOAT=OFF"
fi

# Add CUDA architecture
if [ "$USE_CUDA" = "ON" ]; then
  ARCH=$(detect_cuda_arch)
  echo "🧠 Detected CUDA architecture: $ARCH"
  CMAKE_ARGS+=" -DCMAKE_CUDA_ARCHITECTURES=$ARCH"
fi

# Run CMake
echo "⚙️  Running CMake..."
cmake .. -DUSE_CUDA=$USE_CUDA -DUSE_OPENMP=$USE_OPENMP $CMAKE_ARGS

# Compile
echo "🔨 Building..."
make -j$(nproc)

echo "✅ Build complete."
echo "👉 Run the simulation with: ./bin/monodomain_simulation path/to/config.ini"
