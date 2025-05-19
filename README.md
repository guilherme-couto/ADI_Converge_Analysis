# Finite Difference Numerical Methods for Monodomain

This project provides a flexible and modular framework for solving the **monodomain equation** using a variety of finite difference numerical methods (e.g., FE, OS-ADI, SSI-ADI), with support for **CPU**, **OpenMP**, and **GPU (CUDA)** backends.

## Features

- Modular architecture with support for different cell models and numerical methods.
- Parallel execution using OpenMP and CUDA.
- Configurable simulation parameters via `.ini` files (no need to recompile).
- Custom output saving options (frames, final state, etc.).
- Designed for reproducibility and extensibility.

---

## Project Structure

/
├── include/ # Header files
├── src/ # Source code (.c and .cu)
├── external/ # External libraries (e.g., inih)
├── configs/ # Example .ini configuration files
├── bin/ # Compiled executables
├── CMakeLists.txt # CMake build configuration
├── README.md # You're here!
└── build/ # CMake-generated build directory

---

## Requirements

- CMake ≥ 3.10
- GCC or Clang (for CPU/OpenMP)
- NVIDIA CUDA Toolkit (for GPU support) — optional
- OpenMP (optional but recommended for CPU parallelism)

> 💡 You can disable CUDA or OpenMP support at compile time.

---

## Build Instructions

### 1. Clone the repository

```bash
git clone https://github.com/guilherme-couto/numerical_methods_monodomain.git
cd numerical_methods_monodomain
```

### 2. Create a build directory and run CMake
```bash
mkdir build
cd build
cmake ..                # Default: builds with CUDA and OpenMP
```

Custom build options:
```bash
# Disable CUDA
cmake .. -DUSE_CUDA=OFF

# Disable OpenMP
cmake .. -DUSE_OPENMP=OFF

# Disable both
cmake .. -DUSE_CUDA=OFF -DUSE_OPENMP=OFF
```

### 3. Compile
```bash
make
```
The binary will be placed in the `bin/` directory.

---

## Running the Simulation

The executable requires a path to a `.ini` configuration file that defines all simulation parameters.

```bash
./bin/monodomain_simulation configs/example_config.ini
```

---

## Configuration File
Example `configs/example_config.ini`:
```ini
[simulation]
execution_mode=SERIAL
equation_type=Monodomain
cell_model=AFHN
numerical_method=SSIADI

total_time=1000 # ms
length_x=2.0 # cm
length_y=2.0 # cm

dt=0.01 # ms
dx=0.01 # cm
dy=0.01 # cm
sigma=1.2e-3 # Diffusion coefficient

save_function=save_as_vtk
save_rate=100 # ms
output_dir=outputs/afhn_spiral
save_frames=true
measure_velocity=true

[stim_1]
amplitude=100
start_time=0.0 # ms
duration=2.0 # ms
x_min=0.0 # cm
x_max=0.2 # cm
y_min=0.0 # cm
y_max=2.0 # cm

[stim_2]
amplitude=100
start_time=120.0 # ms
duration=2.0 # ms
x_min=0.0 # cm
x_max=1.0 # cm
y_min=0.0 # cm
y_max=1.0 # cm

```

---

## Technical Notes
- CUDA compilation requires `USE_CUDA=ON`
- OpenMP support is automatic if available
- Outputs are saved in the specified output directory

