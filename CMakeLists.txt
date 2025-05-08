cmake_minimum_required(VERSION 3.10)
project(MonodomainSimulation C CXX CUDA)

# === Build options ===
option(USE_CUDA "Enable CUDA support" ON)
option(USE_OPENMP "Enable OpenMP support" ON)

# === Set folders ===
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(INC_DIR "${CMAKE_SOURCE_DIR}/include")
set(EXT_DIR "${CMAKE_SOURCE_DIR}/external/inih")
set(BIN_DIR "${CMAKE_SOURCE_DIR}/bin")

# === Include paths ===
include_directories(${INC_DIR} ${EXT_DIR})

# === Collect sources ===
file(GLOB SRC_C_FILES "${SRC_DIR}/*.c")
file(GLOB SRC_CU_FILES "${SRC_DIR}/*.cu")
file(GLOB EXT_C_FILES "${EXT_DIR}/*.c")

# === Start with only C sources ===
set(SOURCES ${SRC_C_FILES} ${EXT_C_FILES})

# === If CUDA enabled, add .cu files ===
if(USE_CUDA)
    enable_language(CUDA)
    message(STATUS "CUDA support ENABLED")
    list(APPEND SOURCES ${SRC_CU_FILES})
else()
    message(STATUS "CUDA support DISABLED")
endif()

# === Create executable ===
add_executable(monodomain_simulation ${SOURCES})

# === Set output binary directory ===
set_target_properties(monodomain_simulation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
)

# === Link CUDA libraries if needed ===
if(USE_CUDA)
    target_link_libraries(monodomain_simulation cusparse)
endif()

# === Link OpenMP if enabled ===
if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_FOUND)
        target_link_libraries(monodomain_simulation OpenMP::OpenMP_C)
        message(STATUS "OpenMP support ENABLED")
    endif()
else()
    message(STATUS "OpenMP support DISABLED")
endif()

# === pthread always required if OpenMP is used ===
target_link_libraries(monodomain_simulation pthread m)
