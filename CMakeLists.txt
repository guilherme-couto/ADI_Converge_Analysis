cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0104 NEW)

project(MonodomainSimulation C CXX CUDA)

# === Set Release build and flags ===
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -finline-functions -funroll-loops -flto -fno-fast-math -funsafe-math-optimizations -fassociative-math -freciprocal-math -fno-trapping-math")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -finline-functions -funroll-loops -flto -fno-fast-math -funsafe-math-optimizations -fassociative-math -freciprocal-math -fno-trapping-math")
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -Xcompiler=-march=native -finline-functions -funroll-loops -flto -fno-fast-math -funsafe-math-optimizations -fassociative-math -freciprocal-math -fno-trapping-math")

# === Build options ===
option(USE_CUDA "Enable CUDA support" ON)
option(USE_OPENMP "Enable OpenMP support" ON)
option(USE_FLOAT "Use float precision instead of double" OFF)

# === Directories ===
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(CELL_MODELS_DIR "${SRC_DIR}/cell_models")
set(NUMERICAL_METHODS_DIR "${SRC_DIR}/numerical_methods")
set(INC_DIR "${CMAKE_SOURCE_DIR}/include")
set(EXT_DIR "${CMAKE_SOURCE_DIR}/external/inih")
set(BIN_DIR "${CMAKE_SOURCE_DIR}/bin")

# === Include paths ===
include_directories(${INC_DIR} ${EXT_DIR})
include_directories(${CELL_MODELS_DIR})
file(GLOB_RECURSE CELL_HEADERS "${CELL_MODELS_DIR}/*.h")
file(GLOB_RECURSE CELL_C_FILES "${CELL_MODELS_DIR}/*.c")
file(GLOB_RECURSE CELL_CU_FILES "${CELL_MODELS_DIR}/*.cu")
include_directories(${NUMERICAL_METHODS_DIR})
file(GLOB_RECURSE NUMERICAL_METHODS_HEADERS "${NUMERICAL_METHODS_DIR}/*.h")
file(GLOB_RECURSE NUMERICAL_METHODS_C_FILES "${NUMERICAL_METHODS_DIR}/*.c")
file(GLOB_RECURSE NUMERICAL_METHODS_CU_FILES "${NUMERICAL_METHODS_DIR}/*.cu")

# === Other source files ===
file(GLOB SRC_C_FILES "${SRC_DIR}/*.c")
file(GLOB SRC_CU_FILES "${SRC_DIR}/*.cu")
file(GLOB EXT_C_FILES "${EXT_DIR}/*.c")

set(SOURCES ${SRC_C_FILES} ${EXT_C_FILES} ${CELL_C_FILES} ${NUMERICAL_METHODS_C_FILES})

# === CUDA support ===
if(USE_CUDA)
    enable_language(CUDA)
    message(STATUS "✅ CUDA support ENABLED")
    list(APPEND SOURCES ${SRC_CU_FILES} ${CELL_CU_FILES} ${NUMERICAL_METHODS_CU_FILES})
    set(CUDA_LIBRARIES cusparse)
else()
    message(STATUS "❌ CUDA support DISABLED")
endif()

# === Create executable ===
add_executable(monodomain_simulation ${SOURCES})

# === Output binary ===
set_target_properties(monodomain_simulation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
)

# === CUDA libraries ===
if(USE_CUDA)
    target_link_libraries(monodomain_simulation ${CUDA_LIBRARIES})
endif()

# === OpenMP support ===
if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_FOUND)
        target_link_libraries(monodomain_simulation OpenMP::OpenMP_C)
        message(STATUS "✅ OpenMP support ENABLED")
    endif()
else()
    message(STATUS "❌ OpenMP support DISABLED")
endif()

# === Precision ===
if(USE_FLOAT)
    target_compile_definitions(monodomain_simulation PRIVATE USE_FLOAT)
    message(STATUS "🧮 Precision set to: float")
else()
    target_compile_definitions(monodomain_simulation PRIVATE USE_DOUBLE)
    message(STATUS "🧮 Precision set to: double")
endif()

# === System libraries ===
target_link_libraries(monodomain_simulation pthread m)

# === CUDA architecture (logging only) ===
if(DEFINED CMAKE_CUDA_ARCHITECTURES)
    message(STATUS "🧠 CUDA architecture: ${CMAKE_CUDA_ARCHITECTURES}")
endif()
