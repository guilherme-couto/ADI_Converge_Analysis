# ==== Project directories ====
INCDIR   := include
SRCDIR   := src
BINDIR   := bin
EXTDIR   := external

# ==== Compiler and tools ====
CC       := gcc
NVCC     := nvcc

# ==== Default build settings ====
USE_CUDA   ?= 1  # Set to 0 to disable CUDA compilation
USE_OPENMP ?= 1  # Set to 0 to disable OpenMP support

# ==== Compilation flags ====
CFLAGS  := -O3 -I$(INCDIR) -I$(EXTDIR)
LDFLAGS := -lm

ifeq ($(USE_OPENMP),1)
    CFLAGS  += -fopenmp
    LDFLAGS += -fopenmp -lpthread
endif

ifeq ($(USE_CUDA),1)
    CUDA_LDFLAGS := -lcusparse
endif

# ==== Source and object files ====
C_SRCS     := $(wildcard $(SRCDIR)/*.c)
CU_SRCS    := $(wildcard $(SRCDIR)/*.cu)
C_OBJS     := $(C_SRCS:.c=.o)
CU_OBJS    := $(CU_SRCS:.cu=.o)

# Selectively include CUDA object files
ifeq ($(USE_CUDA),1)
    FINAL_OBJS := $(C_OBJS) $(CU_OBJS)
else
    FINAL_OBJS := $(C_OBJS)
endif

# ==== Target binary ====
TARGET := monodomain_simulation

# ==== Build rules ====
.PHONY: all clean

all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(FINAL_OBJS)
	@mkdir -p $(BINDIR)
	$(NVCC) $^ -o $@ $(LDFLAGS) $(CUDA_LDFLAGS)

# Compile .c files with gcc
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

# Compile .cu files with nvcc (only if CUDA is enabled)
%.o: %.cu
ifeq ($(USE_CUDA),1)
	$(NVCC) -c $< -o $@ -I$(INCDIR) -I$(EXTDIR)
endif

clean:
	rm -f $(SRCDIR)/*.o $(BINDIR)/$(TARGET)

# ==== Usage Instructions ====
# To build with both CUDA and OpenMP:
#     make
#
# To disable CUDA:
#     make USE_CUDA=0
#
# To disable OpenMP:
#     make USE_OPENMP=0
#
# To disable both CUDA and OpenMP:
#     make USE_CUDA=0 USE_OPENMP=0
#
# To run the executable (after build):
#     ./bin/monodomain_simulation path/to/config.ini