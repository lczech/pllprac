# --------------------------------
#   Program
# --------------------------------

PROG = pllprac

# --------------------------------
#   Compiler Options
# --------------------------------

# Debug & Profiling (comment out if not needed)
DBG     = -g -pg

# Warning flags
WARN    = -Wall -Wextra -pedantic-errors

# Compiler flags
STDCC   = g++
MPICC   = mpic++
CCFLAGS = -std=c++11 -O1 $(WARN) $(DBG)
LDFLAGS = -lpll-sse3 -lm

# Includes, if needed
# -I/usr/local/include
# -L/usr/local/lib

# --------------------------------
#   Make rules
# --------------------------------

# Collect all files
SRCS = $(wildcard src/*.cpp)
OBJS = $(SRCS:.cpp=.o)

# Target for non-mpi version
all: CC = ${STDCC}
all: $(PROG)
	@echo "\n========== Done std  =========="

# Target for mpi version
mpi: CC = ${MPICC}
mpi: CCFLAGS += -DUSE_MPI
mpi: $(PROG)
	@echo "\n========== Done mpi  =========="

# Target for Linking
$(PROG): $(OBJS)
	@echo "\n==========  Linking  =========="
	@echo "Objects: $(OBJS)\n"
	$(CC) $(OBJS) $(LDFLAGS) -o $@

# Target for Compiling
.cpp.o:
	@echo "\n========== Compiling =========="
	@echo "File: $< > $@"
#	@echo "Sources: $(SRCS)"
	$(CC) -c $(CCFLAGS) $< -o $@

# Target for Cleaning
clean:
	@echo "\n========== Cleaning  =========="
	rm -fv $(OBJS) $(PROG)
