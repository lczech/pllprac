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
CC      = g++
CCFLAGS = -std=c++11 -O1 $(WARN) $(DBG)
LDFLAGS = -lpll-sse3 -lm

# Includes, if needed
# -I/usr/local/include
# -L/usr/local/lib

# --------------------------------
#   Make rules
# --------------------------------

SRCS = $(wildcard src/*.cpp)
OBJS = $(SRCS:.cpp=.o)

all: $(PROG)
	@echo "\n==========   Done    =========="

$(PROG): $(OBJS)
	@echo "\n==========  Linking  =========="
	@echo "Objects: $(OBJS)\n"
	$(CC) $(OBJS) $(LDFLAGS) -o $@

.cpp.o:
	@echo "\n========== Compiling =========="
	@echo "File: $< > $@\n"
#	@echo "Sources: $(SRCS)"
	$(CC) -c $(CCFLAGS) $< -o $@

clean:
	@echo "\n========== Cleaning  =========="
	rm -fv $(OBJS) $(PROG)
