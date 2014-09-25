# --------------------------------
#   Program
# --------------------------------

PROG = pllprac

# --------------------------------
#   Compiler Options
# --------------------------------

# Use profiling?
#PROF    = -g -pg

# Which warnings?
WARN    = -Wall -Wextra -pedantic-errors

# Compiler flags
CC      = g++
CCFLAGS = -I/usr/local/include -O3 $(WARN) $(PROF)
LDFLAGS = -I/usr/local/include -L/usr/local/lib/ -lpll-sse3 -lm

# --------------------------------
#   Make rules
# --------------------------------

SRCS =$(wildcard src/*.cpp)
OBJS = $(SRCS:.cpp=.o)

all: $(SRCS) $(PROG)
	@echo "\n========== All       =========="

$(PROG): $(OBJS)
	@echo "\n========== Linking   =========="
	@echo "Objects: $(OBJS)\n"
	$(CC) $(OBJS) $(LDFLAGS) -o $@

.cpp.o:
	@echo "\n========== Compiling =========="
	@echo "Sources: $(SRCS)\n"
	$(CC) -c $(CCFLAGS) $< -o $@

clean:
	@echo "\n========== Cleaning  =========="
	rm -fv $(OBJS) $(PROG)
