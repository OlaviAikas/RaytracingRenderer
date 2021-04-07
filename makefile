#Where to look for specific extensions (might not actually be required)
vpath %.h headers
vpath %.cpp source
vpath %.o objects

# Lists of files to compile
SRCS := $(wildcard source/*.cpp)
# NAMES := $(SRCS:source/='')
HEADERS := $(wildcard headers/*.hpp)
OBJS := $(patsubst source/%.cpp,objects/%.o,$(SRCS))

# Dependency files (see §6.4 in http://nuclear.mutantstargoat.com/articles/make/)
# DEP = $(OBJS:.o=.d)  # one dependency file for each source

# Compiler and flags
CXX := g++
CFLAGS := --std=c++11 -O3

# Rules
all: exec

debug: CFLAGS := --std=c++11 -g
debug: clean exec
	
verbose: CFLAGS := -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy \
				   -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs \
				   -Wnoexcept -Woverloaded-virtual -Wredundant-decls \
				   -Wsign-conversion -Wstrict-null-sentinel -Wstrict-overflow=3 \
				   -Wswitch-default -Wundef -Werror -Wno-unused
verbose: clean exec

objects/%.o: %.cpp ${HEADERS}
	$(CXX) $(CFLAGS) -c -o $@ $<

exec: ${OBJS} ${HEADERS}
	$(CXX) $(CFLAGS) -o exec main.cpp ${OBJS}

test: ${OBJS} ${HEADERS}
	$(CXX) $(CFLAGS) -o test test.cpp ${OBJS}

clean:
	rm -f ${OBJS} exec test
