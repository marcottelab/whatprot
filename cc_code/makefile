################################################################################
# Author: Matthew Beauregard Smith                                             #
# Affiliation: The University of Texas at Austin                               #
# Department: Oden Institute and Institute for Cellular and Molecular Biology  #
# PI: Edward Marcotte                                                          #
# Project: Protein Fluorosequencing                                            #
################################################################################

# I have tried to comment this makefile with useful information about what it is
# doing. However, to understand this makefile, the comments will not be enough
# if you do not know the meaning of the following terms.
#   * rule
#   * target
#   * prerequisite
#   * command
#   * recipe
#
# Also, knowledge of the following makefile topics will help with understanding.
#   * macros
#   * functions
#   * pattern rules
#   * automatic variables
#
# Explaining the above terms is beyond the scope of these comments.

# This sets the shell environment. Every makefile should contain this line.
SHELL = /bin/sh

# Here we explicity set the list of suffixes with implicit rules. The first line
# clears the list, the second sets the suffixes.
.SUFFIXES:
.SUFFIXES: .cc .o .d

# Here we set the C++ compiler we plan to use. Other options have been set
# assuming that g++ (the GNU C++ compiler) is being used. If you change this,
# expect things to break.
CC = g++

# Set any flags for libraries.
LDFLAGS = -L./extern/mrmpi-7Apr14/src/ -lmrmpi_serial
LDFLAGS += -L./extern/mrmpi-7Apr14/mpistubs/ -lmpi_stubs
LDFLAGS += -fopenmp
LDFLAGS += -L./extern/levmar-2.6/ -llevmar

LEVMAR = extern/levmar-2.6/liblevmar.a

SRC_DIR = ./src
TEST_DIR = ./test
MAIN = $(SRC_DIR)/main/main.cc
TEST_MAIN = $(TEST_DIR)/test-module.cc
SRCS = $(shell find $(SRC_DIR) $(TEST_DIR) -name '*.cc')
HEADERS = $(shell find $(SRC_DIR) $(TEST_DIR) -name '*.h')

BUILD_DIR = ./build
BUILD_DEBUG_DIR = ./build/debug
BUILD_RELEASE_DIR = ./build/release

DEBUG_OBJS = $(SRCS:./%.cc=$(BUILD_DEBUG_DIR)/%.o)
DEBUG_MAIN_OBJ = $(MAIN:./%.cc=$(BUILD_DEBUG_DIR)/%.o)
DEBUG_TEST_MAIN_OBJ = $(TEST_MAIN:./%.cc=$(BUILD_DEBUG_DIR)/%.o)
DEBUG_DEV_OBJS = $(filter-out %.test.o $(DEBUG_TEST_MAIN_OBJ), $(DEBUG_OBJS))
DEBUG_TEST_OBJS = $(filter-out $(DEBUG_MAIN_OBJ), $(DEBUG_OBJS))

RELEASE_OBJS = $(SRCS:./%.cc=$(BUILD_RELEASE_DIR)/%.o)
RELEASE_MAIN_OBJ = $(MAIN:./%.cc=$(BUILD_RELEASE_DIR)/%.o)
RELEASE_TEST_MAIN_OBJ = $(TEST_MAIN:./%.cc=$(BUILD_RELEASE_DIR)/%.o)
RELEASE_DEV_OBJS = $(filter-out %.test.o $(RELEASE_TEST_MAIN_OBJ), $(RELEASE_OBJS))

OBJS = $(DEBUG_OBJS) $(RELEASE_OBJS)

# Set the output directory.
BIN_DIR = ./bin

# Set the include paths to search for headers.
#   * -iquote is for headers included with quotes
#   * -I is for headers included with angle brackets.
INCLUDE = -iquote$(SRC_DIR)
INCLUDE += -I./extern/mrmpi-7Apr14/src/
INCLUDE += -I./extern/mrmpi-7Apr14/mpistubs/
INCLUDE += -I./extern/boost_1_73_0_subset/
INCLUDE += -I./extern/FakeIt/single_header/boost/
INCLUDE += -I./extern/cxxopts-3.0.0/include/
INCLUDE += -I./extern/cereal-1.3.2/include/
INCLUDE += -I./extern/nlohmann-json-3.11.2/
INCLUDE += -I./extern/levmar-2.6/

# These flags tell g++ to generate dependencies in the form of makefile rules
# as a side effect of compilation. This macro assumes that the target is the
# compilation destination.
#   * '-MMD' does a lot. It combines '-M', '-MM', and '-MD'. It also disables
#     '-E'.
#       * Generate dependencies of the target as makefile rules.
#       * Exclude any dependencies on standard library headers.
#       * Output makefile results into a file. In particular, it takes the file
#         given to the '-o' flag (the compilation destination), removes the last
#         file extension, and puts a '.d' extension in its place.
#       * Compilation continues as normal after generating the makefile info.
#         This is not the default behavior when using most '-M' flags. This is
#         intended to cut down on compilation time.
#   * '-MP' tells g++ to generated so called "phony" targets for header files.
#     This is a weird trick that will stop make from complaining if the header
#     no longer exists because it has been moved or deleted. these targets will
#     have no listed dependencies and no recipe.
#   * '-MT' takes a parameter to set the name of the target in the generated
#     rule. Without this, the compiler will strip off all of the directories
#     from the target name and use the result. '-MT' overrides this behavior.
#     Here we give '$@' so that the generated target has the same name as the
#     target we are currently running.
GENDEPS = -MMD -MP -MT '$@'

# Here we set any additional flags for the compiler.
CFLAGS = -Wall -Wextra -Wno-cast-function-type -Wno-deprecated-copy -Wno-maybe-uninitialized -Wno-unused-parameter -Wpedantic -Werror -fopenmp

# Generate binary. This is first so that it will be the default rule.
all: CFLAGS += -Og -ggdb
all: $(DEBUG_DEV_OBJS) $(LEVMAR)
	@mkdir -p $(BIN_DIR)/debug
	$(CC) $^ $(LDFLAGS) -o $(BIN_DIR)/debug/whatprot

# Generate release version.
release: CFLAGS += -O2 -flto
release: LDFLAGS += $(CFLAGS)
release: $(RELEASE_DEV_OBJS) $(LEVMAR)
	@mkdir -p $(BIN_DIR)/release
	$(CC) $^ $(LDFLAGS) -o $(BIN_DIR)/release/whatprot

# Compile and run tests.
test: CFLAGS += -Og -ggdb
test: $(DEBUG_TEST_OBJS) $(LEVMAR)
	@mkdir -p $(BIN_DIR)/debug
	$(CC) $^ $(LDFLAGS) -o $(BIN_DIR)/debug/test
	$(BIN_DIR)/debug/test

# Cleanup
clean:
	rm -rf $(BIN_DIR)/*
	rm -rf $(BUILD_DIR)/*
	$(MAKE) -C extern/levmar-2.6 cleanall

# Run clang-format
format: $(BUILD_DIR)/format

$(BUILD_DIR)/format: $(SRCS) $(HEADERS)
	clang-format -style=file -i $?
	touch $(BUILD_DIR)/format

# "-include" includes other makefiles.
#   * Each file included is a makefile with one rule for making a .o file with
#     an empty recipe (no commands).
#   * Doing this adds those dependencies to that .o file, without interfering
#     with any recipes we give later to the same target.
#   * Dependencies are generated as a side effect of compilation - see the
#     $(GENDEPS) macro.
#   * Dependencies are generated at the same time the file is compiled, not
#     before. This is ok, The dependencies can't change without also changing
#     the files that were dependencies the most recent time the target was
#     rebuilt.
-include $(OBJS:.o=.d)

# Generate .o and .d files.
$(BUILD_DEBUG_DIR)/%.o: ./%.cc
	@mkdir -p $(shell dirname $@)
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@ $(GENDEPS)

$(BUILD_RELEASE_DIR)/%.o: ./%.cc
	@mkdir -p $(shell dirname $@)
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@ $(GENDEPS)

# Running make in external directories.
$(LEVMAR):
	$(MAKE) -C extern/levmar-2.6 liblevmar.a
