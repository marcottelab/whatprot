################################################################################
# Author: Matthew Beauregard Smith                                             #
# Affiliation: The University of Texas at Austin                               #
# Department: Oden Institute and Institute for Cellular and Molecular Biology  #
# PI: Edward Marcotte                                                          #
# Project: Protein Fluorosequencing                                            #
################################################################################

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

SRC_DIR = ./src
TEST_DIR = ./test
MAIN = $(SRC_DIR)/main/main.cc
TEST_MAIN = $(TEST_DIR)/test-module.cc
SRCS = $(shell find $(SRC_DIR) $(TEST_DIR) -name '*.cc')

BUILD_DIR = ./build

OBJS = $(SRCS:./%.cc=$(BUILD_DIR)/%.o)
DEV_OBJS = $(filter-out %.test.cc $(MAIN:./%.cc=$(BUILD_DIR)/%.o), $(OBJS))
TEST_OBJS = $(filter-out $(TEST_MAIN:./%.cc=$(BUILD_DIR)/%.o), $(OBJS))

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

# Here we set any additional flags for the compiler.
CFLAGS =

# "-include" includes other makefiles.
#   * Each file listed in $(DEPS) is a makefile with one rule for making a .o
#     file with an empty recipe (no commands).
#   * Doing this adds those dependencies to that .o file, without interfering
#     with any recipes we give later to the same target.
#   * Recipes for this are given below with "$(OBJ_DIR)/%.o" as the target.
#   * Dependencies are generated at the same time the files are compiled. That
#     may seem wrong and weird, but it's fine. The dependencies can't change
#     without also changing the files that were previously dependencies.
-include $(OBJS:.o=.d)

# Here we build the .o (object) file and generate the .d (dependency) files.
#   * A lot of "automatic variables" and other symbols have been used. These can
#     be a little confusing to interpret. Here's what they do.
#       * '@' tells make not to echo the command to the command line. We do this
#         for basic maintenance operations to avoid spamming the command line.
#       * '$<' anywhere it appears in the recipe is substituted for the target.
#       * '$@' anywhere it appears in the recipe is substituted for the first
#         dependency. The first dependency wasn't strictly necessary because of
#         "-include $(DEPS)" above, but putting it in allows us to do this.
#       * '$*' anywhere it appears in the recipe will be substituted with
#         whatever was matched by '%' in the target. The '%' character cannot be
#         used in a recipe, but we are allowed to use this instead.
#   * We need to make any directories where our targets are located, otherwise
#     compilation will fail.
#       * mkdir makes a directory, and the -p option tells it to do two very
#         useful things.
#           * Recursively create any parent directories as needed. This allows
#             it to work easily to any depth.
#           * Don't complain if the directory already exists. This allows us to
#             not have to worry about when mkdir needs to be run and when it
#             doesn't.
#       * The shell command 'dirname' strips off the last part of a filepath.
#         Here we strip off the last part of whatever was matched by '%' in the
#         target.
#   * When generating object files, we want to also generate dependencies as a
#     side effect. The options starting with '-M' will do this with the gcc and
#     g++ compilers.
#       * '-MMD' does a lot. It combines '-M', '-MM', and '-MD'. It also
#         disables '-E'.
#           * Generate dependencies of the object as makefile rules.
#           * Exclude any dependencies on standard library headers.
#           * Output makefile results into a file. In particular, it takes the
#             file given to the '-o' flag, which in this case is '$@' (the
#             target), removes the '.o', and puts a '.d' in its place.
#           * Compilation continues as normal after generating the makefile
#             info. This is not the default behavior when using most '-M' flags.
#             This is intended to cut down on compilation time.
#       * '-MP' tells g++ to generated so called "phony" targets for header
#         files. This is a weird trick that will stop make from complaining if
#         the header no longer exists because it has been moved or deleted.
#         these targets will have no listed dependencies and no recipe.
#       * '-MT' takes a parameter to set the name of the target in the
#         generated rule. Without this, the compiler will strip off the
#         directories from the target name and use that. '-MT' overrides this
#         behavior. Here we give '$@' so that the generated target has the same
#         name as the target we are currently running.
$(OBJ_DIR)/%.o: ./%.cc
	@mkdir -p $(shell dirname $@)
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@ -MMD -MP -MT '$@'

# Generate binary.
all: $(DEV_OBJS)
	$(CC) $^ $(LDFLAGS) -o $(BIN_DIR)/classify

# Compile and run tests.
test: $(TEST_OBJS)
	$(CC) $^ $(LDFLAGS) -o $(BIN_DIR)/test
	$(BIN_DIR)/test

# Cleanup
clean:
	rm -f $(BIN_DIR)/*
	rm -rf $(BUILD_DIR)/*
