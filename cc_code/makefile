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

# Set the source directory and find all source files.
#   * The shell function runs the command following it in the shell.
#   * The find command will recursively search all subdirectories of $(SRC_DIR)
#     for files named "*.cc" and return the results.
#   * We do similar work to set $(TEST_SRCS).
SRC_DIR = ./src
TEST_DIR = ./test
MAIN = $(SRC_DIR)/main/main.cc
TEST_MAIN = $(TEST_DIR)/test-module.cc
ALL_SRCS = $(filter-out $(MAIN) $(TEST_MAIN), $(shell find $(SRC_DIR)/ -name '*.cc'))
SRCS = $(filter-out %.test.cc, $(ALL_SRCS))
TEST_SRCS = $(filter %.test.cc, $(ALL_SRCS))

BUILD_DIR = ./build

# Set the object file directory and define all names of .o (object) files.
#   * We use the $(SRCS) variable to define these names.
#   * '%' is a wildcard that can be used to capture the part of the expression
#     we want to keep.
#   * We replace $(SRC_DIR) with $(OBJ_DIR) and .cc with .o.
#   * We do similarly to set $(TEST_OBJS).
OBJ_DIR = $(BUILD_DIR)/obj/dev
OBJS = $(SRCS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
MAIN_OBJ = $(BUILD_DIR)/mains/obj/main.o
TEST_OBJ_DIR = $(BUILD_DIR)/obj/test
TEST_OBJS = $(TEST_SRCS:$(SRC_DIR)/%.test.cc=$(TEST_OBJ_DIR)/%.o)
TEST_MAIN_OBJ = $(BUILD_DIR)/mains/obj/test.o

# Set the make dependency file directory and define all names of .d (dependency)
# files.
#   * The formula for $(DEPS) is very similar to the one for $(OBJS), see above.
#   * These dependency files are included directly by the makefile to determine
#     which dependencies are needed for each object file.
#   * This is related to "-include $(DEPS)" below.
#   * Here we go ahead and put both tests and non-tests into the same variable,
#     because there is no reason we need them separate.
DEP_DIR = $(BUILD_DIR)/autodep/dev
TEST_DEP_DIR = $(BUILD_DIR)/autodep/test
MAIN_DEP = $(BUILD_DIR)/mains/autodep/main.d
TEST_MAIN_DEP = $(BUILD_DIR)/mains/autodep/test.d
DEPS = $(SRCS:$(SRC_DIR)/%.cc=$(DEP_DIR)/%.d)
DEPS += $(TEST_SRCS:$(SRC_DIR)/%.test.cc=$(TEST_DEP_DIR)/%.d)
DEPS += $(MAIN_DEP)
DEPS += $(TEST_MAIN_DEP)

# Set the output directory.
BIN_DIR = ./bin

# Set the include paths to search for headers.
#   * -iquote is for headers included with quotes
#   * -I is for headers included with angle brackets.
INCLUDE = -iquote$(SRC_DIR)
INCLUDE += -I./extern/mrmpi-7Apr14/src/
INCLUDE += -I./extern/mrmpi-7Apr14/mpistubs/
TEST_INCLUDE = $(INCLUDE)
TEST_INCLUDE += -I./extern/boost_1_73_0_subset/
TEST_INCLUDE += -I./extern/FakeIt/single_header/boost/

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
-include $(DEPS)

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
#   * The mv command puts the .d (dependency) file in its proper place. When
#     building the .o (objects) and .d (dependencies) at the same time, gcc and
#     g++ have no option to put them in different directories, so we have to
#     manually move the file if we want it somewhere different.
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(OBJ_DIR)/$(shell dirname $*)
	@mkdir -p $(DEP_DIR)/$(shell dirname $*)
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@ -MMD -MP -MT '$@'
	@mv $(OBJ_DIR)/$*.d $(DEP_DIR)/$*.d

$(MAIN_OBJ): $(MAIN)
	@mkdir -p $(shell dirname $(MAIN_OBJ))
	@mkdir -p $(shell dirname $(MAIN_DEP))
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@ -MMD -MP -MT '$@'
	@mv $(MAIN_OBJ:.o=.d) $(MAIN_DEP)

$(TEST_OBJ_DIR)/%.o: $(SRC_DIR)/%.test.cc
	@mkdir -p $(TEST_OBJ_DIR)/$(shell dirname $*)
	@mkdir -p $(TEST_DEP_DIR)/$(shell dirname $*)
	$(CC) -c $(CFLAGS) $(TEST_INCLUDE) $< -o $@ -MMD -MP -MT '$@'
	@mv $(TEST_OBJ_DIR)/$*.d $(TEST_DEP_DIR)/$*.d

$(BUILD_DIR)/mains/obj/test.o: $(TEST_MAIN)
	@mkdir -p $(shell dirname $(TEST_MAIN_OBJ))
	@mkdir -p $(shell dirname $(TEST_MAIN_DEP))
	$(CC) -c $(CFLAGS) $(TEST_INCLUDE) $< -o $@ -MMD -MP -MT '$@'
	@mv $(TEST_MAIN_OBJ:.o=.d) $(TEST_MAIN_DEP)

# Generate binary.
all: $(OBJS) $(MAIN_OBJ)
	$(CC) $^ $(LDFLAGS) -o $(BIN_DIR)/classify

# Compile and run tests.
test: $(OBJS) $(TEST_OBJS) $(TEST_MAIN_OBJ)
	$(CC) $^ $(LDFLAGS) -o $(BIN_DIR)/test
	$(BIN_DIR)/test

# Cleanup
clean:
	rm -f $(BIN_DIR)/*
	rm -rf $(BUILD_DIR)/*
