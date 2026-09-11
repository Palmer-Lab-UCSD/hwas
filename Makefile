#
# 2025 Palmer Lab
#
######################################################################
# machine dependent options
######################################################################
# CXXFLAGS note: Remember that -g flag is for generating source-level 
# debug info.

# OBJ_OUTPUT_OPTIONS: compiler options.  Clang, and I presume 
# also gcc, support the creation of dependency files (-MMD) and (-MP) phony
# targets required for constructing an object file.  The (-o) and $@
# are the standard output file designation submitted to the compiler,
# and $@ is an automatic variable storing the rules target.
# library archive program
#

ifneq ($(shell which clang++),)
CXX					= clang++
# CXXFLAGS			= -pedantic -fsanitize=address
CXXFLAGS			= -pedantic
else ifneq ($(shell which g++),)
CXX					= g++
CXXFLAGS			= -Wpedantic -Wextra 
else
$(error "Couldn't establish either clang or gcc compiler availability")
endif


CXXFLAGS			+= -g -std=c++17 -Wall -Werror

ifndef VIM
CXXFLAGS += -fdiagnostics-color=always
endif

# Recall that -c flag prevents the compiler linking object files
OBJ_OUTPUT_OPTIONS 	= -c -MMD -MP -o $@
AR 					= ar
AR_FLAGS 			= crs

LOCAL_LIB			= $(HOME)/local/lib
LOCAL_LD			= $(HOME)/local/include

HEADER_DIR 			= $(PWD)/inst/include

CXXLDFLAGS = $(addprefix -I, $(HEADER_DIR) $(LOCAL_LD) $(CXXLD))
CXXLIBFLAGS = $(addprefix -L, $(LOCAL_LIB) $(CXXLIB))

######################################################################
# DIRECTORIES
######################################################################

BUILD_DIR = build
SRC_DIR = src
SRC_FILES = $(SRC_DIR)/grm.cpp $(SRC_DIR)/bcfio.cpp
OBJ_FILES = $(subst $(SRC_DIR), $(BUILD_DIR), $(SRC_FILES:.cpp=.o))
APP_DEPS = $(OBJ_FILES:.o=.d)

TEST_DATA_DIR = $(PWD)/inst/exdata

TEST_DIR = tests/cpp
TEST_SRC = $(wildcard $(TEST_DIR)/test_*.cpp)
TEST_OBJ = $(subst $(TEST_DIR), $(BUILD_DIR), $(TEST_SRC:.cpp=.o))

######################################################################
# DO NOT EDIT BELOW
######################################################################


.PHONY: hwas
hwas: $(OBJ_FILES)

-include $(APP_DEPS)

$(OBJ_FILES): $(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(OBJ_OUTPUT_OPTIONS) $<

$(BUILD_DIR):
	mkdir $@

######################################################################
# Test Build Rules
######################################################################

.PHONY: tests
tests: $(BUILD_DIR)/runtests 
	./$(BUILD_DIR)/runtests

$(BUILD_DIR)/runtests: $(TEST_DIR)/main.cpp $(TEST_OBJ) $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(CXXLIBFLAGS) -o $@ $^ -lgtest -lhts

$(BUILD_DIR)/test_%.o: $(TEST_DIR)/test_%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(OBJ_OUTPUT_OPTIONS) $<


######################################################################
# static analysis
######################################################################

.PHONY: static
static: $(SRC_FILES)
	clang --analyze $(CXXFLAGS) $(CXXLDFLAGS) $^

RCPP_HEADER = $(shell R --no-echo -e 'cat(system.file("include", package="Rcpp"))')
RCPP_EIG_HEADER = $(shell R --no-echo -e 'cat(system.file("include", package="RcppEigen"))')

RSYS_LDFLAGS = $(shell R CMD config --cppflags)

.PHONY: static_rbcfio
static_rbcfio: src/hwas_bcfio.cpp
	clang --analyze -std=c++17 -I$(RCPP_HEADER) $(RSYS_LDFLAGS) -I$(LOCAL_LD) -I$(HEADER_DIR) $<


.PHONY: static_rgrm
static_rgrm: src/hwas_grm.cpp
	clang --analyze -std=c++17 \
		-I$(RCPP_HEADER) \
		$(RSYS_LDFLAGS) \
		-I$(LOCAL_LD) \
		-I$(HEADER_DIR) \
		$<

.PHONY: static_rpgsim
static_rpgsim: src/hwas_rpgsim.cpp
	clang --analyze -std=c++17 \
		-I$(RCPP_HEADER) \
		-I$(RCPP_EIG_HEADER) \
		$(RSYS_LDFLAGS) \
		-I$(LOCAL_LD) \
		-I$(HEADER_DIR) \
		$^


