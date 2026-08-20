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
CXXFLAGS			= -pedantic # -Wextra
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
SRC_FILES = src/bcfio.cpp 
OBJ_FILES = $(subst src, $(BUILD_DIR), $(SRC_FILES:.cpp=.o))
APP_DEPS = $(OBJ_FILES:.o=.d)

TEST_DATA_DIR = $(PWD)/inst/exdata

TEST_DIR = tests/cpp
TEST_SRC = $(wildcard $(TEST_DIR)/test_*.cpp)
TEST_OBJ = $(subst $(TEST_DIR), $(BUILD_DIR), $(TEST_SRC:.cpp=.o))

######################################################################
# DO NOT EDIT BELOW
######################################################################

-include $(APP_DEPS)

.PHONY: hwas_cpp
hwas: build/bcfio.o 

build/bcfio.o: src/bcfio.cpp | $(BUILD_DIR)
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



# $(TEST_TARGET_PRG): $(TEST_DIR)/main.cpp $(TEST_OBJS) $(APP_OBJS) | $(TARGET)
# 	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(CXXLIBFLAGS) -o $@ $^ -lgtest -lhts
# 
# $(BUILD_DIR)/test_%.o: $(TEST_DIR)/test_%.cpp
# 	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(OBJ_OUTPUT_OPTIONS) $<
# 
# TEST_GRM_PRG = $(BUILD_DIR)/test_grm
# test_grm: $(TEST_GRM_PRG)
# 	./$(TEST_GRM_PRG)
# 
# $(TEST_GRM_PRG): $(BUILD_DIR)/test_grm.o $(BUILD_DIR)/grm.o | $(BUILD_DIR)
# 	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(CXXLIBFLAGS) -o $@ $^ -lgtest -lgtest_main
# 
# data: | $(TEST_DATA_DST)
# 
# $(BUILD_DIR)/geno_test_data%: $(TEST_DIR)/geno_test_data%
# 	rsync -avz $< $(BUILD_DIR)/
# 
# # tests: $(BUILD_DIR)/test_log #$(BUILD_DIR)/test_argparse
# # 
# # $(BUILD_DIR)/test_log: $(BUILD_DIR)/test_log.o $(BUILD_DIR)/logger.o ~/.local/lib/libgtest.a
# # 	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(CXXLIBFLAGS) -o $@ $^
# # 
# # $(BUILD_DIR)/test_argparse: $(BUILD_DIR)/test_argparse.o \
# # 	$(BUILD_DIR)/argparse.o \
# # 	~/.local/lib/libgtest.a
# # 	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) -I$(LOCAL_INCLUDE) -L$(LOCAL_LIB) -o $@ $^
# 
# # $(TEST_OBJS): $(TEST_SRC)
# #	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) -I$(LOCAL_INCLUDE) $(OBJ_OUTPUT_OPTIONS) $<
# #
# # $(BUILD_DIR)/test_log.o: $(TEST_DIR)/test_log.cpp
# # 	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(OBJ_OUTPUT_OPTIONS) $^
# 
# 
# # $(TEST_OBJS): $(TEST_SRC)
# #	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) -I$(LOCAL_INCLUDE) -L$(LOCAL_LIB) -o $@ $^
# 
# # $(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(CXXLIBFLAGS) $(OBJ_OUTPUT_OPTIONS) $^
# 
# 
# ######################################################################
# # 
# ######################################################################
# 
# check:
# 	./$(TEST_TARGET_PRG)
# 
# ######################################################################
# # 
# ######################################################################
# 
# 
# -include $(APP_DEPS)
# -include $(TEST_DEPS)
# 
# .PHONY: help
# help:
# 	-@echo "build grm"
# 	-@echo "2025 Palmer Lab"
# 	-@echo ""
# 	-@echo "make grm executable"
# 	-@echo "make libargparse"


######################################################################
# install
######################################################################

# install:
# 	dir_header=$${prefix%/}/include/stitchr; \
# 	if [ ! -d $${dir_header} ]; then \
# 		mkdir -p $${dir_header}; \
# 	fi; \
# 	for hfile in $$(ls $(HEADER_DIR)); do \
# 		cp $$hfile $${dir_header}/$${hfile}; \
# 	done; \
# 	 \
# 	dir_lib=$${prefix%/}/lib; \
# 	if [ ! -d $${dir_lib} ]; then \
# 		mkdir -p $${dir_lib}; \
# 	fi; \
# 	for libfile in $$(ls $(BUILD_DIR)/*.a); do
# 		cp $$libfile $${dir_lib}/$${libfile}; \
# 	done; \
# 	 \
# 	dir_bin = $${prefix%/}/bin; \
# 	if [ ! -d $${dir_bin} ]; then \
# 	   mkdir -p $${dir_bin}; \
# 	fi; \
# 	cp $(TARGET) $${dir_bin}/$(notdir $(TARGET))
# 


