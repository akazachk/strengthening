### Makefile for strengthening project
# If there are errors, debug with
#   make --just-print
# or check defined variables with
#   make --print-data-base
# You can also use
#   make --warn-undefined-variables
### Shell type ###
# REMEMBER: use hard tabs only in a makefile
HOSTNAME := $(shell hostname)
UNAME := $(shell uname)
ARCH := $(shell uname -m)
ifeq ($(ARCH),x86_64)
	ARCH=x86-64
endif
ifeq ($(UNAME),Linux)
  CC     = g++
  SYSTEM = ${ARCH}_linux
endif
ifeq ($(UNAME),Darwin)
  CC     = clang++
  SYSTEM = ${ARCH}_osx
endif
RM = rm -f

### Build type ###
# Choose 'debug' or 'release'
# Can also be chosen through make "BUILD_CONFIG=XX" from command line 
# Or one can call make debug or make release directly
BUILD_CONFIG = release
BUILD_CONFIG = debug

### Variables user should set ###
ifeq ($(REPOS_DIR),)
	REPOS_DIR=${PWD}/..
endif
PROJ_DIR=${PWD}
COIN_VERSION = trunk
ifeq (${COIN_OR_HOME},)
	COIN_OR = $(PROJ_DIR)/lib/Cbc-$(COIN_VERSION)
else
	COIN_OR = ${COIN_OR_HOME}
endif
EIG_LIB = $(REPOS_DIR)/eigen
VPC_DIR = ${REPOS_DIR}/vpc

ifeq ($(USER),otherperson)
  #COIN_OR = enter/dir/here
  #EIG_LIB = enter/dir/here

  # Optional (for testing branch and bound or enabling certain functions):
  #GUROBI_DIR = enter/dir/here
  #GUROBI_LINK="gurobi80"
  #CPLEX_DIR = enter/dir/here
endif

# rossobianco
ifeq ($(USER),kazaalek)
  GUROBI_LINK = gurobi91
	GUROBI_DIR = ${GUROBI_HOME}
  #CPLEX_DIR = /home/ibm/cplex/20.1/cplex
	CPLEX_DIR = ${CPLEX_HOME}
	COIN_OR = /local_workspace/$(USER)/coin-or/Cbc-$(COIN_VERSION)
	EIG_LIB = ${HOME}/repos/eigen
endif

# rupert, HiPerGator
ifeq ($(USER),akazachkov)
	ifeq ($(UNAME),Linux)
	  not_hpg =
		ifeq ($(HOSTNAME),ISE-D41L3Q3)
      # w401
			not_hpg = true
		endif
		ifeq ($(HOSTNAME),rupert0)
			not_hpg = true
		endif
		ifeq ($(HOSTNAME),rupert1)
			not_hpg = true
		endif
		ifeq ($(HOSTNAME),rupert2)
			not_hpg = true
		endif
		ifeq ($(HOSTNAME),rupert3)
			not_hpg = true
		endif
		ifeq ($(HOSTNAME),rupert4)
			not_hpg = true
		endif
		ifeq ($(HOSTNAME),rupert5)
			not_hpg = true
		endif
		ifeq ($(HOSTNAME),rupert6)
			not_hpg = true
		endif

		ifdef not_hpg
			GUROBI_LINK = gurobi100
			GUROBI_DIR = ${GUROBI_HOME}
			CPLEX_DIR = ${CPLEX_HOME}
		else
			# HiPerGator
			COIN_OR = ${HOME}/repos/coin-or/Cbc-$(COIN_VERSION)
			GUROBI_LINK = gurobi95
			GUROBI_DIR = ${GUROBI_LOCAL}
			CPLEX_DIR = ${CPLEX_HOME}
			CONDA_LIB = ${HOME}/.conda/envs/vpc/lib
		endif
	endif
	# MacStudio
  ifeq ($(UNAME),Darwin)
    GUROBI_LINK = gurobi110
    GUROBI_DIR = ${GUROBI_HOME}
    CPLEX_DIR = ${CPLEX_HOME}
  endif
endif

ifeq ($(USER),akazachk)
	# ComputeCanada
  ifeq ($(UNAME),Linux)
	  COIN_OR = ${HOME}/projects/def-alodi/$(USER)/coin-or/Cbc-$(COIN_VERSION)
    GUROBI_LINK = gurobi91
    GUROBI_DIR = ${GUROBI_LOCAL}
    CPLEX_DIR = ${CPLEX_HOME}
  endif
	# Mac
  ifeq ($(UNAME),Darwin)
    GUROBI_LINK = gurobi10
    GUROBI_DIR = ${GUROBI_HOME}
    CPLEX_DIR = ${CPLEX_HOME}
		#COIN_OR = $(PROJ_DIR)/../coin-or/Cbc-$(COIN_VERSION)
		#COIN_OR = $(PROJ_DIR)/../vpc/lib/Cbc-$(COIN_VERSION)
  endif
endif

# Options for solvers
USE_COIN   = 1
USE_CLP    = 1
USE_EIGEN  = 1
USE_CBC    = 1
USE_GUROBI = 1
USE_CPLEX  = 0
USE_CLP_SOLVER = 1
USE_CPLEX_SOLVER = 0
CALC_COND_NUM = 0

# Concerning executable
EXECUTABLE_STUB = main
SRC_DIR = src
SOURCES = main.cpp
DIR_LIST = $(SRC_DIR) $(SRC_DIR)/cut $(SRC_DIR)/utility

# Code version
CODE_VERSION    = $(shell git log -1 --pretty=format:"%H")
VPC_VERSION     = $(shell git -C ${VPC_DIR} log -1 --pretty=format:"%H")
VPC_CBC_VERSION = $(shell git -C ${COIN_OR}/Cbc log -1 --pretty=format:"%H")
VPC_CLP_VERSION = $(shell git -C ${COIN_OR}/Clp log -1 --pretty=format:"%H")

SOURCES += \
		cut/CglAdvCut.cpp \
    cut/disjcuts.cpp \
    cut/gmic.cpp \
    cut/regularity.cpp \
    cut/strengthen.cpp \
		utility/eigen.cpp \
		utility/verify.cpp

# VPC directories
VPC_SRC_DIR = ${VPC_DIR}/src
VPC_DIR_LIST = branch cut disjunction utility

VPC_SOURCES += \
		branch/CbcBranchStrongDecision.cpp \
		branch/CbcCompareBFS.cpp \
		branch/OsiChooseStrongCustom.cpp \
    utility/nbspace.cpp \
    utility/OsiProblemData.cpp \
		utility/SolverHelper.cpp \
		utility/utility.cpp \
		utility/VPCSolverInterface.cpp \
		branch/VPCEventHandler.cpp \
		cut/CglVPC.cpp \
		cut/CutHelper.cpp \
    cut/PRLP.cpp \
    disjunction/Disjunction.cpp \
    disjunction/PartialBBDisjunction.cpp \
    disjunction/SplitDisjunction.cpp \
    disjunction/VPCDisjunction.cpp

# For running tests (need not include these or main if releasing code to others)
DIR_LIST += $(SRC_DIR)/test
SOURCES += test/analysis.cpp test/BBHelper.cpp
ifeq ($(CALC_COND_NUM), 1)
	VPC_SOURCES += test/condition_number.cpp
endif

### Set build values based on user variables ###
ifeq ($(BUILD_CONFIG),debug)
  # "Debug" build - no optimization, include debugging symbols, and keep inline functions
	SOURCES += utility/strengthening_debug.cpp
  VPC_SOURCES += utility/vpc_debug.cpp utility/debug.cpp
  OUT_DIR = Debug
  DEBUG_FLAG = -g3
  OPT_FLAG = -O0
  DEFS = -DTRACE -DPRINT_LP_WITH_CUTS -DVPC_DEBUG -DCODE_VERSION="\#${CODE_VERSION}" -DVPC_VERSION="\#${VPC_VERSION}"
  # message-length sets line wrapping for error messages; 0 = no line wrapping
	# PIC stands for "position-independent code" to generate machine code that
	# can be loaded at different memory addresses, such as by using relative rather than absolute jumps
	# which is needed for shared libraries
  EXTRA_FLAGS = -fmessage-length=0 -fPIC
  ifeq ($(CC),g++)
    ifneq ($(USE_CPLEX),1)
      EXTRA_FLAGS += -fkeep-inline-functions 
    endif
  endif
endif
ifeq ($(BUILD_CONFIG),release)
  # "Release" build - maximum optimization, no debug symbols
  OUT_DIR = Release
  DEBUG_FLAG = 
  OPT_FLAG = -O3
  DEFS = -DCODE_VERSION="\#${CODE_VERSION}" -DVPC_VERSION="\#${VPC_VERSION}"
  EXTRA_FLAGS = -fmessage-length=0 -ffast-math -fPIC
endif
ifeq ($(USE_COIN),1)
  DEFS += -DUSE_COIN
endif
ifeq ($(USE_CLP),1)
  DEFS += -DUSE_CLP
  DEFS += -DVPC_CLP_VERSION="\#${VPC_CLP_VERSION}"
endif
ifeq ($(USE_CLP_SOLVER),1)
  DEFS += -DUSE_CLP_SOLVER
endif
ifeq ($(USE_CBC),1)
  DEFS += -DUSE_CBC
  DEFS += -DVPC_CBC_VERSION="\#${VPC_CBC_VERSION}"
endif
ifeq ($(USE_EIGEN),1)
  DEFS += -DUSE_EIGEN
endif
ifeq ($(USE_GUROBI),1)
  DEFS += -DUSE_GUROBI
  SOURCES += test/GurobiHelper.cpp
  GUROBI_INC="${GUROBI_DIR}/include"
  GUROBI_LIB="${GUROBI_DIR}/lib"
endif
ifeq ($(USE_CPLEX),1)
  DEFS += -DIL_STD -DUSE_CPLEX
  SOURCES += test/CplexHelper.cpp
  CPLEX_INC = $(CPLEX_DIR)/include
  CPLEX_LIB = $(CPLEX_DIR)/lib/$(SYSTEM)/static_pic
endif
ifeq ($(USE_CPLEX_SOLVER),1)
  DEFS += -DUSE_CPLEX_SOLVER
endif
ifeq ($(COIN_VERSION),2.10)
  DEFS += -DCBC_VERSION_210PLUS -DSAVE_NODE_INFO
endif
ifeq ($(COIN_VERSION),trunk)
  DEFS += -DCBC_VERSION_210PLUS -DCBC_TRUNK -DSAVE_NODE_INFO
endif
ifeq ($(CALC_COND_NUM),1)
	DEFS += -DCALC_COND_NUM
endif

EXECUTABLE = $(OUT_DIR)/$(EXECUTABLE_STUB)

# It is important that the only thing that changes about these directories
# is OUT_DIR depending on whether it is the debug or release build
# This is because later (building dependencies, archive file, cleaning, etc.)
# depends on this fact (when doing *_debug targets)
OBJ_DIR = $(OUT_DIR)/$(SRC_DIR)
OUT_DIR_LIST = $(addprefix $(OUT_DIR)/,$(DIR_LIST))
OBJECTS = $(SOURCES:.cpp=.o)
OUT_OBJECTS = $(addprefix $(OBJ_DIR)/,$(OBJECTS))

VPC_OBJ_DIR = $(OUT_DIR)/vpc
VPC_OUT_DIR_LIST = $(addprefix $(OUT_DIR)/vpc/,$(VPC_DIR_LIST))
VPC_OBJECTS = $(VPC_SOURCES:.cpp=.o)
VPC_OUT_OBJECTS = $(addprefix $(VPC_OBJ_DIR)/,$(VPC_OBJECTS))

# Set includes
APPLINCLS = -Iinclude -Iinclude/test
## TEMPORARY CHANGE TO PARENT VPC (NO WIFI)
APPLINCLS += -I${VPC_DIR}/include -I${VPC_DIR}/include/common

APPLLIB = -lm -lz -lbz2 -lreadline
ifeq ($(UNAME),Linux)
	ifeq ($(USER),akazachkov)
		ifneq (${CONDA_LIB},)
			APPLLIB += -L${CONDA_LIB}
	  endif
	endif
endif

# Linker
CFLAGS = -Wall -MMD -MP
CFLAGS += -m64 $(DEBUG_FLAG) $(OPT_FLAG) $(EXTRA_FLAGS)
CXXFLAGS = $(CFLAGS) -std=c++17
#CXXFLAGS = $(CFLAGS) -std=c++17 -Wextra -Wpedantic
CXXLINKFLAGS += -std=c++17
ifeq ($(CC),clang++)
  CXXFLAGS += -Wno-gnu-zero-variadic-macro-arguments
  #CXXFLAGS += -stdlib=libc++ 
  #CXXLINKFLAGS += -stdlib=libc++ 
  APPLLIB += -framework Accelerate
endif
ifeq ($(CC),g++)
  ifneq (${ENV_BLAS_LIB},)
    APPLLIB += -L${ENV_BLAS_LIB} -lblas
  endif
  ifneq (${ENV_LAPACK_LIB},)
    APPLLIB += -L${ENV_LAPACK_LIB} -llapack
  endif
endif

# Set up COIN-OR stuff
ifeq ($(USE_COIN),1)
	# If not defined for the environment, define CBC / BCP here
	ifeq ($(BUILD_CONFIG),debug)
		CBC = $(COIN_OR)/buildg
	endif
	ifeq ($(BUILD_CONFIG),release)
		CBC = $(COIN_OR)/build
	endif
	CBClib = $(CBC)/lib
	# When switching from svn to coinbrew, the new include directory is coin-or not coin
	ifeq ($(COIN_VERSION),trunk)
		CBCinc = $(CBC)/include/coin-or
  else
		CBCinc = $(CBC)/include/coin
	endif
	APPLINCLS += -isystem $(CBCinc)
	APPLLIB += -L$(CBClib)
  CXXLINKFLAGS += -Wl,-rpath $(CBClib)
	ifeq ($(USE_CBC),1)
		#APPLLIB += -lCbcSolver
		APPLLIB += -lCbc
	endif
  ifeq ($(USE_CLP),1)
    APPLLIB += -lOsiClp
  endif
  APPLLIB += -lCgl
  APPLLIB += -lOsi
  ifeq ($(USE_CLP),1)
    APPLLIB += -lClp
  endif
  APPLLIB += -lCoinUtils
endif
ifeq ($(USE_GUROBI),1)
  APPLINCLS += -isystem ${GUROBI_INC}
  APPLLIB   += -L${GUROBI_LIB} -lgurobi_c++ -l${GUROBI_LINK} -lm
endif
ifeq ($(USE_CPLEX),1)
  APPLINCLS      += -isystem "$(CPLEX_INC)"
  APPLLIB        += -L${CPLEX_LIB} -lilocplex -lcplex -lm -lpthread -ldl
  CXXLINKFLAGS   += -Wl,-rpath $(CPLEX_LIB)
  #APPLLIB       += -lOsiCpx
endif
ifeq ($(USE_EIGEN),1)
  APPLINCLS += -I$(EIG_LIB) -ISpectra
endif

### Targets ###
all: | directories $(EXECUTABLE)
debug: FORCE
	@$(MAKE) "BUILD_CONFIG=debug"
release: FORCE
	@$(MAKE) "BUILD_CONFIG=release"

$(EXECUTABLE): $(OUT_OBJECTS) $(VPC_OUT_OBJECTS)
		@echo ' '
		@echo 'Building target: $@'
		@echo 'Invoking' $(CC) 'linker'
		$(CC) $(DEFS) $(CXXLINKFLAGS) $(APPLINCLS) -o $@ $^ $(APPLLIB)

# Uses pattern matching to determine where source is located and where to put object
# By this rule, if OBJ_DIR=Debug/src, SRC_DIR=src, and
# the target is Debug/src/main.o, then source is automatically defined as Debug/src/main.cpp
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
		@echo ' '
		@echo 'Building target: $@'
		@echo 'Invoking' $(CC) 'compiler'
		$(CC) $(CXXFLAGS) $(DEFS) $(APPLINCLS) -c $< -o $@ 
		@echo 'Finished building: $@'
$(VPC_OBJ_DIR)/%.o: $(VPC_SRC_DIR)/%.cpp
		@echo ' '
		@echo 'Building target: $@'
		@echo 'Invoking' $(CC) 'compiler'
		$(CC) $(CXXFLAGS) $(DEFS) $(APPLINCLS) -c $< -o $@ 
		@echo 'Finished building: $@'

### Archive file ###
AR = ar
AR_FLAGS = -rcs
LIB_DIR=lib
archive_%:
	@$(MAKE) archive "BUILD_CONFIG=$*"

archive: | lib_directory $(OUT_DIR)/lib$(EXECUTABLE_STUB).a
	  @cp ${OUT_DIR}/lib${EXECUTABLE_STUB}.a ${LIB_DIR}/lib${EXECUTABLE_STUB}.a

$(OUT_DIR)/lib$(EXECUTABLE_STUB).a: $(OUT_OBJECTS)
		@echo ' '
		@echo 'Making archive file: $@'
		@echo 'Invoking archiver'
		$(AR) $(AR_FLAGS) $@ $^
		@echo 'Finished making archive'

### Shared object library file ###
SO_TARGET_STUB = $(LIB_DIR)/lib$(EXECUTABLE_STUB)
ifeq ($(UNAME),Linux)
  SO_TARGET = $(SO_TARGET_STUB).so
  SO_FLAGS = -shared
endif
ifeq ($(UNAME),Darwin)
  SO_TARGET = $(SO_TARGET_STUB).dylib
  #SO_FLAGS += -undefined dynamic_lookup
  SO_FLAGS += -dynamiclib
endif
shared_lib_%:
	@$(MAKE) shared_lib "BUILD_CONFIG=$*"

shared_lib: $(SO_TARGET)

$(SO_TARGET): $(OUT_OBJECTS)
		@echo ' '
		@echo 'Making shared object library file: $@'
		$(CC) $(DEFS) $(CXXLINKFLAGS) $(APPLINCLS) $(SO_FLAGS) -o $@ $^ $(APPLLIB)
		@echo 'Finished making shared object library file'

### Dependencies ###
# Dependencies (the -include says to ignore errors)
DEPENDENCIES = $(OUT_OBJECTS:.o=.d)
-include $(DEPENDENCIES)

### Phony ###
#.PHONY = \
#		all \
#		clean clean_debug distclean_debug clean_release distclean_release \
#		directories dir_debug dir_release dir_lib_debug dir_lib_release \
#		print print_dep \
#		archive_debug archive_release archive

.PHONY = all clean directories distclean doxygen print test

### Docs ###
doxygen: FORCE
	@doxygen

### Cleaning ###
clean_%: FORCE
	@$(MAKE) clean "BUILD_CONFIG=$*"
distclean_%: FORCE
	@$(MAKE) distclean "BUILD_CONFIG=$*"

clean: FORCE
	@$(RM) $(OUT_OBJECTS) $(EXECUTABLE) $(OUT_DIR)/lib$(EXECUTABLE_STUB).a

distclean: FORCE
	@$(RM) $(OUT_OBJECTS) $(EXECUTABLE) $(DEPENDENCIES) $(OUT_DIR)/lib$(EXECUTABLE_STUB).a $(LIB_DIR)/lib$(EXECUTABLE_STUB).a

### Making directories that you need ###
MKDIR_P = mkdir -p

dir_%: FORCE
	@$(MAKE) directories "BUILD_CONFIG=$*"
directories: $(OUT_DIR_LIST) $(VPC_OUT_DIR_LIST)
$(OUT_DIR_LIST):
	$(MKDIR_P) $(OUT_DIR_LIST)
$(VPC_OUT_DIR_LIST):
	$(MKDIR_P) $(VPC_OUT_DIR_LIST)

dir_lib_%: FORCE
	@$(MAKE) dir_lib "BUILD_CONFIG=$*"
dir_lib: FORCE
	$(MKDIR_P) $(LIB_DIR)
lib_directory: $(LIB_DIR)
$(LIB_DIR):
	$(MKDIR_P) $(LIB_DIR)

test_%: FORCE
	@$(MAKE) test "BUILD_CONFIG=$*"
test: FORCE
	@$(EXECUTABLE) -f test/bm23.mps -d 2 --strengthen=1 --gomory=-1

print: FORCE
	$(info UNAME: ${UNAME})
	$(info CC: ${CC})
	$(info COIN_OR: ${COIN_OR})
	$(info CPLEX_DIR: ${CPLEX_DIR})
	$(info GUROBI_DIR: ${GUROBI_DIR})
	$(info GUROBI_LINK: ${GUROBI_LINK})
	$(info USE_COIN: ${USE_COIN})
	$(info USE_CLP: ${USE_CLP})
	$(info USE_CBC: ${USE_CBC})
	$(info USE_GUROBI: ${USE_GUROBI})
	$(info USE_CPLEX: ${USE_CPLEX})
	$(info OUT_DIR: ${OUT_DIR})
	$(info DEBUG_FLAG: ${DEBUG_FLAG})
	$(info OPT_FLAG: ${OPT_FLAG})
	$(info DEFS: ${DEFS})
	$(info EXTRA_FLAGS: ${EXTRA_FLAGS})
	$(info LIB_DIR: ${LIB_DIR})
	$(info SOURCES: ${SOURCES})
	$(info OBJ_DIR: ${OBJ_DIR})
	$(info OUT_DIR_LIST: ${OUT_DIR_LIST})
	$(info OUT_OBJECTS: ${OUT_OBJECTS})
	$(info VPC_SOURCES: ${VPC_SOURCES})
	$(info VPC_OBJ_DIR: ${VPC_OBJ_DIR})
	$(info VPC_OUT_DIR_LIST: ${VPC_OUT_DIR_LIST})
	$(info VPC_OBJECTS: ${VPC_OBJECTS})
	$(info VPC_OUT_OBJECTS: ${VPC_OUT_OBJECTS})

print_dep: FORCE
	$(info DEPENDENCIES: ${DEPENDENCIES})

FORCE: 
