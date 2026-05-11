#==============================================================================
#
#   Project      : LaMEM
#   License      : MIT, see LICENSE file for details
#   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
#
#==============================================================================
# LaMEM compilation command:
#  make mode=deb all (compile debug version of LaMEM and put in /bin/deb)
#  make mode=opt all (compile optimized version of LaMEM and put in /bin/opt)
#  make all          (compile optimized version of LaMEM and put in /bin/opt)
#==============================================================================

# define compilation mode
# deb - debug 
# opt - optimized (default)

mode = opt

# define PETSc installation directories:
# PETSC_DEB = /directory/where/petsc/debug/is/installed
# PETSC_OPT = /directory/where/petsc/optimized/is/installed

ifeq ($(mode), deb)
   ifeq (${PETSC_DEB},)
      $(error Environmental variable PETSC_DEB must be set to PETSc debug installation directory)
   endif
   PETSC_DIR = ${PETSC_DEB}
else ifeq ($(mode), opt)
   ifeq (${PETSC_OPT},)
      $(error Environmental variable PETSC_OPT must be set to PETSc optimized installation directory)
   endif
   PETSC_DIR = ${PETSC_OPT}
else
   $(error Unknown compilation mode specified)
endif

# include PETSc environmental variables:
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# define platform
PLATFORM = $(shell uname -m)

# check if platform is supported
ifneq ($(findstring $(PLATFORM), x86_64 arm64 aarch64 armv6l armv7l i686 powerpc64le amd64), $(PLATFORM))
   $(error Unsupported platform)
endif

# define MPI compiler and linker
MPI_COMPILER = ${CXX} ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}
LINKER       = ${CXXLINKER}

# extract actual compiler
MPI_COMPILER_SHOW := $(shell \
	$(MPI_COMPILER) --showme:command 2>/dev/null || \
	$(MPI_COMPILER) --show 2>/dev/null || \
	$(MPI_COMPILER) -show 2>/dev/null )

CXX_COMPILER := $(word 1, $(filter-out ccache,$(MPI_COMPILER_SHOW)))

# set dependency flag (should be supported by all compilers)
DEPEN_FLAGS = -MM

# set common compiler flags (should be supported by all compilers)
LAMEM_FLAGS = -std=c++17 -Wall -Wextra -Wconversion \
              -Wpointer-arith -Wcast-align -Wwrite-strings -Wformat=2 \
              -Wundef -Wnon-virtual-dtor -Wno-unused-result -Wunused-but-set-variable

# tune compiler flags
ifeq ($(CXX_COMPILER), icpx)
   LAMEM_FLAGS += -Wimplicit-fallthrough -Wno-sign-conversion \
   -Wno-float-conversion -Wno-deprecated-declarations -Wno-vla-cxx-extension
else ifeq ($(CXX_COMPILER), g++)
   LAMEM_FLAGS += -Wduplicated-cond -Wduplicated-branches\
                   -Wlogical-op -Wnull-dereference -Wuseless-cast 
else ifeq ($(CXX_COMPILER), clang++)
   LAMEM_FLAGS += -Wimplicit-fallthrough \
                  -Wshorten-64-to-32 -Wnewline-eof
else
    $(error Unsupported compiler)
endif

# add debuggiing flag if in debug mode
ifeq ($(mode), deb)
   LAMEM_FLAGS += -g
endif

# Windows requires this linker flag
UNAME := $(shell uname)
ifneq ($(filter MSYS%,$(UNAME)),)
   CLIB_FLAGS = -lssp
endif

#==============================================================================

# define list of LaMEM library source files
# list files to be excluded after filter-out
CSRC = $(filter-out LaMEM.cpp, $(wildcard *.cpp))

#==============================================================================

# generate lists of library object files:
COBJ := $(addprefix ../lib/${mode}/, $(notdir $(CSRC:.cpp=.o)))

# generate lists of library dependency files:
CDEP := $(addprefix ../dep/${mode}/, $(notdir $(CSRC:.cpp=.d)))

# get library and executable objects:
LAMEM = ../bin/${mode}/LaMEM
LAMEM_OBJ = ../lib/${mode}/LaMEM.o
LAMEM_DEP = ../dep/${mode}/LaMEM.d
LAMEM_LIB = ../lib/${mode}/liblamem.a

# if we build LaMEM for BinaryBuilder we link to a PETSc build that has a different name
ifeq ($(LAMEM_BINARYBUILDER), true)
   PETSC_LIB := $(filter-out -lpetsc, $(PETSC_LIB))
   PETSC_LIB += -lpetsc_double_real_Int32
endif
#==============================================================================

# target list

.PHONY: builddir buildlib lamemlib linkexe lamem print clean_all doc

# default target (build executable)

.DEFAULT_TARGET: all

all : lamem

#==============================================================================

# directory make rules

../bin :
	@if [ ! -d $@ ]; then mkdir $@; fi

../lib :
	@if [ ! -d $@ ]; then mkdir $@; fi

../dep :
	@if [ ! -d $@ ]; then mkdir $@; fi

../bin/${mode} :
	@if [ ! -d $@ ]; then mkdir $@; fi

../lib/${mode} :
	@if [ ! -d $@ ]; then mkdir $@; fi

../dep/${mode} :
	@if [ ! -d $@ ]; then mkdir $@; fi


# build directories

builddir : ../bin ../lib ../dep ../bin/${mode} ../lib/${mode} ../dep/${mode}

#==============================================================================

# force full rebuild when Makefiles are modified

../dep/$(mode)/makefile.stat : Makefile
	@rm -rf ../lib/$(mode)/*
	@rm -rf ../bin/$(mode)/*
	@rm -rf ../dep/$(mode)/*
	@echo 'makefile timestamp status file' > ../dep/$(mode)/makefile.stat

#==============================================================================

# build LaMEM library

lamemlib : builddir buildlib

buildlib : ${LAMEM_LIB}

${LAMEM_LIB} : ../dep/$(mode)/makefile.stat ${COBJ}
	@echo "............................................."
	@echo ".......... Building LaMEM Library ..........."
	@echo "............................................."
	ar cr $@ ${COBJ}
	ranlib $@

# create a dynamic library
dylib: ${COBJ}
	$(MPI_COMPILER) -shared -fPIC -o ../lib/$(mode)/LaMEMLib.dylib ${COBJ}  ${PETSC_LIB} 

#==============================================================================

# link LaMEM executable

lamem : lamemlib linkexe

linkexe : ${LAMEM}

${LAMEM} : ${LAMEM_LIB} ${LAMEM_OBJ}
	@echo "............................................."
	@echo "......... Linking LaMEM Executable .........."
	@echo "............................................."
	${LINKER} ${LAMEM_OBJ} ${LAMEM_LIB} ${PETSC_LIB} ${CLIB_FLAGS} -o $@

#==============================================================================

# pattern rules for automatic generation of object & dependency files
# insert full path to object files in dependency files with sed command

../lib/${mode}/%.o : %.cpp
	${MPI_COMPILER} ${LAMEM_FLAGS} -c $< -o $@
	@${MPI_COMPILER} ${DEPEN_FLAGS} ${LAMEM_FLAGS} $< -o ../dep/${mode}/$*.d.tmp
	@sed '1s,^,../lib/${mode}/,' < ../dep/${mode}/$*.d.tmp > ../dep/${mode}/$*.d
	@rm -f ../dep/${mode}/$*.d.tmp

#==============================================================================

# include available dependency files
-include ${CDEP} ${LAMEM_DEP}

#==============================================================================

print :
	@echo "............................................."
	@echo "........ Environmental variables ............"
	@echo "............................................."
	@echo "mode         : " ${mode}
	@echo "............................................."
	@echo "PLATFORM     : " ${PLATFORM}
	@echo "............................................."
	@echo "PETSC_DIR    : " ${PETSC_DIR}
	@echo "............................................."
	@echo "MPI_COMPILER : " ${MPI_COMPILER}
	@echo "............................................."
	@echo "CXX_COMPILER : " ${CXX_COMPILER}
	@echo "............................................."
	@echo "LINKER       : " ${LINKER}	
	@echo "............................................."
	@echo "DEPEN_FLAGS  : " ${DEPEN_FLAGS}
	@echo "............................................."
	@echo "LAMEM_FLAGS  : " ${LAMEM_FLAGS}
	@echo "............................................."
	@echo "LAMEM_LIB    : " ${LAMEM_LIB}
	@echo "............................................."
	@echo "PETSC_LIB    : " ${PETSC_LIB}
	@echo "............................................."
	@echo "CLIB_FLAGS   : " ${CLIB_FLAGS}
	@echo "............................................."

#==============================================================================

clean_all :
	@echo "............................................."
	@echo ".......... Performing full clean ............"
	@echo "............................................."
	@rm -rf ../lib/$(mode)/*
	@rm -rf ../bin/$(mode)/*
	@rm -rf ../dep/$(mode)/*

#==============================================================================
