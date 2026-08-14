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
#  make mode=deb surf=scape (compile debug version of LaMEM with fastscape)
#  surf_mode = 1 in .dat file; don't use fastscape; surf_mode = 2 in .dat file; use fastscape 
#==============================================================================

# define compilation mode
# deb - debug 
# opt - optimized (default)

mode = opt
surf = none

# define PETSc installation directories:
# PETSC_DEB = /directory/where/petsc/debug/is/installed
# PETSC_OPT = /directory/where/petsc/optimized/is/installed
# FASTSCAPELIB = /directory/where/FastScapeLib/is/located

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

ifeq ($(surf), scape)
   ifeq (${FASTSCAPE_LIB},)
      ${error Environmental variable FASTSCAPE_LIB must be set to installation directory}
   endif
   FASTSCAPE_DIR = ${FASTSCAPE_LIB}
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

CXX_COMMAND  := $(word 1, $(filter-out ccache,$(MPI_COMPILER_SHOW)))
CXX_COMPILER := $(shell readlink -f $$(which $(CXX_COMMAND)))

# set dependency flag (should be supported by all compilers)
DEPEN_FLAGS = -MM

# set common compiler flags (should be supported by all compilers)
LAMEM_FLAGS = -std=c++17 -Wall -Wextra -Wconversion \
              -Wpointer-arith -Wcast-align -Wwrite-strings -Wformat=2 \
              -Wundef -Wnon-virtual-dtor -Wno-unused-result -Wunused-but-set-variable

# tune compiler flags
ifeq ($(findstring clang++, $(CXX_COMPILER)), clang++)
   LAMEM_FLAGS += -Wimplicit-fallthrough \
                  -Wshorten-64-to-32 -Wnewline-eof
else ifeq ($(findstring icpx, $(CXX_COMPILER)), icpx)
   LAMEM_FLAGS += -Wimplicit-fallthrough -Wno-sign-conversion \
                  -Wno-float-conversion -Wno-deprecated-declarations -Wno-vla-cxx-extension
else ifeq ($(findstring g++, $(CXX_COMPILER)), g++)
   LAMEM_FLAGS += -Wduplicated-cond -Wlogical-op -Wnull-dereference
else
    $(error Unsupported compiler: $(CXX_COMPILER))
endif

# fastscape part
ifeq ($(surf), scape)
LAMEM_FLAGS += -DWITH_FASTSCAPE
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

# Define FastScapeFortranLib
# FS_LIB = -L${FASTSCAPE_DIR} -lfastscapelib_fortran -lgfortran
# (embed an rpath so the resulting binary finds libfastscapelib_fortran at
# runtime without requiring LD_LIBRARY_PATH/DYLD_FALLBACK_LIBRARY_PATH to be set)
FS_LIB = ${CC_LINKER_SLFLAG}${FASTSCAPE_DIR} -L${FASTSCAPE_DIR} -lfastscapelib_fortran
#====================================================

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

.PHONY: builddir buildlib lamemlib linkexe lamem print clean_all doc format check checkformat checkgetrestore

# default target (build executable)

.DEFAULT_TARGET: all

all : lamem

#==============================================================================

# format source files
#
# NOTE: different astyle versions format the code differently, so the version
# is pinned here and verified before formatting/checking. Use the same one as
# the CI (the version shipped with Ubuntu 22.04/24.04). Homebrew currently
# installs a much newer astyle, which reformats files that this version
# considers correct (and vice versa).

ASTYLE_VERSION = 3.1

checkastyle:
	@if ! command -v astyle > /dev/null 2>&1; then \
		echo "ERROR: astyle not found; LaMEM formatting requires astyle ${ASTYLE_VERSION}."; \
		exit 1; \
	fi
	@if ! astyle --version 2>&1 | grep -q "Version ${ASTYLE_VERSION}"; then \
		echo "............................................."; \
		echo "ERROR: wrong astyle version."; \
		echo "  required : ${ASTYLE_VERSION}"; \
		echo "  found    : `astyle --version 2>&1`"; \
		echo "Different astyle versions format the code differently, so using"; \
		echo "another version would reformat unrelated files and break the CI."; \
		echo "............................................."; \
		exit 1; \
	fi

format: checkastyle
	@echo "Formatting source files ..."
	@astyle --options=.astylerc *.cpp *.h
	@echo "Done."

#==============================================================================

# check that source files are formatted according to .astylerc
# (fails if `make format` would modify any file; leaves the tree untouched)

checkformat:
	@echo "Checking source formatting ..."
	@if astyle --dry-run --options=.astylerc *.cpp *.h | grep -q '^Formatted'; then \
		echo "............................................."; \
		echo "ERROR: source files are not properly formatted."; \
		astyle --dry-run --options=.astylerc *.cpp *.h | grep '^Formatted'; \
		echo "Run 'make format' locally and commit the changes."; \
		echo "............................................."; \
		exit 1; \
	else \
		echo "All source files are properly formatted."; \
	fi

#==============================================================================

# check that every PETSc VecGetArray*/DMDAVecGetArray* call (and DOF/Write
# variants) is matched by the corresponding Restore call within its function
# (leaves the tree untouched; see scripts/petsc_getrestore_check.jl)

GETRESTORE_CHECK = ../scripts/petsc_getrestore_check.jl

checkgetrestore:
	@echo "Checking PETSc Get/Restore Array pairing ..."
	@julia ${GETRESTORE_CHECK} .

#==============================================================================

# run all static checks

check: checkformat checkgetrestore

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
ifeq ($(surf), scape)
	$(MPI_COMPILER) -shared -fPIC -o ../lib/$(mode)/LaMEMLib.dylib ${COBJ}  ${PETSC_LIB} ${FS_LIB}
else
	$(MPI_COMPILER) -shared -fPIC -o ../lib/$(mode)/LaMEMLib.dylib ${COBJ}  ${PETSC_LIB} 
endif
#==============================================================================

# link LaMEM executable

lamem : lamemlib linkexe

linkexe : ${LAMEM}

${LAMEM} : ${LAMEM_LIB} ${LAMEM_OBJ}
	@echo "............................................."
	@echo "......... Linking LaMEM Executable .........."
	@echo "............................................."
ifeq ($(surf), scape)
	${LINKER} ${LAMEM_OBJ} ${LAMEM_LIB} ${PETSC_LIB} ${CLIB_FLAGS} ${FS_LIB} -o  $@
else
	${LINKER} ${LAMEM_OBJ} ${LAMEM_LIB} ${PETSC_LIB} ${CLIB_FLAGS} -o  $@
endif

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
ifeq ($(surf), scape)
	@echo "............................................."
	@echo "FASTSCAPE_LIB   : " ${FS_LIB}
	@echo "............................................."
endif

#==============================================================================

clean_all :
	@echo "............................................."
	@echo ".......... Performing full clean ............"
	@echo "............................................."
	@rm -rf ../lib/$(mode)/*
	@rm -rf ../bin/$(mode)/*
	@rm -rf ../dep/$(mode)/*

#==============================================================================
