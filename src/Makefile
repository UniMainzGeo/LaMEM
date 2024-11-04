#==============================================================================
#
#   Project      : LaMEM
#   License      : MIT, see LICENSE file for details
#   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
#
#==============================================================================

# Include platform-specific constants

include Makefile.in

#====================================================

# Include PETSc library
#
# Environmental variables for PETSc installation directories:
# PETSC_DEB = /directory/where/petsc/debug/is/installed
# PETSC_OPT = /directory/where/petsc/optimized/is/installed
# FASTSCAPELIB = /directory/where/FastScapeLib/is/located
#
# LaMEM compilation command:
#  make mode=deb all (compile debug version of LaMEM and put in /bin/deb)
#  make mode=opt all (compile optimized version of LaMEM and put in /bin/opt)
#  make all          (compile optimized version of LaMEM and put in /bin/opt)
#  make mode=optFS all surface=SURFACE=2 (compile optimized version of LaMEM with fastscape and put in /bin/optFS)
#  surf_mode = 1; don't use fastscape; surf_mode = 2; use fastscape

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
else ifeq ($(mode), optFS)
# ifneq (, $(filter,${FASTSCAPE_LIB} ${PETSC_OPT}))
ifeq (${FASTSCAPE_LIB},)
${error Environmental variable PETSC_OPT and FASTSCAPE_LIB must be set to installation directory, respectively}
else ifeq (${PETSC_OPT},)
${error Environmental variable PETSC_OPT and FASTSCAPE_LIB must be set to installation directory, respectively}
endif
PETSC_DIR = ${PETSC_OPT}
FASTSCAPE_DIR = ${FASTSCAPE_LIB}
else
$(error Unknown compilation mode specified)
endif

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Define PETSc-based C++ compiler command
CCOMPILER = ${CXX} ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

# Define FastScape-based fortran compiler command
GF = gfortran
# Define FastScapeFortranLib 
FS_LIB = -L${FASTSCAPE_DIR} -lfastscapelib_fortran

#====================================================

# Define list of LaMEM library source files
# List files to be excluded after filter-out
CSRC = $(filter-out LaMEM.cpp, $(wildcard *.cpp))
FS = $(wildcard *.f90)

#====================================================
ifeq ($(mode), optFS)
# Generate lists of library object files:
COBJ := $(addprefix ../lib/${mode}/, $(notdir $(CSRC:.cpp=.o) $(FS:.f90=.o)))

# Generate lists of library dependency files:
CDEP := $(addprefix ../dep/${mode}/, $(notdir $(CSRC:.cpp=.d) $(FS:.f90=.d)))

else
COBJ := $(addprefix ../lib/${mode}/, $(notdir $(CSRC:.cpp=.o)))
CDEP := $(addprefix ../dep/${mode}/, $(notdir $(CSRC:.cpp=.d)))
endif

# Get library and executable objects:
LaMEM = ../bin/${mode}/LaMEM
LaMEM_OBJ = ../lib/${mode}/LaMEM.o
LaMEM_DEP = ../dep/${mode}/LaMEM.d
LaMEM_LIB = ../lib/${mode}/liblamem.a

# If we build LaMEM for BinaryBuilder we link to a PETSc build that has a different name
ifeq ($(LAMEM_BINARYBUILDER), true)
	PETSC_LIB := $(filter-out -lpetsc, $(PETSC_LIB))
	PETSC_LIB += -lpetsc_double_real_Int32
endif
#====================================================

# Target list

.PHONY: builddir buildlib lamemlib linkexe lamem print clean_all doc 

# Default target (build executable)

.DEFAULT_TARGET: all

all : lamem

#====================================================

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

#====================================================

# force full rebuild when Makefiles are modified

../dep/$(mode)/makefile.stat : Makefile Makefile.in
	@rm -rf ../lib/$(mode)/*
	@rm -rf ../bin/$(mode)/*
	@rm -rf ../dep/$(mode)/*
	@echo 'makefile timestamp status file' > ../dep/$(mode)/makefile.stat

#====================================================

# Build LaMEM library

lamemlib : builddir buildlib

buildlib : ${LaMEM_LIB}

${LaMEM_LIB} : ../dep/$(mode)/makefile.stat ${COBJ}
	@echo "............................................."
	@echo ".......... Building LaMEM Library ..........."
	@echo "............................................."
	ar cr $@ ${COBJ}
	ranlib $@

# Create a dynamic library
dylib: ${COBJ}
ifeq ($(mode), optFS)
	$(CCOMPILER) -shared -fPIC -o ../lib/$(mode)/LaMEMLib.dylib ${COBJ}  ${PETSC_LIB} ${FS_LIB}
else
	$(CCOMPILER) -shared -fPIC -o ../lib/$(mode)/LaMEMLib.dylib ${COBJ}  ${PETSC_LIB} 
endif
#====================================================

# Link LaMEM executable

lamem : lamemlib linkexe

linkexe : ${LaMEM}

${LaMEM} : ${LaMEM_LIB} ${LaMEM_OBJ}
	@echo "............................................."
	@echo "......... Linking LaMEM Executable .........."
	@echo "............................................."
ifeq ($(mode), optFS)
	${CXXLINKER} ${LaMEM_OBJ} ${LaMEM_LIB} ${PETSC_LIB} ${CLIB_FLAGS} ${FS_LIB} -o  $@
else
	${CXXLINKER} ${LaMEM_OBJ} ${LaMEM_LIB} ${PETSC_LIB} ${CLIB_FLAGS} -o  $@
endif

#====================================================

# Pattern rules for automatic generation of object & dependency files
# Insert full path to object files in dependency files with sed command
# NOTE: IBM XL compiler generates dependency as a by-product of compilation
ifeq ($(mode), optFS)
../lib/${mode}/%.o : %.f90
	${GF} ${FS_LIB}  -c $< -o $@ 
endif

# define SURFACE 1 in Cpp script when there isn't a value setting in bash
surface=SURFACE=1

../lib/${mode}/%.o : %.cpp 
	${CCOMPILER} ${LAMEM_FLAGS} -D$(surface) -c $< -o $@ 

ifeq ($(PLATFORM), ppc64)
	@mv -f ../lib/${mode}/$*.d ../dep/${mode}/$*.d.tmp
else
	@${CCOMPILER} ${DEPEN_FLAGS} ${LAMEM_FLAGS} $< -o ../dep/${mode}/$*.d.tmp
endif
	@sed '1s,^,../lib/${mode}/,' < ../dep/${mode}/$*.d.tmp > ../dep/${mode}/$*.d
	@rm -f ../dep/${mode}/$*.d.tmp

#====================================================

# Include available dependency files
-include ${CDEP} ${LaMEM_DEP}

#====================================================

print :
	@echo "............................................."
	@echo "........ Environmental variables ............"
	@echo "............................................."
ifeq ($(PLATFORM), x86_64)
ifneq ($(shell mpicc --version 2>&1 | grep -c pgcc), 0)
	@echo "COMPILER    :  PORTLAND GROUP"
else ifneq ($(shell mpicc --version 2>&1 | grep -c gcc), 0)
	@echo "COMPILER    :  GNU "
else ifneq ($(shell mpicc --version 2>&1 | grep -c icc), 0)
	@echo "COMPILER    :  INTEL "
else
	@echo "COMPILER    :  UNKNOWN "
endif
endif
ifeq ($(PLATFORM), ppc64)
	@echo "COMPILER    :  IBM "
endif
	@echo "............................................."
	@echo "mode        : " ${mode}
	@echo "............................................."
	@echo "PETSC_DIR   : " ${PETSC_DIR}
	@echo "............................................."
	@echo "DEPEN_FLAGS : " ${DEPEN_FLAGS}
	@echo "............................................."
	@echo "LAMEM_FLAGS : " ${LAMEM_FLAGS}
	@echo "............................................."
	@echo "LAMEM_LIB   : " ${LAMEM_LIB}
	@echo "............................................."
	@echo "CCOMPILER   : " ${CCOMPILER}
	@echo "............................................."
	@echo "CLINKER     : " ${CLINKER}
	@echo "............................................."
	@echo "PETSC_LIB   : " ${PETSC_LIB}
	@echo "............................................."
ifeq ($(mode), optFS)
	@echo "............................................."
	@echo "FASTSCAPE_LIB   : " ${FS_LIB}
	@echo "............................................."
endif
	@echo "CLIB_FLAGS  : " ${CLIB_FLAGS}
	@echo "............................................."

#====================================================

clean_all :
	@echo "............................................."
	@echo ".......... Performing full clean ............"
	@echo "............................................."
	@rm -rf ../lib/$(mode)/*
	@rm -rf ../bin/$(mode)/*
	@rm -rf ../dep/$(mode)/*
#====================================================

# Create automatic documentation

doc:
	PDFLATEX=$(PDFLATEX) BIBTEX=$(BIBTEX) ../doc/Manual/./CreateDevelDoc.sh

#====================================================
