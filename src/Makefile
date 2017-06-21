	#==============================================================================
#
#     Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
#     All rights reserved.
#
#     This software was developed at:
#
#          Institute of Geosciences
#          Johannes-Gutenberg University, Mainz
#          Johann-Joachim-Becherweg 21
#          55128 Mainz, Germany
#
#     project:    LaMEM
#     filename:   Makefile
#
#     LaMEM is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published
#     by the Free Software Foundation, version 3 of the License.
#
#     LaMEM is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#     See the GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
#
#
#     Contact:
#         Boris Kaus       [kaus@uni-mainz.de]
#         Anton Popov      [popov@uni-mainz.de]
#
#
#     Main development team:
#          Anton Popov      [popov@uni-mainz.de]
#          Boris Kaus       [kaus@uni-mainz.de]
#          Tobias Baumann
#          Adina Pusok
#          Arthur Bauville
#
#==============================================================================

# Include platform-specific constants

include Makefile.in

#====================================================

# Include PETSc library
ifeq ($(mode), deb)
	PETSC_DIR = ${PETSC_DEB}
endif

ifeq ($(mode), opt)
	PETSC_DIR = ${PETSC_OPT}
endif

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Define PETSc-based C++ compiler command
CCOMPILER = ${CXX} ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

#====================================================

# Define list of LaMEM library source files
# List files to be excluded after filter-out
CSRC = $(filter-out LaMEM.c, $(wildcard *.c))

#====================================================

# Generate lists of library object files:
COBJ := $(addprefix ../lib/${mode}/, $(notdir $(CSRC:.c=.o)))

# Generate lists of library dependency files:
CDEP := $(addprefix ../dep/${mode}/, $(notdir $(CSRC:.c=.d)))

# Get library and executable objects:
LaMEM = ../bin/${mode}/LaMEM
LaMEM_OBJ = ../lib/${mode}/LaMEM.o
LaMEM_DEP = ../dep/${mode}/LaMEM.d
LaMEM_LIB = ../lib/${mode}/liblamem.a

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

#====================================================

# Link LaMEM executable

lamem : lamemlib linkexe

linkexe : ${LaMEM}

${LaMEM} : ${LaMEM_LIB} ${LaMEM_OBJ}
	@echo "............................................."
	@echo "......... Linking LaMEM Executable .........."
	@echo "............................................."
	${CLINKER} ${LaMEM_OBJ} ${LaMEM_LIB} ${PETSC_LIB} -o $@

#====================================================

# Pattern rules for automatic generation of object & dependency files
# Insert full path to object files in dependency files with sed command
# NOTE: IBM XL compiler generates dependency as a by-product of compilation

../lib/${mode}/%.o : %.c
	${CCOMPILER} ${LAMEM_FLAGS} -c $< -o $@
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
	@echo "CCOMPILER   : " ${CCOMPILER}
	@echo "............................................."
	@echo "CLINKER     : " ${CLINKER}
	@echo "............................................."
	@echo "PETSC_LIB   : " ${PETSC_LIB}
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

