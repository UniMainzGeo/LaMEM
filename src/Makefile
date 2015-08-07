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
#     filename:   nlsolve.c
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

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

# Define PETSc-based C-compiler command
CCOMPILER = ${PCC} ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS}

#====================================================

# Define list of LaMEM library source files
CSRC = \
	LaMEMLib.c \
	Parsing.c \
	fdstag/fdstag.c \
	fdstag/lsolve.c \
	fdstag/nlsolve.c \
	fdstag/paraViewOutBin.c \
	fdstag/paraViewOutSurf.c \
	fdstag/paraViewOutMark.c \
	fdstag/advect.c \
	fdstag/marker.c \
	fdstag/cvi.c \
	fdstag/AVD.c \
	fdstag/constEq.c \
	fdstag/dfzero.c \
	fdstag/JacRes.c \
	fdstag/JacResTemp.c \
	fdstag/interpolate.c \
	fdstag/outFunct.c \
	fdstag/scaling.c \
	fdstag/matrix.c \
	fdstag/multigrid.c \
	fdstag/check_fdstag.c \
	fdstag/bc.c \
	fdstag/break.c \
	fdstag/surf.c \
	fdstag/input.c \
	fdstag/matProps.c \
	fdstag/tssolve.c \
	fdstag/tools.c \
	fdstag/sub_comm.c \
	fdstag/pc_semiredundant.c \
	fdstag/objFunct.c \
	fdstag/AVDView.c
	
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

# Create dependencies
OBJDIR = \
	../lib \
	../lib/${mode} \
	../bin \
	../bin/${mode} \
	../dep \
	../dep/${mode}

#====================================================

# Default target (build everything)

.DEFAULT_TARGET: all
.PHONY : all
all : lamem

#====================================================
# Build dependencies
.PHONY : depend
depend : ${OBJDIR}

${OBJDIR} :
	@mkdir -p ${OBJDIR}

#====================================================
# Build LaMEM library
.PHONY : lamemlib
lamemlib : ${LaMEM_LIB}

${LaMEM_LIB} : ${OBJDIR} ${COBJ}
	@echo "............................................."
	@echo ".......... Building LaMEM Library ..........."
	@echo "............................................."
	ar cr $@ ${COBJ}
	ranlib $@

#====================================================

# Link LaMEM executable
.PHONY : lamem
lamem : ${LaMEM}

${LaMEM} : ${LaMEM_LIB} ${LaMEM_OBJ}
	@echo "............................................."
	@echo "......... Linking LaMEM Executable .........."
	@echo "............................................."
	${CLINKER} ${LaMEM_OBJ} ${LaMEM_LIB} ${PETSC_LIB} -o $@

#====================================================

# Pattern rules for automatic generation of object & dependency files
# Insert full path to object files in dependency files with sed command
# NOTE: IBM XL compiler generates dependency as a by-product of compilation

../lib/${mode}/%.o : %.c Makefile Makefile.in
	${CCOMPILER} ${LAMEM_FLAGS} -c $< -o $@ 
ifeq ($(PLATFORM), ppc64)
	@mv -f ../lib/${mode}/$*.d ../dep/${mode}/$*.d.tmp
else
	@${CCOMPILER} ${DEPEN_FLAGS} ${LAMEM_FLAGS} $< -o ../dep/${mode}/$*.d.tmp
endif
	@sed '1s,^,../lib/${mode}/,' < ../dep/${mode}/$*.d.tmp > ../dep/${mode}/$*.d
	@rm -f ../dep/${mode}/$*.d.tmp

../lib/${mode}/%.o : fdstag/%.c Makefile Makefile.in
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

.PHONY : print
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

.PHONY : clean_all
clean_all :
	@echo "............................................."
	@echo ".......... Performing full clean ............"
	@echo "............................................."
	@rm -f ../lib/$(mode)/*
	@rm -rf ../bin/$(mode)/*
	@rm -f ../dep/$(mode)/*
	@rm -f *.o *.a LaMEM
	@rm -f fdstag/*.o
	@rm -f ../lib/*.a
	@rm -f ../dep/*.d

#====================================================

# Create automatic documentation
.PHONY : doc
doc:
	PDFLATEX=$(PDFLATEX) BIBTEX=$(BIBTEX) ../doc/Manual/./CreateDevelDoc.sh

#====================================================
