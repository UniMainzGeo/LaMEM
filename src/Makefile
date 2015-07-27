#====================================================
#
# Makefile for LaMEM
# developed since 2011 at the Johannes Gutenberg University of Mainz
#
# by 
#	Anton Popov
#	Boris Kaus
#	Tobias Baumann
#	Adina Püsök
#	Arthur Bauville
#	Naiara Fernandez	
#
# Revision by Anton Popov
# Revision by Tobias Baumann
#
# NOTE: See Makefile.in for platform dependent definitions
#
#====================================================

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

# Define list of LAMEM library source files
CSRC = \
	LaMEMLib.c \
	LaMEMLib_FDSTAG.c \
	LaMEMLib_FDSTAG_private.c \
	LaMEM_Particles.c \
	LaMEM_Temperature.c \
	StokesOperators.c \
	Elements.c \
	Quadrature.c \
	Solvers.c \
	Breakpoint.c \
	Output.c \
	ParaViewOutput.c \
	NonDimensionalisation.c \
	Material.c \
	Mesh.c \
	Assembly.c \
	Assembly_FDSTAG.c \
	Utils.c \
	Parsing.c \
	Attributes.c \
	LaMEMVelPressureDA.c \
	LaMEMVelPressureDA_Q2PM1.c \
	LaMEMVelPressureDA_Q1P0.c \
	LaMEMVelPressureDA_Q1Q1.c \
	LaMEMVelPressureDA_FDSTAG.c \
	LaMEM_Initialize.c \
	LaMEM_AnalyticalSolutions.c \
	ApplyBoundaryConditions.c \
	AVDPhaseViewer.c \
	GetGravityField.c \
	GetSurfaceVelocity.c \
	LaMEM_FE_ErosionCode.c \
	Addition_FDSTAG.c \
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

#====================================================

# Default target (build everything)

.DEFAULT_TARGET: all
.PHONY : all
all : lamem


#====================================================


#====================================================

# Build LaMEM library
.PHONY : lamemlib
lamemlib : ${LaMEM_LIB}

${LaMEM_LIB} : ${COBJ}
	@mkdir -p ../bin
	@mkdir -p ../bin/opt
	@mkdir -p ../bin/deb
	@mkdir -p ../dep
	@mkdir -p ../dep/opt
	@mkdir -p ../dep/deb
	@mkdir -p ../lib
	@mkdir -p ../lib/opt
	@mkdir -p ../lib/deb
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
# Alternative generation of dependency files as a separate rule.
# Has a disadvantage that deeper modification is necessary to
# the files produced by compiler, i.e. dependency itself should
# be included in the targets. Also, the files are always built,
# even when they are not required.
#../dep/${mode}/%.d : %.c
#	@${CCOMPILER} ${DEPEN_FLAGS} ${LAMEM_FLAGS} $< -o $@.tmp
#	@sed -e '1s,^,../lib/${mode}/,' -e 's,:, $@ : ,' < $@.tmp > $@
#	@rm -f $@.tmp
#====================================================
# This test is usually unnecessary
# make fails at include stage if wrong mode is given 
#ifneq ($(mode), deb)
#ifneq ($(mode), opt)
#	@echo "............................................."
#	@echo "Error! mode="$(mode)" is invalid"
#	@echo "Should be either mode=deb or mode=opt"
#	@echo "............................................."
#	@exit 1
#endif
#endif
#====================================================
