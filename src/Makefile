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
#
# LaMEM compilation command:
#  make mode=deb all (compile debug version of LaMEM and put in /bin/deb)
#  make mode=opt all (compile optimized version of LaMEM and put in /bin/opt)
#  make all          (compile optimized version of LaMEM and put in /bin/opt)

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

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Define PETSc-based C++ compiler command
CCOMPILER = ${CXX} ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

#====================================================

# Define list of LaMEM library source files
# List files to be excluded after filter-out
CSRC = $(filter-out LaMEM.cpp LaMEM_C.cpp, $(wildcard *.cpp))

# Add C wrapper source
#CSRC += $(C_WRAPPER_SOURCES)

# Objects generated from C_headers/*.cpp
#C_HEADERS_OBJ := $(addprefix ../lib/${mode}/, $(notdir $(wildcard C_headers/*.cpp:.cpp=.o)))
C_HEADERS_OBJ := $(patsubst C_headers/%.cpp, ../lib/${mode}/%.o, $(wildcard C_headers/*.cpp))

#====================================================

# Generate lists of library object files (exclude C_HEADERS_OBJ):
COBJ := $(addprefix ../lib/${mode}/, $(notdir $(filter-out C_headers/*.cpp, $(CSRC:.cpp=.o))))

# Generate lists of library dependency files:
CDEP := $(addprefix ../dep/${mode}/, $(notdir $(CSRC:.cpp=.o:.o=.d)))

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
dylib: ${COBJ} ${C_HEADERS_OBJ}
	@echo "............................................."
	@echo ".......... Building Dynamic Library ........"
	@echo "............................................."
ifeq ($(shell uname), Darwin)
	$(CCOMPILER) -dynamiclib -undefined dynamic_lookup -o ../lib/$(mode)/liblamem.dylib ${COBJ} ${C_HEADERS_OBJ} ${PETSC_LIB}
else
	$(CCOMPILER) -shared -fPIC -o ../lib/$(mode)/liblamem.so ${COBJ} ${C_HEADERS_OBJ} ${PETSC_LIB}
endif
#====================================================

# Link LaMEM executable

lamem : lamemlib linkexe

linkexe : ${LaMEM}

${LaMEM} : ${LaMEM_LIB} ${LaMEM_OBJ}
	@echo "............................................."
	@echo "......... Linking LaMEM Executable .........."
	@echo "............................................."
	${CXXLINKER} ${LaMEM_OBJ} ${LaMEM_LIB} ${PETSC_LIB} ${CLIB_FLAGS} -o $@

#====================================================

# Pattern rules for automatic generation of object & dependency files
# Insert full path to object files in dependency files with sed command
# NOTE: IBM XL compiler generates dependency as a by-product of compilation

../lib/${mode}/%.o : %.cpp
	${CCOMPILER} ${LAMEM_FLAGS} -c $< -o $@
ifeq ($(PLATFORM), ppc64)
	@mv -f ../lib/${mode}/$*.d ../dep/${mode}/$*.d.tmp
else
	@${CCOMPILER} ${DEPEN_FLAGS} ${LAMEM_FLAGS} $< -o ../dep/${mode}/$*.d.tmp
endif
	@sed '1s,^,../lib/${mode}/,' < ../dep/${mode}/$*.d.tmp > ../dep/${mode}/$*.d
	@rm -f ../dep/${mode}/$*.d.tmp

# Pattern rule for C_headers files
../lib/${mode}/%.o : C_headers/%.cpp
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
	@echo "LAMEM_LIB   : " ${LAMEM_LIB}
	@echo "............................................."
	@echo "CCOMPILER   : " ${CCOMPILER}
	@echo "............................................."
	@echo "CLINKER     : " ${CLINKER}
	@echo "............................................."
	@echo "PETSC_LIB   : " ${PETSC_LIB}
	@echo "............................................."
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
	@rm -f $(C_HEADERS_OBJ)

#====================================================

debug_sources:
	@echo "C_WRAPPER_SOURCES: $(C_WRAPPER_SOURCES)"
	@echo "CSRC: $(CSRC)"
	@echo "COBJ: $(COBJ)"

debug_headers:
	@echo "C_HEADERS_OBJ: $(C_HEADERS_OBJ)"

debug_wildcard:
	@echo "C_headers/*.cpp: $(wildcard C_headers/*.cpp)"

# Create automatic documentation

doc:
	PDFLATEX=$(PDFLATEX) BIBTEX=$(BIBTEX) ../doc/Manual/./CreateDevelDoc.sh

#====================================================