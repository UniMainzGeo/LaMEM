#################################################################################
#
#       Petsc Makefile
#
#################################################################################

include ${PETSC_DIR}/conf/base

# compliler flags
CFLAGS=		-g
INCLUDE=	-I./include
PETSC_INCS=	${PETSC_INCLUDE}
MPI_INCS=	-I/${MPI_DIR}/include

PETSC_LIBS=	${PETSC_LIB}
MPI_LIBS=	-L/${MPI_DIR}/lib -lmpich -lpmpich

# Libraries for linking
LFLAGS= $(PETSC_LIBS) $(MPI_LIBS) -lm


OBJ	=	verify-da-decomp

all: ${OBJ}

${OBJ}: ${OBJ}.o
	-${CLINKER} -o ${OBJ} ${OBJ}.o ${PETSC_LIB}

clean_all:
	rm ${OBJ}
	rm ${OBJ}.o

