#====================================================
#
# Makefile for reparitioning subroutine
#
#====================================================

# Include platform-specific constants

include ../../src/Makefile.in

#====================================================

# Include PETSc files

#====================================================

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


all:: repart

repart: DomainRepartitioning.o
	@echo "............................................."
	@echo "..... Linking Repartition Executable ........"
	@echo "............................................."
	${CLINKER} -o Repartition DomainRepartitioning.o ${PETSC_LIB}
	mv -f Repartition ../../bin

#====================================================

# Pattern rule for object file generation

%.o: %.cpp
	${PCXX} ${PCXX_FLAGS} ${CFLAGS} ${CCPPFLAGS} ${WRN_FLAGS_CXX} ${OPT_FLAGS} ${LAMEM_FLAGS} ${LAMEM_INCS} -c $< -o $@ 

#====================================================
clean_obj :
	@echo "............................................."
	@echo "........... Cleaning Object Files ..........."
	@echo "............................................."
	rm -f *.o
clean_all :
	@echo "............................................."
	@echo "........... Cleaning Object Files ..........."
	@echo "............................................."
	rm -f *.o
	@echo "............................................."
	@echo "........ Cleaning Executable ................"
	@echo "............................................."
	rm -f LaMEM ../bin/Repartition

#====================================================

