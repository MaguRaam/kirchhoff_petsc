CFLAGS     = -Wall
CPPFLAGS   =
LIBFILES   = 
TARGET     = hype
SRC        = $(wildcard *.c)
OBJ        = $(SRC:.c=.o)
CLEANFILES = $(TARGET)
LOCDIR     = $(CURDIR)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

all: $(TARGET)

$(TARGET): $(OBJ)  chkopts
	-${CLINKER} -o $(TARGET) $(OBJ) ${PETSC_LIB}
	${RM} $(OBJ)
	
distclean :
	$(RM) $(TARGET) *.vts *.png *.dat *.dat.info *.bin *.bin*