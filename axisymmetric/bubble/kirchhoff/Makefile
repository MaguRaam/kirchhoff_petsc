include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


TARGET = kirchhoff

all: $(TARGET)

$(TARGET) : $(TARGET).o
	-${CXXLINKER} -o $(TARGET) $(TARGET).o  ${PETSC_LIB}
	$(RM) $(TARGET).o

distclean :
	$(RM) $(TARGET) *.dat *.info *.png
