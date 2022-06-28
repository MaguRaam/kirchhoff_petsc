CC = g++
CXXFLAGS = -std=c++14
INCLUDE = -I /home/magu/library/boost_1_77_0
DEPS =
OBJ = area.o

%.o: %.cpp $(DEPS)
	$(CC) $(CXXFLAGS) -c $< $(INCLUDE)

a.out: $(OBJ)
	$(CC) $(CXXFLAGS) $^

phony: clean distclean

clean:
	@rm a.out $(OBJ)

distclean:
	@rm *.png *.dat
