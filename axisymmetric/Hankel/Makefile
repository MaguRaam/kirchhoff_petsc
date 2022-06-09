CC = g++
CXXFLAGS = -std=c++14
INCLUDE = -I /home/magu/library/boost_1_77_0


Hankel: Hankel.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) $^ -o $@

Kirchhoff: Kirchhoff.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) $^ -o $@
	
test_Kirchhoff: test_Kirchhoff.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) $^ -o $@



phony: clean distclean

clean:
	@rm Hankel Kirchhoff test_Kirchhoff

distclean:
	@rm *.png *.dat
