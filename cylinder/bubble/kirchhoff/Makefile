CC = g++
CXXFLAGS = -std=c++14 -g
INCLUDE = -I /home/magu/library/boost_1_77_0
TARGET = main

$(TARGET): $(TARGET).cpp
	$(CC) $(CXXFLAGS) $^ $(INCLUDE)

phony: clean

clean:
	@rm $(TARGET)
