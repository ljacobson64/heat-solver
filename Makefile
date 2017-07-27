CXX = g++
CXXFLAGS = -O2

ASTYLE_DIR = /opt/astyle/3.0.1

all: clean build run plot

build: clean
	$(CXX) $(CXXFLAGS) -c Array.cpp
	$(CXX) $(CXXFLAGS) heat.cpp Array.o -o heat

run:
	./heat

plot:
	python plot.py

format:
	$(ASTYLE_DIR)/bin/astyle --options=$(ASTYLE_DIR)/file/google.ini \
	                         --verbose --formatted *.cpp *.hpp

clean:
	rm -f *.o *.txt heat
