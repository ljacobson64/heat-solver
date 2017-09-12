CXX = g++
CXXFLAGS = -O2

ASTYLE_DIR = /opt/astyle/3.0.1

all: clean build

build: clean
	$(CXX) $(CXXFLAGS) heat_solver.cpp Array.cpp -o heat_solver

run:
	./heat_solver

plot:
	python plot.py

format:
	$(ASTYLE_DIR)/bin/astyle --options=$(ASTYLE_DIR)/file/google.ini \
	                         --verbose --formatted *.cpp *.hpp

clean:
	rm -f *.o *.txt heat_solver
