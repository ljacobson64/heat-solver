CXX = g++
CXXFLAGS = -O2

ASTYLE_DIR = $$HOME/astyle

all: clean build run

build: clean
	$(CXX) $(CXXFLAGS) -c Array.cpp
	$(CXX) $(CXXFLAGS) heat.cpp Array.o -o heat

run:
	./heat

format:
	$(ASTYLE_DIR)/astyle --options=$(ASTYLE_DIR)/google.ini \
	                     --verbose \
	                     --formatted \
	                     *.cpp *.hpp

clean:
	rm -f *.o heat
