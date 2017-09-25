CXX = g++
CXXFLAGS = -O2

ASTYLE_DIR = /opt/astyle/3.0.1

SOURCES = Array.cpp BoundaryCond.cpp Common.cpp heat_solver.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = heat-solver

build: realclean $(OBJECTS) $(EXECUTABLE)
	make -s clean

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXECUTABLE):
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

clean:
	rm -f $(OBJECTS)

realclean: clean
	rm -f $(EXECUTABLE)

format:
	$(ASTYLE_DIR)/bin/astyle --options=$(ASTYLE_DIR)/file/google.ini \
	                         --verbose --formatted *.cpp *.hpp
