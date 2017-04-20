CXX = g++
CXXFLAGS = -O2

heat: clean
	$(CXX) $(CXXFLAGS) -c Array.cpp
	$(CXX) $(CXXFLAGS) heat.cpp Array.o -o heat

clean:
	rm -f *.o heat
