CATCH = /usr/local/include/catch/single_include
EIGEN = /usr/local/include/eigen

CXX      = g++
CC       = $(CXX)
CPPFLAGS = -g -O2
CFLAGS   = -g -O2
CXXFLAGS = -c -std=c++11 
DFLAGS   = -D UNIT_TEST

EXEC         = dream
SOURCES      = $(wildcard *.cpp)
OBJECTS      = $(SOURCES:.cpp=.o)

$(EXEC): $(OBJECTS) $(OBJECTS_LIB)
	$(CXX) $(OBJECTS) $(CFLAGS) -o $(EXEC)

%.o: %.cpp
	$(CXX) $(DFLAGS) $(CPPFLAGS) $(CXXFLAGS) -I $(EIGEN) -I $(CATCH) $< -o $@

clean:
	rm -f $(EXEC) $(OBJECTS)
