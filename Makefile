CXX = g++
CXXFLAGS = -O3
MPICXX= mpicxx
MPICXXFLAGS=-O3 -D__MPI

LIBDIR = /home/jzw2/mylib/OpenBLAS-0.3.23
LIBS = -lopenblas

SRCDIR = source
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = abasuc

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS)  $(LIBS) -o $@

mpi: $(OBJECTS)
	$(MPICXX) $(MPICXXFLAGS) $(OBJECTS) $(LIBS) -o $@
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(SRCDIR)/*.o $(EXECUTABLE) mpi

.PHONY: all clean
