COMPILER_OPTIONS := -g -Wall
CXXFLAGS := -std=c++11 $(COMPILER_OPTIONS)
objects = $(patsubst %.cpp, %, $(wildcard *.cpp))

all: $(objects)

clean: 
	rm  $(objects)


