CC = g++
FLAGS = -g -O2 -std=c++11 -Wall
LIBS = -lgsl -lstdc++ -lm
INCLUDE = -I./inc
OBJECTS = src/migdalMC.o src/Target.o

default: migdalMC

migdalMC: $(OBJECTS)
	$(CC) $(FLAGS) $^ -o $@ $(INCLUDE) $(LIBS)

src/%.o: src/%.cpp
	$(CC) $< $(FLAGS) $(INCLUDE) -c -o $@

clean:
	-rm src/*.o
	-rm -f ./migdalMC
