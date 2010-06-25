all: kdtree.so kdtree3.so move_particles.so

BOOST = $(wildcard ~)/boost_1_38_0

move_particles.so: move_particles.cpp kdtree.so materials.cpp particles.cpp\
		   mesh.cpp kdtree.o kdtree3.o particles.hpp statistics.cpp statistics.hpp
	g++ -g -fPIC -O2 -c -Wall mesh.cpp
	g++ -g -fPIC -O2 -c -Wall materials.cpp
	g++ -g -fPIC -O2 -c -Wall particles.cpp
	g++ -g -fPIC -O2 -c -Wall move_particles.cpp
#	g++ -g -fPIC -O2 -c -Wall statistics.cpp -I$(BOOST) #killing for now
	ld -bundle -flat_namespace -undefined suppress \
			-o move_particles.so move_particles.o \
			 kdtree.o kdtree3.o materials.o

kdtree.o: kdtree.c  
	gcc -fPIC -O2 -c -Wall kdtree.c

kdtree3.o: kdtree3.c  
	gcc -fPIC -O2 -c -Wall kdtree3.c


kdtree.so: kdtree.o test_kdtree.c 
	gcc -Wall -o test_kdtree test_kdtree.c kdtree.o
	ld -bundle -flat_namespace -undefined suppress -o kdtree.so kdtree.o

kdtree3.so: kdtree3.o
	ld -bundle -flat_namespace -undefined suppress -o kdtree3.so kdtree3.o


driftscatter.so: driftscatter.c 
	gcc -fPIC -O2 -c -Wall driftscatter.c
	ld -bundle -flat_namespace -undefined suppress \
		-o driftscatter.so driftscatter.o
