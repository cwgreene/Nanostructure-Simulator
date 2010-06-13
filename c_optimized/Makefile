all: kdtree.so move_particles.so

BOOST = $(wildcard ~)/boost_1_38_0

move_particles.so: move_particles.cpp kdtree.so materials.cpp particles.cpp\
		   mesh.cpp kdtree.o particles.hpp statistics.cpp statistics.hpp
	g++ -g -fPIC -O2 -c -Wall mesh.cpp
	g++ -g -fPIC -O2 -c -Wall move_particles.cpp
	g++ -g -fPIC -O2 -c -Wall materials.cpp
	g++ -g -fPIC -O2 -c -Wall particles.cpp
	g++ -g -fPIC -O2 -c -Wall statistics.cpp -I$(BOOST)
	ld -bundle -flat_namespace -undefined suppress \
			-o move_particles.so move_particles.o \
						kdtree.o materials.o mesh.o\
						particles.o statistics.o

kdtree.o: kdtree.c  
	gcc -fPIC -O2 -c -Wall kdtree.c

kdtree.so: kdtree.o test_kdtree.c 
	gcc -Wall -o test_kdtree test_kdtree.c kdtree.o
	ld -bundle -flat_namespace -undefined suppress -o kdtree.so kdtree.o

driftscatter.so: driftscatter.c 
	gcc -fPIC -O2 -c -Wall driftscatter.c
	ld -bundle -flat_namespace -undefined suppress \
		-o driftscatter.so driftscatter.o
