all: kdtree.so move_particles.so

move_particles.so: move_particles.cpp kdtree.so materials.cpp particles.cpp\
		   mesh.cpp kdtree.o particles.hpp
	g++ -fPIC -O2 -c -Wall mesh.cpp
	g++ -fPIC -O2 -c -Wall move_particles.cpp
	g++ -fPIC -O2 -c -Wall materials.cpp
	g++ -fPIC -O2 -c -Wall particles.cpp
	ld -bundle -flat_namespace -undefined suppress \
			-o move_particles.so move_particles.o \
						kdtree.o materials.o mesh.o\
						particles.o

kdtree.o: kdtree.c  
	gcc -fPIC -O2 -c -Wall kdtree.c

kdtree.so: kdtree.o test_kdtree.c 
	gcc -Wall -o test_kdtree test_kdtree.c kdtree.o
	ld -bundle -flat_namespace -undefined suppress -o kdtree.so kdtree.o

driftscatter.so: driftscatter.c 
	gcc -fPIC -O2 -c -Wall driftscatter.c
	ld -bundle -flat_namespace -undefined suppress \
		-o driftscatter.so driftscatter.o
