all: kdtree.so kdtree3.so move_particles.so

BOOST = `pkg-config --cflags boost`
EIGEN = `pkg-config --cflags eigen2`
EIGEN3 = /Users/christophergreene/src/eigen
CPPOPTIONS =-g -fPIC -O2 -c -Wall -x c++ 

Face: Face.hpp
	 g++ $(CPPOPTIONS) -I $(EIGEN) Face.hpp

tests: Boundary.hpp Boundary_Test.cpp
	g++ $(CPPOPTIOnS) $(EIGEN) Boundary_Test.cpp -o tests/boundarytest

testGeometry: 
	 g++ -g -O2 -Wall -x c++ $(EIGEN) testGeometry.cpp -o testGeometry

testPolytope: testPolytope.cpp Polytope.hpp Face.hpp
	 g++ -g -O2 -Wall -x c++ $(EIGEN) testPolytope.cpp -o testPolytope

move_particles.so: move_particles.cpp kdtree.o materials.cpp particles.cpp\
		   mesh.cpp kdtree.o particles.hpp statistics.cpp\
		   statistics.hpp mesh.hpp Face.hpp Polytope.hpp particles.hpp Boundary.hpp
#	g++ -g -fPIC -O2 -c -Wall -I $(EIGEN) mesh.cpp
	g++ -mstackrealign -m32 -fPIC -02 -c -Wall $(EIGEN) materials.cpp
#	g++ -g -fPIC -O2 -c -Wall -I $(EIGEN) particles.cpp
	g++ -mstackrealign -m32 -fPIC -02 -c -Wall  $(EIGEN) move_particles.cpp
#	g++ -g -fPIC -O2 -c -Wall statistics.cpp -I$(BOOST) #killing for now
	ld -bundle -flat_namespace -undefined suppress \
			-o move_particles.so move_particles.o \
			 kdtree.o kdtree3.o materials.o

kdtree.o: kdtree.c  
	gcc -m32 -fPIC -O2 -c -Wall kdtree.c

kdtree3.o: kdtree3.c  
	gcc -m32 -fPIC -O2 -c -Wall kdtree3.c


kdtree.so: kdtree.o test_kdtree.c 
	gcc -m32 -Wall -o test_kdtree test_kdtree.c kdtree.o
	ld -bundle -flat_namespace -undefined suppress -o kdtree.so kdtree.o

kdtree3.so: kdtree3.o
	ld -bundle -flat_namespace -undefined suppress -o kdtree3.so kdtree3.o


driftscatter.so: driftscatter.c 
	gcc -fPIC -O2 -c -Wall driftscatter.c
	ld -bundle -flat_namespace -undefined suppress \
		-o driftscatter.so driftscatter.o
