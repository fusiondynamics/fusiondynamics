all: test_langevin

test_langevin: langevin.o
	g++ langevin.o test_langevin.cxx -o test_langevin

langevin.o:
	g++ -c ../langevin/langevin.cpp -o langevin.o

clean:
	rm -rf *.o

mrproper: clean
	rm -rf test_langevin
