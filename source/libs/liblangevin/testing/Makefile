all: new_test_langevin

new_test_langevin: langevin.o
	g++ langevin.o new_test_langevin.cxx -o new_test_langevin

langevin.o:
	g++ -c ../langevin/langevin.cpp -o langevin.o

clean:
	rm -rf *.o

mrproper: clean
	rm -rf new_test_langevin
