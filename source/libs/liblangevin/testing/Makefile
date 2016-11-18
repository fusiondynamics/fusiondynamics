all: turfu_langevin 

turfu_langevin: langevin.o
	g++ langevin.o turfu_langevin.cxx -o turfu_langevin

langevin.o:
	g++ -c ../langevin/langevin.cpp -o langevin.o

clean:
	rm -rf *.o

mrproper: clean
	rm -rf turfu_langevin
