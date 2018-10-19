cc = g++
run: main.o  header.h
	$(cc)  main.o -o  run

main.o: main.cpp header.h
	$(cc) -c main.cpp -O2 -std=c++11

clean:
	\rm *.*~ *~ *.o

cleanrun:
	\rm run
