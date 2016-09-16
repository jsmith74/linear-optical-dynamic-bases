all: LOCircuitDesign

LOCircuitDesign: main.o LOTransform.o AncillaAugment.o BFGS_Optimization.o PUA.o MeritFunction.o ChaseAlgorithm.o
	g++ -O3 -o LOCircuitDesign LOTransform.o main.o AncillaAugment.o BFGS_Optimization.o PUA.o MeritFunction.o ChaseAlgorithm.o

ChaseAlgorithm.o: ChaseAlgorithm.cpp
	g++ -O3 -c ChaseAlgorithm.cpp

main.o: main.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN main.cpp

MeritFunction.o: MeritFunction.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN MeritFunction.cpp

PUA.o: PUA.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN PUA.cpp

AncillaAugment.o: AncillaAugment.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN AncillaAugment.cpp

BFGS_Optimization.o: BFGS_Optimization.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN BFGS_Optimization.cpp

LOTransform.o: LOTransform.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN LOTransform.cpp

clean:
	rm *o LOCircuitDesign *.dat *.out
