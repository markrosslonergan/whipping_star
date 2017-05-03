all:	
	g++ -g -std=c++11 -c SBNconfig.c -o SBNconfig.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNspec.c -o SBNspec.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c prob.c  -o prob.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNosc.c  -o SBNosc.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNchi.c -o SBNchi.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c test.cxx -o test.o   -I $$ROOTSYS/include 
	g++ -g -std=c++11  -Wall -fPIC  -Wno-deprecated-declaratios tinyxml/* test.o prob.o SBNconfig.o SBNspec.o SBNosc.o SBNchi.o -o test -lgsl -lgslcblas -L/home/mark/programs/root/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic 

