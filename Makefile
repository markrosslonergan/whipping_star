mc:	
	g++ -g -std=c++11 -c SBNcovar.c -o SBNcovar.o  -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNconfig.c -o SBNconfig.o  -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNfit3pN.c -o SBNfit3pN.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNspec.c -o SBNspec.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c prob.c  -o prob.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNosc.c  -o SBNosc.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNchi.c -o SBNchi.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNfit.c -o SBNfit.o  -I $$ROOTSYS/include
	g++ -g -std=c++11 -c testmc.cxx -o testmc.o   -I $$ROOTSYS/include 
	g++ -g -std=c++11  -Wall -fPIC  -Wno-deprecated-declaratios tinyxml/* SBNfit3pN.o testmc.o prob.o SBNfit.o SBNconfig.o SBNspec.o SBNosc.o SBNchi.o SBNcovar.o -o testmc -lgsl -lgslcblas -L/home/mark/programs/root/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMathMore -lThread -lMultiProc -pthread -lm -ldl -rdynamic 

all:	
	g++ -g -std=c++11 -c SBNcovar.c -o SBNcovar.o  -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNconfig.c -o SBNconfig.o  -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNfit3pN.c -o SBNfit3pN.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNspec.c -o SBNspec.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c prob.c  -o prob.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNosc.c  -o SBNosc.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNchi.c -o SBNchi.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c SBNfit.c -o SBNfit.o  -I $$ROOTSYS/include
	g++ -g -std=c++11 -c examples.cxx -o examples.o   -I $$ROOTSYS/include 
	g++ -g -std=c++11  -Wall -fPIC  -Wno-deprecated-declaratios tinyxml/* SBNfit3pN.o examples.o prob.o SBNfit.o SBNconfig.o SBNspec.o SBNosc.o SBNchi.o SBNcovar.o -o examples -lgsl -lgslcblas -L/home/mark/programs/root/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMathMore -lThread -lMultiProc -pthread -lm -ldl -rdynamic 

