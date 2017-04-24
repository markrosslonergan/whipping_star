#ifndef SBNSPEC_H_
#define SBNSPEC_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNconfig.h"
#include <TH1D.h>

class SBNspec :public SBNconfig{
	
	//Have a creation version, that creates an appropiate root one?
	public:
	std::vector<TH1D > hist;	
//	std::vector< std::vector< std::vector<std::vector < TH1D *> >>> mhists;
	std::vector<double > fullVec;
	std::vector<double > compVec;

	SBNspec(const char *); //Load in config file

	void randomScale();

	// Todo: most important
	int calcFullVector();
	int compressVector();


};




#endif
