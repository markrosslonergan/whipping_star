#ifndef SBNSPEC_H_
#define SBNSPEC_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNconfig.h"
#include <TH1D.h>

class SBNspec :public SBNconfig{
	
	public:
	std::vector<TH1D > hist;	
	
	SBNspec(const char *); //Load in config file



};




#endif
