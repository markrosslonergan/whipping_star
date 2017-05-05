#ifndef SBNSPEC_H_
#define SBNSPEC_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNconfig.h"
#include <TH1D.h>
#include <string>
#include <TF1.h>

class SBNspec : public SBNconfig{
	
	//Have a creation version, that creates an appropiate root one?
	public:
	std::vector<TH1D > hist;	
//	std::vector< std::vector< std::vector<std::vector < TH1D *> >>> mhists;
	std::vector<double > fullVec;
	std::vector<double > compVec;


	SBNspec(const char *, std::string); //Load in config file
	SBNspec() {};

	int randomScale(); //mainly a debugging function, just randomly scales each hist by 0-2
	int poissonScale(); //Scales every histogram by a poissonian random number

	int ScaleAll(double);
	int Scale(std::string name, double val);
	int Scale(std::string name, TF1 *);
	int NormAll(double);
	int Norm(std::string name, double val);


	//Done: most important
	int calcFullVector();
	int compressVector();

	int printFullVec();
	int printCompVec();

	int Add(SBNspec*);

};



#endif
