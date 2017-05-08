#ifndef SBNSPEC_H_
#define SBNSPEC_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNconfig.h"
#include <TH1D.h>
#include <string>
#include <TF1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
//#include <TROOT.h>

#include <ctime>
#include <TFile.h>
#include "params.h"
#include <TRandom3.h>
namespace sbn{



class SBNspec : public SBNconfig{
	
	public:
	std::vector<TH1D > hist;	
	std::vector<double > fullVec;
	std::vector<double > compVec;

	SBNspec(const char *, std::string); //Load in root file and config file
	SBNspec(std::string, int); //Load in config file, create EMPTY hists, with optional numbering 
	SBNspec() {};

	int randomScale(); //mainly a debugging function, just randomly scales each hist by 0-2
	int poissonScale(); //Scales every histogram by a poissonian random number

	int ScaleAll(double);
	int Scale(std::string name, double val);
	int Scale(std::string name, TF1 *);
	int NormAll(double);
	int Norm(std::string name, double val);


	int calcFullVector();
	int compressVector();

	int printFullVec();
	int printCompVec();
	int writeOut(std::string);

	int Add(SBNspec*);

};



}

#endif
