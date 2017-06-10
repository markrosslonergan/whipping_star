#ifndef SBGEN_H_
#define SBGEN_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNconfig.h"
#include "SBNspec.h"
#include <TH1D.h>
#include <string>
#include <TF1.h>
#include <THStack.h>
#include <TTree.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <array>
#include <map>

#include <ctime>
#include <TFile.h>
#include "params.h"
#include <TRandom3.h>

namespace sbn{

struct myarray{

	double data[100];
};

//This is the basic class that holds all spectral information in whatever reco or true variable you have decided you want in the xml files.
// Inherits from SBNconfig as thats how its all configured/kept equal! :
class SBgeN : public SBNspec{
	
	public:

	SBgeN(std::string xmlname);
	~SBgeN();

	std::vector<std::vector<int> > vars_i;
	std::vector<std::vector<double> > vars_d;
	std::vector<std::vector< myarray> >  vars_dA;

	std::vector<TFile *> files;	
	std::vector<TTree *> trees;	

	std::vector<int> nentries;
	int Nfiles;

	virtual bool eventSelection(int file);
	virtual int fillHistograms(int file, int uni, double wei);
	virtual int tidyHistograms();


	std::vector<std::map<std::string, double* >> vmapD;
	std::vector<std::map<std::string, int* >> vmapI;
	std::vector<std::map<std::string, myarray* >> vmapDA;

	int doMC();

};



}

#endif
