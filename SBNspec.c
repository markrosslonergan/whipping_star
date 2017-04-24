#include "SBNspec.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <TFile.h>
#include "params.h"
#include <TH1D.h>

SBNspec::SBNspec(const char * name){

	SBNconfig();

	char namei[200];
	sprintf(namei,"%s_bkg.root",name);	
	TFile f(namei);



	for(auto f: filenames){ 	//Loop over all filenames that should be there, and load up the histograms.
		hist.push_back(*((TH1D*)f.Get(f.c_str()))); 
	}


	std::cout<<hist[8].GetBinContent(2)<<std::endl;


f.Close();



}

