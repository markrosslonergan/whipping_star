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
	sprintf(namei,"%s.root",name);	
	TFile f(namei);



	for(auto fn: fullnames){ 	//Loop over all filenames that should be there, and load up the histograms.
		std::cout<<"Attempting to load: "<<fn.c_str()<<" from: "<<namei<<std::endl;
		hist.push_back(*((TH1D*)f.Get(fn.c_str()))); 
	}
	std::cout<<"Finished loading all spectra histograms"<<std::endl;






f.Close();
}//end constructor

void SBNspec::randomScale(){

	for(auto h: hist){
	

	}

}


int SBNspec::calcFullVector(){
	fullVec.clear();

	for(auto h: hist){
		//std::cout<<"Hist size: "<<h.GetSize()-2<<std::endl;
		for(int i = 1; i <= h.GetSize()-2; i++){
			//std::cout<<h.GetBinContent(i)<<" ";
			fullVec.push_back(h.GetBinContent(i));
		}	
	}

return 1;
}

int SBNspec::compressVector(){
	compVec.clear();


	for(int im = 0; im < Nmode; im++){
		for(int id =0; id < Ndet; id++){
			int edge = id*Tdet + Tmode*im;
			for(int ic = 0; ic < Nchan; ic++){
				
				for(int sc = 0; sc < Chan[ic]; sc++){
				

				}	
			}
		}
	}
return 1;
}









