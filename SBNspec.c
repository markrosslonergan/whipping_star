#include "SBNspec.h"
#include <ctime>
#include <TFile.h>
#include "params.h"
#include <TH1D.h>
#include <TRandom3.h>

SBNspec::SBNspec(const char * name, std::string whichxml) : SBNconfig(whichxml) {


	char namei[200];
	sprintf(namei,"%s.root",name);	
	TFile f(namei);



	for(auto fn: fullnames){ 	//Loop over all filenames that should be there, and load up the histograms.
		//std::cout<<"Attempting to load: "<<fn.c_str()<<" from: "<<namei<<std::endl;
		hist.push_back(*((TH1D*)f.Get(fn.c_str()))); 
	}


	f.Close();
}//end constructor


int SBNspec::Add(SBNspec *in){

	if(xmlname != in->xmlname){ std::cout<<"ERROR: SBNspec::Add, trying to add differently configured SBNspecs!"<<std::endl; exit(EXIT_FAILURE);}

	for(int i=0; i< hist.size(); i++){
		hist[i].Add( &(in->hist[i]));
	}	


	return 1;
}


int SBNspec::randomScale(){
	TRandom3 *rangen    = new TRandom3(0);

	for(auto& h: hist){
			h.Scale(rangen->Uniform(0,2));

	}
return 1;
}

int SBNspec::ScaleAll(double sc){
	for(auto& h: hist){
			h.Scale(sc);

	}
return 1;
}

int SBNspec::Scale(std::string name, double val){
	for(auto& h: hist){
		std::string test = h.GetName();
		
			if(test.find(name)!=std::string::npos ){
				//std::cout<<name<<". found in: "<<test<<" at "<<test.find(name)<<std::endl;
				h.Scale(val);
			}

	}
return 1;
}

int SBNspec::NormAll(double n){

	for(auto& h: hist) {
		h.Scale(n/h.GetSumOfWeights());
	}
	return 1;
}

int SBNspec::Norm(std::string name, double val){
	for(auto& h: hist){
		std::string test = h.GetName();
		
			if(test.find(name)!=std::string::npos ){
				//std::cout<<name<<". found in: "<<test<<" at "<<test.find(name)<<std::endl;
				h.Scale(val/h.GetSumOfWeights());
			}

	}
return 1;
}

int SBNspec::calcFullVector(){
	fullVec.clear();

	for(auto& h: hist){
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
	//This needs to be confirmed and checked. Looks good, mark 24th april

	for(int im = 0; im < Nmode; im++){
		for(int id =0; id < Ndet; id++){
			int edge = id*Tdet + Tmode*im; // This is the starting index for this detector


			for(int ic = 0; ic < Nchan; ic++){
				int corner=edge;
				for(int j=0; j< Bins[ic]; j++){
					double tempval=0;

					for(int sc = 0; sc < Chan[ic]; sc++){
						tempval += fullVec[j+sc*Bins[ic]+corner];
						edge +=1;	//when your done with a channel, add on every bin you just summed
					}
					compVec.push_back(tempval);
				}
				
				
				
									


			}
		}
	}
return 1;
}

int SBNspec::printFullVec(){
	for(double d: fullVec){
		std::cout<<d<<" ";
	}
	std::cout<<std::endl;
return 1;
}

int SBNspec::printCompVec(){ 
	for(double d: compVec){
		std::cout<<d<<" ";
	}
	std::cout<<std::endl;
return 1;
}






