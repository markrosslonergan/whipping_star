#include "SBNconfig.h"
#include "TH1D.h"
#include "TFile.h"

SBNconfig::SBNconfig(){
	Ndet = 3;
	Nchan = 2;
	Nmode = 2;
	
	Chan[0] = 7;
	Chan[1] = 2;
	Bins[0] = 11;
	Bins[1] = 19;

	mname.push_back("nu");
	mname.push_back("nubar");

	dname.push_back("SBND");
	dname.push_back("uBooNE");
	dname.push_back("ICARUS");

	cname.push_back("elike");
	cname.push_back("mlike");

	scname.resize(Nchan);
	scname[0].push_back("fulloscnu"); 
	scname[0].push_back("fulloscnuebar"); 
	scname[0].push_back("intrinsic");
	scname[0].push_back("mismuon");
	scname[0].push_back("misphoton");
	scname[0].push_back("dirt");
	scname[0].push_back("cosmic");
	
	scname[1].push_back("intrinsic");
	scname[1].push_back("misncpion");

	Tdet = 0;
	for(int i =0; i< Nchan; i++){
		Tdet += Chan[i]*Bins[i];	
	}		
	Tmode = Tdet*Ndet;
	Tall = Nmode*Tmode;

	std::cout<<"Tdet: "<<Tdet<<" Tmode: "<<Tmode<<" Tall: "<<Tall<<std::endl;


	/*
	char namei[200];
	sprintf(namei,"uBooNE_bkg.root");	
	TFile f2(namei);


	TFile *f = new TFile("test.root","RECREATE");
	f->cd();
	*/




	std::string tempn;
	int indexcount = 0;
	for(int im = 0; im < Nmode; im++){
		for(int id =0; id < Ndet; id++){
			for(int ic = 0; ic < Nchan; ic++){
				for(int sc = 0; sc < Chan[ic]; sc++){
					tempn = mname[im] +"_" +dname[id]+"_"+cname[ic]+"_"+scname[ic][sc];
					std::cout<<tempn<<" from: "<<indexcount<<" to "<<indexcount+Bins[ic]-1<<std::endl;
					
					fullnames.push_back(tempn);
					std::vector<int> tvec = {indexcount, indexcount+Bins[ic]-1}; 

					mapIndex[tempn] = 	tvec;				
					indexcount = indexcount + Bins[ic];
	

					/*				
					if(Bins[ic]==19){
						TH1D * tem =  ((TH1D*)f2.Get("dis_muon_sin")); 
						tem->SetName(tempn.c_str());
						tem->Write();
					}
					else if(Bins[ic]==11){

						TH1D * tem =  ((TH1D*)f2.Get("fullosc_nue_sin")); 
						tem->SetName(tempn.c_str());
						tem->Write();
					}
					*/
				

				}	
			}
		}
	}

/*
	f->Close();
	f2.Close();
*/
	





//std::cout<<"Map test: "<<mapIndex["nubar_SBND_elike_fulloscnu"][1]<<" "<<mapIndex["nubar_SBND_elike_fulloscnu"][0]<<std::endl;
}//end constructor


