#include "SBNconfig.h"
#include "TH1D.h"
#include "TFile.h"

SBNconfig::SBNconfig(){

	CorrMatRoot ="rootfiles/covariance_matrices_xcheck_690x690.root";
	CorrMatName ="TMatrixT<float>;7";
	

	Ndet = 3;
	Nchan = 2;
	Nmode = 2;
	
	Chan[0] = 7;
	Chan[1] = 2;

	Bins[0] = 11;
	Bins[1] = 19;

	mname.push_back("nu");
	mname.push_back("nubar");

	mBool.push_back(1);
	mBool.push_back(1);

	dname.push_back("SBND");
	dname.push_back("uBooNE");
	dname.push_back("ICARUS");

	dBool.push_back(1);
	dBool.push_back(1);
	dBool.push_back(1);

	cname.push_back("elike");
	cname.push_back("mlike");

	cBool.push_back(1);
	cBool.push_back(1);

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

	scBool.resize(Nchan);
	for(int i=0;i<7;i++){ scBool[0].push_back(1);}
	scBool[1].push_back(1);
	scBool[1].push_back(1);

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
				
			//std::cout<<tempn<<" from: "<<indexcount<<" to "<<indexcount+Bins[ic]-1<<" mBool: "<<mBool[im]<<" dBool: "<<dBool[id]<<" cBool: "<<cBool[ic]<<std::endl;
					// This is where you choose NOT to use some fields	
					if(mBool[im] && dBool[id] && cBool[ic] && scBool[ic][sc]){					
						fullnames.push_back(tempn);

						for(int k = indexcount; k < indexcount+Bins[ic]; k++){
								useBins.push_back(k);
							
						}
					}
					
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

	//For here on down everything is derivable, above is just until I actually get config working.
	Nmode = 0;
	for(bool i:mBool){	if(i) Nmode++;	}

	Ndet = 0;
	for(bool i:dBool){	if(i) Ndet++;	}
	
	Nchan = 0;
	for(bool i:cBool){	if(i) Nchan++;	}

	for(int i=0; i< Nchan; i++){
		Chan[i] = 0;
		for(bool j: scBool[i]){ if(j) Chan[i]++;}
	}



	Tdet = 0;
	TdetComp = 0;
	for(int i =0; i< Nchan; i++){
		Tdet += Chan[i]*Bins[i];
		TdetComp += Bins[i];	
	}		
	Tmode = Tdet*Ndet;
	TmodeComp = TdetComp*Ndet;
	Tall = Nmode*Tmode;
	TallComp = Nmode*TmodeComp;


	std::cout<<"Ndet: "<<Ndet<<" Nmode: "<<Nmode<<" Nchan: "<<Nchan<<" Nsc[0]: "<<Chan[0]<<" Nsc[1]: "<<Chan[1]<<std::endl;
	std::cout<<"Tdet: "<<Tdet<<" Tmode: "<<Tmode<<" Tall: "<<Tall<<" Tcomp: "<<TallComp<<std::endl;




/*
	f->Close();
	f2.Close();
*/
	





//std::cout<<"Map test: "<<mapIndex["nubar_SBND_elike_fulloscnu"][1]<<" "<<mapIndex["nubar_SBND_elike_fulloscnu"][0]<<std::endl;
}//end constructor


