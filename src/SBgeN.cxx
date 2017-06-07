#include "SBgeN.h"

using namespace sbn;

SBgeN::SBgeN(std::string whichxml) : SBNspec(whichxml,-1) {

	//Load all files as per xml
	std::vector<TFile *> files;	
	std::vector<TTree *> trees;	

	std::cout<<"SBgeN::SBgeN || Opening all files in xml: "<<whichxml<<std::endl;
	for(auto &fn: multisim_file){
		files.push_back(new TFile(fn.c_str()));
	}

	std::cout<<"SBgeN::SBgeN || Getting all TTrees from Files."<<std::endl;
	for(int i=0; i<multisim_name.size(); i++){
		trees.push_back((TTree*)files.at(i)->Get(multisim_name[i].c_str()) );
	}

	std::vector<int> nentries;
	for(auto &t: trees){
		nentries.push_back(5000);
		//nentries.push_back(t->GetEntries());
	}
	int Nfiles = files.size();

	//load all doubles, and int.. as mentioned in xml
	vars_i = std::vector<std::vector<int>>(Nfiles   ,    std::vector<int>(branch_names_int.at(0).size(),0));
	vars_d = std::vector<std::vector<double>>(Nfiles   , std::vector<double>(branch_names_double.at(0).size(),0.0));

	std::cout<<"SBgeN::SBgeN || Starting Branch Address Assignment."<<std::endl;
	for(int i=0; i< Nfiles; i++){
		for(auto &bfni: branch_names_int){
			for(int k=0; k< bfni.size();k++){
				trees.at(i)->SetBranchAddress(bfni[k].c_str(), &(vars_i.at(i).at(k)));
			}
		}
		for(auto &bfnd: branch_names_double){
			for(int k=0; k< bfnd.size();k++){
				trees.at(i)->SetBranchAddress(bfnd[k].c_str(), &(vars_d.at(i).at(k)));
			}
		}
	}

	std::cout<<"SBgeN::SBgeN || Starting file and event loops."<<std::endl;
	for(int j=0;j<Nfiles;j++){
		double pot_factor = pot.at(j)/(pot_scaling.at(j) * (double)nentries.at(j));

		for(int i=0; i< nentries.at(j); i++){
			if(i%2500==0)std::cout<<"SBgeN::SBgeN || Event: "<<i<<" of "<<nentries[j]<<" from File: "<<multisim_file[j]<<" POT factor: "<<pot_factor<<std::endl;

			trees.at(j)->GetEntry(i);
			//here we put low level selection criteria, for example nuance interaction 1001 == CCQE, its a virtual bool.
			if( this->eventSelection(j) ){
				//This is the part where we will every histogram in this Universe
				this->fillHistograms(j, -1, 1);
			}//end of low level selection
	
		}//end of entry loop
	}//end of file loop 


	/***************************************************************
	 *		Now some clean-up and Writing
	 * ************************************************************/

	this->writeOut("gen.root");

	for(auto &f: files){
		f->Close();
	}

}

/***************************************************************
 *		Some virtual functions for selection and histogram filling
 * ************************************************************/

bool SBgeN::eventSelection(int which_file){
	//from here have access to vars_i  and vars_d  to make a selection

	bool ans = false;
	if(which_file==0){
		if(vars_i.at(which_file)[3] == 11 || vars_i.at(which_file)[3] == -11  ){
			ans = true;
		}
	} else if (which_file==1){
		if(vars_i.at(which_file)[3] == 13 || vars_i.at(which_file)[3] == -13  ){
			ans = true;
		}


	}

	return ans;
}

int SBgeN::fillHistograms(int file, int uni, double wei){
	//Fill the histograms
	//
	/*	TRandom3 rangen(0);
		double sigma = 0;
		if(file==0){
		sigma=;
		}else if(file==1){
		sigma=;
		}	
		double en = rangen.Gaus( vars_d.at(file), sigma );
		*/	
	double en = vars_d.at(file)[0];
	hist.at(file).Fill(en, wei);

	return 0;

}
