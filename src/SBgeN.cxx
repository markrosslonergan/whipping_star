#include "SBgeN.h"

using namespace sbn;

SBgeN::SBgeN(std::string whichxml) : SBNspec(whichxml,-1) {

	//Load all files as per xml
	std::cout<<"SBgeN::SBgeN || Opening all files in xml: "<<whichxml<<std::endl;
	for(auto &fn: multisim_file){
		files.push_back(new TFile(fn.c_str()));
		std::cout<<"SBgeN::SBgeN || Files: "<<fn<<std::endl;
	}
	Nfiles = files.size();
	std::cout<<"SBgeN::SBgeN || Loaded "<<Nfiles<<" in total (multisum_file.size()): "<<multisim_file.size()<<std::endl;

	std::cout<<"SBgeN::SBgeN || Getting all TTrees from Files."<<std::endl;
	for(int i=0; i<multisim_name.size(); i++){
		trees.push_back((TTree*)files.at(i)->Get(multisim_name[i].c_str()) );
	}

	for(auto &t: trees){
	//	nentries.push_back(100000);
		nentries.push_back(t->GetEntries());
	}


	//load all doubles, and int.. as mentioned in xml
	//WARNING!! CURRENTLY ASSUMES THAT EACHFILE HAS SAME BRANCHES. NEED TO FIX THAT!
	vars_i = std::vector<std::vector<int>>(Nfiles   ,    std::vector<int>(branch_names_int.at(0).size(),0));
	vars_d = std::vector<std::vector<double>>(Nfiles   , std::vector<double>(branch_names_double.at(0).size(),0));
	vars_dA = std::vector<std::vector< myarray >>(Nfiles   , std::vector<myarray>(branch_names_double_array.at(0).size()));
	vars_iA = std::vector<std::vector< myarrayInt >>(Nfiles   , std::vector<myarrayInt>(branch_names_int_array.at(0).size()));

	vmapD.resize(Nfiles);
	vmapI.resize(Nfiles);
	vmapDA.resize(Nfiles);
	vmapIA.resize(Nfiles);


	std::cout<<"SBgeN::SBgeN || Setting up Variable Maps."<<std::endl;
	// Create a map for convienance, from variable name in XML to actual variable name. Curently file is still needed... 
	for(int f=0; f<Nfiles; f++){
		for(int i=0; i < branch_names_double.at(f).size(); i++){
			vmapD.at(f)[  branch_names_double.at(f).at(i) ] = &vars_d.at(f).at(i) ;
		}
		for(int i=0; i < branch_names_int.at(f).size(); i++){
			vmapI.at(f)[branch_names_int.at(f).at(i) ] = &vars_i.at(f).at(i) ;
		}
		for(int i=0; i < branch_names_double_array.at(f).size(); i++){
			vmapDA.at(f)[branch_names_double_array.at(f).at(i)] = &vars_dA.at(f).at(i) ;
		}
		for(int i=0; i < branch_names_int_array.at(f).size(); i++){
			vmapIA.at(f)[branch_names_int_array.at(f).at(i)] = &vars_iA.at(f).at(i) ;
		}


	}



	std::cout<<"SBgeN::SBgeN || Starting Branch Address Assignment."<<std::endl;

	for(int i=0; i< Nfiles; i++){
			for(int k=0; k< branch_names_int.at(i).size();k++){
				trees.at(i)->SetBranchAddress(branch_names_int.at(i).at(k).c_str(), &(vars_i.at(i).at(k)));
				std::cout<<"SBgeN::SBgeN || Setting Integer Branch: "<<branch_names_int.at(i).at(k)<<std::endl;
			}
			for(int k=0; k< branch_names_int_array.at(i).size();k++){
				trees.at(i)->SetBranchAddress( branch_names_int_array.at(i).at(k).c_str(), vars_iA.at(i).at(k).data  );
				std::cout<<"SBgeN::SBgeN || Setting Integer Array Branch: "<<branch_names_int_array.at(i).at(k)<<std::endl;
			}


			for(int k=0; k< branch_names_double.at(i).size();k++){
				trees.at(i)->SetBranchAddress(branch_names_double.at(i).at(k).c_str(), &(vars_d.at(i).at(k)));
				std::cout<<"SBgeN::SBgeN || Setting Double Branch: "<<branch_names_double.at(i).at(k)<<std::endl;
			}
		int jj=0;
			for(int k=0; k< branch_names_double_array.at(i).size();k++){
				if( branch_names_double_array_dimension.at(jj).at(k) == 3){
					trees.at(i)->SetBranchAddress(branch_names_double_array.at(i).at(k).c_str(), vars_dA.at(i).at(k).data3);
					std::cout<<"SBgeN::SBgeN || Setting Double Array Branch Dim 3: "<<branch_names_double_array.at(i).at(k)<<std::endl;
				}else{
					trees.at(i)->SetBranchAddress(branch_names_double_array.at(i).at(k).c_str(), vars_dA.at(i).at(k).data);
					std::cout<<"SBgeN::SBgeN || Setting Double Array Branch "<<branch_names_double_array.at(i).at(k)<<std::endl;
				}

				jj++;
			}
	}


	std::cout<<"SBgeN::SBgeN || Setting up detectors SBNdet."<<std::endl;
	for(int i=0; i< detector_names.size(); i++){
		if(detector_bool[i]){
			if(detector_names[i]=="SBND")       detectors.push_back(new SBNdet(DET_SBND,0) 	 );
			if(detector_names[i]=="uBooNE")     detectors.push_back(new SBNdet(DET_UBOONE,1) );
			if(detector_names[i]=="ICARUS")     detectors.push_back(new SBNdet(DET_ICARUS,2) );
		}

	}

	rangen = new TRandom3();



	std::cout<<"SBgeN::SBgeN || Constructor Done."<<std::endl;
}

int SBgeN::doMC(){
	return doMC("NULL");
}

int SBgeN::doMC(std::string nam){

	std::cout<<"SBgeN::doMC || We have "<<Nfiles<<" Files "<<std::endl;
	for(int j=0;j<Nfiles;j++){
		std::cout<<"SBgeN::doMC || Starting file and event loops. On File: "<<multisim_file.at(j)<<" || "<<multisim_name.at(j)<<" #eventts: "<<nentries.at(j)<<std::endl;
		double pot_factor = 1;//pot.at(j)/(pot_scaling.at(j) * (double)nentries.at(j));

		for(int i=0; i< nentries.at(j); i++){
			//if(i%250000==0)std::cout<<"SBgeN::doMC || Event: "<<i<<" of "<<nentries.at(j)<<" POT factor: "<<pot_factor<<" on File "<<j<<std::endl;

			trees.at(j)->GetEntry(i);
			//here we put low level selection criteria, for example nuance interaction 1001 == CCQE, its a virtual bool.
			{
			//if( this->eventSelection(j) ){
				//This is the part where we will every histogram in this Universe
				this->fillHistograms(j, -1, 1);
			}//end of low level selection

		}//end of entry loop	


	double scaled_events = (hist.at(map_hist[multisim_name.at(j)])).GetSumOfWeights();
	std::cout<<"SBgeN::doMC || on\t"<<multisim_name.at(j)<<"\t\tEntries:\t"<<nentries.at(j)<<" Scaled Events:\t"<<scaled_events<<" ratio to entries:\t"<<scaled_events/((double)nentries.at(j))<<"\tFor(x10):\t"<<(scaled_events/((double)nentries.at(j)))/0.1<<" ND: "<<(10*scaled_events/((double)nentries.at(j)))/0.1<<std::endl;

	}//end of file loop 



	std::cout<<"SBgeN::doMC || Finished event loop, running tidyHistograms"<<std::endl;
	this->tidyHistograms();

	std::cout<<"SBgeN::doMC || Finished, writing out "<<std::endl;

	/***************************************************************
	 *		Now some clean-up and Writing
	 * ************************************************************/
	if(nam!="NULL"){
		this->writeSpec(nam);
	}

	std::cout<<"SBgeN::doMC || Done. "<<std::endl;
	return 0;
}


SBgeN::~SBgeN(){
	std::cout<<"SBgeN::~SBgeN || Destructor Starting"<<std::endl;
	for(auto &f: files){
		f->Close();
		delete f;
	}
	for(auto &t: trees){
		//delete t;
	}
	

	std::cout<<"SBgeN::~SBgeN || Destructor Done"<<std::endl;
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

int SBgeN::tidyHistograms(){

	return 0;
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
