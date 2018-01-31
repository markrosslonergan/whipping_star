#include "SBNosc.h"
using namespace sbn;

SBNosc::SBNosc(const char * name, std::string whichxml) : SBNspec(name, whichxml) {
	workingModel.zero();
	mStepSize = 0.04;
	which_mode = BOTH_ONLY;
}

SBNosc::SBNosc(const char * name, std::string whichxml, neutrinoModel in) : SBNosc(name, whichxml) {
	load_model(in);

}

int SBNosc::load_model(neutrinoModel in){
	workingModel = in;
	calcMassSplittings();
	return 0;
}

/*************************************************
 * for a given workingModel. 
 * Calculate how many mass splittings, and which "type" it is,
 *
 * **********************************************/

int SBNosc::calcMassSplittings(){
	mass_splittings.clear();

	double fix41=round(log10((workingModel.dm41Sq))/mStepSize)*mStepSize;
	double fix51=round(log10((workingModel.dm51Sq))/mStepSize)*mStepSize;
	double fix61=round(log10((workingModel.dm61Sq))/mStepSize)*mStepSize;
	
	double round54 = round(log10(fabs(workingModel.dm54Sq))/mStepSize)*mStepSize;
	double round64 = round(log10(fabs(workingModel.dm64Sq))/mStepSize)*mStepSize;
	double round65 = round(log10(fabs(workingModel.dm65Sq))/mStepSize)*mStepSize;

	if(workingModel.numsterile == 1)
	{
		mass_splittings.push_back( std::make_pair(fix41,41 ));
	}
       	else if(workingModel.numsterile ==2)
	{
		mass_splittings.push_back( std::make_pair (fix41,41));
		mass_splittings.push_back( std::make_pair (fix51,51));
		mass_splittings.push_back( std::make_pair (round54,54));
	}	
	else if(workingModel.numsterile ==3)
	{
		mass_splittings.push_back( std::make_pair (fix41,41));
		mass_splittings.push_back( std::make_pair (fix51,51));
		mass_splittings.push_back( std::make_pair (fix61,61));
	
		mass_splittings.push_back( std::make_pair (round54,54));
		mass_splittings.push_back( std::make_pair (round64,64));
		mass_splittings.push_back( std::make_pair (round65,65));
	}

	


	return 0;
}


int SBNosc::OscillateThis(){
		this->calcFullVector();
		this->compressVector();

	calcMassSplittings();

	for(auto ms: mass_splittings){
			char namei[200];
			
			sprintf(namei, "%sprecomp/SBN_SIN_%2.2f", data_path.c_str(), ms.first );
			SBNspec single_frequency(namei , xmlname , false);
			
			sprintf(namei, "%sprecomp/SBN_SINSQ_%2.2f", data_path.c_str(),ms.first );
			SBNspec single_frequency_square(namei , xmlname ,false);


			if(has_been_scaled){
				single_frequency.Scale(scale_hist_name, scale_hist_val);
				single_frequency_square.Scale(scale_hist_name, scale_hist_val);
			}

			double prob_mumu, prob_ee, prob_mue, prob_mue_sq, prob_muebar, prob_muebar_sq;

			int which_dm = ms.second;

			switch (which_mode)
				{
					case APP_ONLY: //Strictly nu_e app only
						prob_mumu =0;
						prob_ee   =0;
						prob_mue = workingModel.oscAmp(2,1,which_dm,1);
						prob_mue_sq = workingModel.oscAmp(2,1,which_dm,2);
						prob_muebar = workingModel.oscAmp(-2,-1,which_dm,1);	
						prob_muebar_sq = workingModel.oscAmp(-2,-1,which_dm,2);				
						break;
					case DIS_ONLY: //Strictly nu_mu dis only
						prob_mumu = workingModel.oscAmp(2,2,which_dm,2);
						prob_ee = 0;
						prob_mue = 0;
						prob_mue_sq =0;
						prob_muebar =0;
						prob_muebar_sq =0;				
						break;
					case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
						prob_mumu = workingModel.oscAmp(2,2,which_dm,2);
						prob_ee = workingModel.oscAmp(1,1,which_dm,2);
						prob_mue = workingModel.oscAmp(2,1,which_dm,1);
						prob_mue_sq = workingModel.oscAmp(2,1,which_dm,2);
						prob_muebar = workingModel.oscAmp(-2,-1,which_dm,1);	
						prob_muebar_sq = workingModel.oscAmp(-2,-1,which_dm,2);			
						break;
					case WIERD_ONLY: // A strange version where nu_e can appear but not disapear 
						prob_mumu =workingModel.oscAmp(2,2,which_dm,2);
						prob_ee = 0.0;
						prob_mue = workingModel.oscAmp(2,1,which_dm,1);
						prob_mue_sq = workingModel.oscAmp(2,1,which_dm,2);
						prob_muebar = workingModel.oscAmp(-2,-1,which_dm,1);	
						prob_muebar_sq = workingModel.oscAmp(-2,-1,which_dm,2);			
						break;
					case DISE_ONLY: // A strange version where nu_e can appear but not disapear 
						prob_mumu = 0.0;
						prob_ee = workingModel.oscAmp(1,1,which_dm,2);
						prob_mue = 0;
						prob_mue_sq = 0;
						prob_muebar = 0;	
						prob_muebar_sq = 0;			
						break;

				}



			//std::cout<<"mm: "<<prob_mumu<<" ee: "<<prob_ee<<" mue: "<<prob_mue<<" mueSQ: "<<prob_mue_sq<<" mubar: "<<prob_muebar<<" muebarSQ: "<<prob_muebar_sq<<std::endl;
			
			
			for(int i=0; i<num_channels; i++){
				for(int j=0; j<num_subchannels.at(i); j++){
					int osc_pattern = subchannel_osc_patterns.at(i).at(j);
					double osc_amp = 0;
					double osc_amp_sq = 0;
					switch(osc_pattern){
						case 11:
							osc_amp_sq = prob_ee;
							break;
						case -11:
							osc_amp_sq = prob_ee;
							break;
						case 22:
							osc_amp_sq = prob_mumu;
							break;
						case -22:
							osc_amp_sq = prob_mumu;
							break;
						case 21:
							osc_amp = prob_mue;
							osc_amp_sq = prob_mue_sq;
							break;
						case -21:
							osc_amp = prob_muebar;
							osc_amp_sq = prob_muebar_sq;
							break;
						case 0:
						default:
							break;
					}
				
					single_frequency.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp );
					single_frequency_square.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp_sq );
				}
			}



			/*single_frequency.Scale("elike_fulloscnue", prob_mue);	
			single_frequency.Scale("elike_fulloscbarnue", prob_muebar);	
			single_frequency.Scale("elike_intrinsic", prob_ee);
			single_frequency.Scale("elike_mismuon", prob_mumu);	
			single_frequency.Scale("elike_misphoton",0.0);
			single_frequency.Scale("elike_dirt",0.0);
			single_frequency.Scale("elike_cosmic",0.0);
			single_frequency.Scale("mlike_intrinsic", prob_mumu);
			single_frequency.Scale("mlike_misncpion",0.0);

			single_frequency_square.Scale("elike_fulloscnue", prob_mue_sq); 
			single_frequency_square.Scale("elike_fulloscbarnue", prob_muebar_sq);
			single_frequency_square.Scale("elike_intrinsic", 0.0);
			single_frequency_square.Scale("elike_mismuon", 0.0);	
			single_frequency_square.Scale("elike_misphoton",0.0);
			single_frequency_square.Scale("elike_dirt",0.0);
			single_frequency_square.Scale("elike_cosmic",0.0);
			single_frequency_square.Scale("mlike_intrinsic", 0.0);
			single_frequency_square.Scale("mlike_misncpion",0.0);
			*/

			this->Add(&single_frequency);
			this->Add(&single_frequency_square);	


	}//Done looping over

	this->calcFullVector();
	this->compressVector();


	return 0;
};

std::vector<double> SBNosc::Oscillate(double scale){

	std::vector<double> tmp = this->Oscillate();
	for(auto & v: tmp){
		v=v*scale;
	}

	return tmp;
}


std::vector<double> SBNosc::Oscillate(){

		this->calcFullVector();
		this->compressVector();

	std::vector<double > temp = compVec;


	for(auto ms: mass_splittings){
			char namei[200];
			
			sprintf(namei, "%sprecomp/SBN_SIN_%2.2f",data_path.c_str(), ms.first );
			SBNspec single_frequency(namei , xmlname, false);		
	
			sprintf(namei, "%sprecomp/SBN_SINSQ_%2.2f",data_path.c_str(), ms.first );
			SBNspec single_frequency_square(namei , xmlname, false);

			if(has_been_scaled){
				single_frequency.Scale(scale_hist_name, scale_hist_val);
				single_frequency_square.Scale(scale_hist_name, scale_hist_val);
			}


			double prob_mumu, prob_ee, prob_mue, prob_mue_sq, prob_muebar, prob_muebar_sq;

			int which_dm = ms.second;

			switch (which_mode)
				{
					case APP_ONLY: //Strictly nu_e app only
						prob_mumu =0;
						prob_ee   =0;
						prob_mue = workingModel.oscAmp(2,1,which_dm,1);
						prob_mue_sq = workingModel.oscAmp(2,1,which_dm,2);
						prob_muebar = workingModel.oscAmp(-2,-1,which_dm,1);	
						prob_muebar_sq = workingModel.oscAmp(-2,-1,which_dm,2);				
						break;
					case DIS_ONLY: //Strictly nu_mu dis only
						prob_mumu = workingModel.oscAmp(2,2,which_dm,2);
						prob_ee = 0;
						prob_mue = 0;
						prob_mue_sq =0;
						prob_muebar =0;
						prob_muebar_sq =0;				
						break;
					case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
						prob_mumu = workingModel.oscAmp(2,2,which_dm,2);
						prob_ee = workingModel.oscAmp(1,1,which_dm,2);
						prob_mue = workingModel.oscAmp(2,1,which_dm,1);
						prob_mue_sq = workingModel.oscAmp(2,1,which_dm,2);
						prob_muebar = workingModel.oscAmp(-2,-1,which_dm,1);	
						prob_muebar_sq = workingModel.oscAmp(-2,-1,which_dm,2);			
						break;
					case WIERD_ONLY: // A strange version where nu_e can appear but not disapear 
						prob_mumu =workingModel.oscAmp(2,2,which_dm,2);
						prob_ee = 0.0;
						prob_mue = workingModel.oscAmp(2,1,which_dm,1);
						prob_mue_sq = workingModel.oscAmp(2,1,which_dm,2);
						prob_muebar = workingModel.oscAmp(-2,-1,which_dm,1);	
						prob_muebar_sq = workingModel.oscAmp(-2,-1,which_dm,2);			
						break;
					case DISE_ONLY: // A strange version where nu_e can appear but not disapear 
						prob_mumu = 0.0;
						prob_ee = workingModel.oscAmp(1,1,which_dm,2);
						prob_mue = 0;
						prob_mue_sq = 0;
						prob_muebar = 0;	
						prob_muebar_sq = 0;			
						break;

				}

	
			for(int i=0; i<num_channels; i++){
				for(int j=0; j<num_subchannels.at(i); j++){
					int osc_pattern = subchannel_osc_patterns.at(i).at(j);
					double osc_amp = 0;
					double osc_amp_sq = 0;
					switch(osc_pattern){
						case 11:
							osc_amp_sq = prob_ee;
							break;
						case -11:
							osc_amp_sq = prob_ee;
							break;
						case 22:
							osc_amp_sq = prob_mumu;
							break;
						case -22:
							osc_amp_sq = prob_mumu;
							break;
						case 21:
							osc_amp = prob_mue;
							osc_amp_sq = prob_mue_sq;
							break;
						case -21:
							osc_amp = prob_muebar;
							osc_amp_sq = prob_muebar_sq;
							break;
						case 0:
						default:
							break;
					}
				
					single_frequency.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp);
					single_frequency_square.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp_sq );
					
				}
			}


		//	std::cout<<"mm: "<<prob_mumu<<" ee: "<<prob_ee<<" mue: "<<prob_mue<<" mueSQ: "<<prob_mue_sq<<" mubar: "<<prob_muebar<<" muebarSQ: "<<prob_muebar_sq<<std::endl;
			/*
			single_frequency.Scale("elike_fulloscnue", prob_mue);	
			single_frequency.Scale("elike_fulloscbarnue", prob_muebar);	
			single_frequency.Scale("elike_intrinsic", prob_ee);
			single_frequency.Scale("elike_mismuon", prob_mumu);	
			single_frequency.Scale("elike_misphoton",0.0);
			single_frequency.Scale("elike_dirt",0.0);
			single_frequency.Scale("elike_cosmic",0.0);
			single_frequency.Scale("mlike_intrinsic", prob_mumu);
			single_frequency.Scale("mlike_misncpion",0.0);

			single_frequency_square.Scale("elike_fulloscnue", prob_mue_sq); 
			single_frequency_square.Scale("elike_fulloscbarnue", prob_muebar_sq);
			single_frequency_square.Scale("elike_intrinsic", 0.0);
			single_frequency_square.Scale("elike_mismuon", 0.0);	
			single_frequency_square.Scale("elike_misphoton",0.0);
			single_frequency_square.Scale("elike_dirt",0.0);
			single_frequency_square.Scale("elike_cosmic",0.0);
			single_frequency_square.Scale("mlike_intrinsic", 0.0);
			single_frequency_square.Scale("mlike_misncpion",0.0);
			*/


			single_frequency.calcFullVector();
			single_frequency.compressVector();
			

			single_frequency_square.calcFullVector();
			single_frequency_square.compressVector();
	
			for(int i=0;i<temp.size(); i++){
				temp[i] += single_frequency.compVec[i];
				temp[i] += single_frequency_square.compVec[i];
			}

	}//Done looping over

	return temp;
};


/*************************************************
 * for a given workingModel. 
 *  precomute the sin and sin^2 for this,
 *
 * **********************************************/





int SBNosc::setMode(int in){
	which_mode = in;

return in;
}

void SBNosc::setAppMode(){
	setMode(APP_ONLY);
}

void SBNosc::setDisMode(){
	setMode(DIS_ONLY);
}

void SBNosc::setBothMode(){
	setMode(BOTH_ONLY);
}

void SBNosc::setWierdMode(){
	setMode(WIERD_ONLY);
}

void SBNosc::setDisEMode(){
	setMode(DISE_ONLY);
}

