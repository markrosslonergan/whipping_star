#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNdet.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBgeN.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovar.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;


struct myphoton{
	TLorentzVector lorentz;
	int isPion;
	double convlength;
	double *convpoint;

};




class SBgeNicarus : public SBgeN{
	public:
		SBgeNicarus(std::string xmlname):SBgeN(xmlname){};
		bool eventSelection(int file);
		int fillHistograms(int file, int uni, double wei);
		int tidyHistograms();
		~SBgeNicarus(){};

		double vertexEnergy(int file);


};


bool SBgeNicarus::eventSelection(int file){

	return true;
}

int SBgeNicarus::fillHistograms(int file, int uni, double wei){
	double vertex_pos[3] = {0,0,0};

	//file 0 is general one
	if(file==0){

		double Enu_true = *vmapD[file]["Enu"];//   vars_d.at(file).at(0);
		double El_true =  *vmapD[file]["El"];

		int CC = *vmapI[file]["CC"];
		int NC = *vmapI[file]["NC"];
		int Npi0dph = *vmapI[file]["Npi0dph"];
		int Nph = *vmapI[file]["Nph"];

		double *Eph =  vmapDA[file]["Eph"]->data; 
		double *Epi0dph =  vmapDA[file]["Epi0dph"]->data; 

		double (*pph)[3] = vmapDA[file]["pph"]->data3; 
		double (*ppi0dph)[3] =  vmapDA[file]["ppi0dph"]->data3; 

		/********************************************************************
		 *			Starting main cutflow
		 * *****************************************************************/


		//	detectors[0]->random_pos(rangen,vertex_pos); //Assign a random position for vertex in detector 
		double osclen = detectors[0]->osc_length(rangen);


		//Is there a visible vertex and how much energy is there!
		double Enu_reco = 0;
		bool vis_vertex = false;
		double vertex_energy = 	this->vertexEnergy(file);

		//Check if we actually have a "visibe vertex"
		if(vertex_energy > vertex_thresh)
		{
			vis_vertex = true;
		} else 
		{
			vis_vertex = false;
		}


		/************************************************************************************************
		 *				CC Intrinsic Nu_e 
		 * **********************************************************************************************/		

		/*
		   if((PDGnu == 12 && CC == 1)) // || PDGnu== -12) && CC == 1)// && Nph == 0 && Npi0dph == 0 && Npi0 ==0)// && Npim == 0 && Npip == 0)
		   {


		   double El_smear = smear_energy(El_true, EMsmear, rangen);
		   double Enu_reco;

		   if(El_smear > EM_thresh)
		   {
		//double prob = workingModel.oscProbSin(Enu,0.001*(osclen+(vertex_pos[2]/1000.0)));
		//double probsq = workingModel.oscProbSinSq(Enu,0.001*(osclen+(vertex_pos[2]/1000.0)));
		Enu_reco = El_smear + vertex_energy;
		hist.at( map_hist["nu_SBND_elike_intrinsic"]).Fill(Enu_reco*weight);
		} //end 200Mev smeared cut	


		}//end nu_e cc cut
		*/

		/************************************************************************************************
		 *				NC photon bit, based off yeon-jae's photon-flow 
		 * **********************************************************************************************/		
		int num_revertex = 100;
		double vertex_pos[3];


		//first go though individual photons

		//photon flow
		std::vector<myphoton> gammas;

		if (Nph!=0){
			for(int j=0; j<Nph;j++ ){

				myphoton temp_vec;
				temp_vec.lorentz.SetPxPyPzE(pph[j][0],pph[j][1],pph[j][2],Eph[j]);

				double temp_energy_smeared = detectors[0]->smear_energy(Eph[j], 0.15, rangen);

				if (temp_energy_smeared > 0.05) {

					temp_vec.lorentz.SetE(temp_energy_smeared);
					gammas.push_back(temp_vec);
					gammas.back().isPion=0;
				}

			}
		}


		//pion flow
		//pi0 gammas
		if (Npi0dph!=0){

			for(int j=0;j<Npi0dph;j++){
				myphoton temp_vec;
				temp_vec.lorentz.SetPxPyPzE(ppi0dph[j][0],ppi0dph[j][1],ppi0dph[j][2],Epi0dph[j]);

				double temp_energy_smeared = detectors[0]->smear_energy(Epi0dph[j], 0.15, rangen);


				if (temp_energy_smeared > 0.05) {

					temp_vec.lorentz.SetE(temp_energy_smeared);
					gammas.push_back(temp_vec);
					gammas.back().isPion=1;
				}

			}

		}

		detectors[0]->random_pos(rangen,vertex_pos); //Assign a random position for vertex in detector 

		for(auto & pho: gammas){
		}




		if( gammas.size() == 2){

			double ii=(gammas.at(0).lorentz+gammas.at(1).lorentz).M();
			if(gammas.at(0).isPion && gammas.at(1).isPion){
				//std::cout<<ii<<" "<<gammas.at(0).lorentz.E()<<" "<<gammas.at(0).isPion<<" "<<gammas.at(1).lorentz.E()<<" "<<gammas.at(1).isPion<<std::endl;
				hist.at(map_hist["nu_uBooNE_ncinvar_ncpi0"]).Fill((gammas.at(0).lorentz+gammas.at(1).lorentz).M());
			}else{
				std::cout<<ii<<" "<<gammas.at(0).lorentz.E()<<" "<<gammas.at(0).isPion<<" "<<gammas.at(1).lorentz.E()<<" "<<gammas.at(1).isPion<<std::endl;
				hist.at(map_hist["nu_uBooNE_ncinvar_nc2gamma"]).Fill(ii);
			}

		}

		for(int iv=0; iv<num_revertex; iv++){
		}//END REVERTEX


		//will be fullosc sample	
	}else if(file ==1) {




	}


	return 0;
}

int SBgeNicarus::tidyHistograms(){


	//hist.at(map_hist["nu_SBND_elike_misphoton"]).SetBinContent(1,10);



	return 0;
}
















/*************************************************************
 *************************************************************
 *		BEGIN example.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{


	std::string xml = "default.xml";
	int iarg = 0;
	opterr=1;
	int index; 
	int test_mode=0;
	std::string filename = "default.root";
	/*************************************************************
	 *************************************************************
	 *		Command Line Argument Reading
	 ************************************************************
	 ************************************************************/

	const struct option longopts[] = 
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"test",		required_argument,	0, 't'},
		{"file",		required_argument,	0, 'f'},
		{0,			no_argument, 		0,  0},
	};


	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:t:f:", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 'f':
				filename = optarg;//`strtof(optarg,NULL);
				break;

			case 't':
				test_mode = strtof(optarg,NULL);
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}

	}

	std::cout<<"GENTEST: starting"<<std::endl;
	SBgeNicarus testGen(xml);
	std::cout<<"GENTEST: doMC"<<std::endl;
	testGen.doMC();

}



double SBgeNicarus::vertexEnergy(int file){
	double Ehad = 0;

	int Np = *vmapI[file]["Np"];
	double *Ep =  vmapDA[file]["Ep"]->data;

	int Npim = *vmapI[file]["Npim"];
	double *Epim =  vmapDA.at(file).at("Epim")->data;

	int Npip = *vmapI[file]["Npip"];
	double *Epip =  vmapDA[file]["Epip"]->data;

	int No = *vmapI[file].at("No");
	double *Eo =  vmapDA[file].at("Eo")->data;
	int *pdgo =  vmapIA[file].at("pdgo")->data; 


	if(Np!=0){
		double p_kin_true = 0;
		double p_kin_smeared = 0;

		for(int j=0; j<Np; j++)
		{
			p_kin_true = Ep[j]-MPROTON;
			p_kin_smeared = detectors[0]->smear_energy(p_kin_true,psmear,rangen);
			if(p_kin_smeared>p_thresh)
			{
				Ehad += p_kin_smeared;
			}
		}
	} //end proton addition 

	if(Npip!=0){
		double pip_kin_true = 0;
		double pip_kin_smeared = 0;
		for(int j=0; j<Npip;j++)
		{
			pip_kin_true = Epip[j]-MPION;
			pip_kin_smeared = detectors[0]->smear_energy(pip_kin_true,pismear,rangen);
			if(pip_kin_smeared>pip_thresh)
			{
				Ehad += pip_kin_smeared+MPION;
			}


		} 
	}//end pi+addition

	if(Npim!=0){
		double pim_kin_true = 0;
		double pim_kin_smeared = 0;
		for(int j=0; j<Npim;j++)
		{
			pim_kin_true = Epim[j]-MPION;
			pim_kin_smeared = detectors[0]->smear_energy(pim_kin_true,pismear,rangen);
			if(pim_kin_smeared > pim_thresh)
			{
				Ehad += pim_kin_smeared+MPION;
			}

		} 
	}//end piminus addition

	if(No!=0){

		for(int j=0; j<No;j++)
		{

			if(pdgo[j]==321 || pdgo[j]==-321 || pdgo[j]==311){
				//std::cout<<pdgo[j]<<std::endl;	
				Ehad += detectors[0]->smear_energy(Eo[j]-MKAON,pismear,rangen)+MKAON;
			}

			if( pdgo[j]==3222 || pdgo[j]==3112 || pdgo[j]==3122){	

				Ehad += detectors[0]->smear_energy(Eo[j]-MSIGMA,pismear,rangen)+MSIGMA;
			}

		} 
	}//end Other 



	return Ehad;



}
