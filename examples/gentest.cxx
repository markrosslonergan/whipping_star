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
	TLorentzVector lorentz_t;
	int isPion;
	double convlength;
	double *convpoint;

};



class SBgeNC : public SBgeN{
	public:
		SBgeNC(std::string xmlname):SBgeN(xmlname){
		
		outstream.open ("NC_pi0_uBooNE.dat");
		outstream<<"# (1 ) True incoming neutrino energy"<<std::endl;
		outstream<<"# (2 ) True outgoing neutrino energy"<<std::endl;
		outstream<<"# (3 ) True Energy Photon 1 "<<std::endl;
		outstream<<"# (4 ) Reconstructed Energy Photon 1 "<<std::endl;
		outstream<<"# (5 ) True Energy Photon 2 "<<std::endl;
		outstream<<"# (6 ) Reconstructed Energy Photon 3 "<<std::endl;	
		outstream<<"# (7 ) True angle between photons "<<std::endl;	
		outstream<<"# (8 ) Reconstructed angle between photons "<<std::endl;
		outstream<<"# (9 ) True Momentum in X, Photon 1 "<<std::endl;
		outstream<<"# (10) True Momentum in Y, Photon 1 "<<std::endl;
		outstream<<"# (11) True Momentum in Z, Photon 1 "<<std::endl;
		outstream<<"# (12) True Momentum in X, Photon 2 "<<std::endl;
		outstream<<"# (13) True Momentum in Y, Photon 2 "<<std::endl;
		outstream<<"# (14) True Momentum in Z, Photon 2 "<<std::endl;
		outstream<<"# (15) Reconstructed sqrt(Invariant Mass) of photon pair"<<std::endl;
		outstream<<"# (16) Reconstructed vertex energy, 0 if no reconstructed vertex"<<std::endl;
		outstream<<"# (17) Event weight, to be applied to all events"<<std::endl;
		outstream<<"# (18) PDG code of parent neutrino"<<std::endl;


		};
		bool eventSelection(int file);
		int fillHistograms(int file, int uni, double wei);
		int tidyHistograms();
		~SBgeNC(){};


		std::ofstream outstream;



		int total_NC;
		int total_CC;

		double vertexEnergy(int file);


};

bool SBgeNC::eventSelection(int file){

	return true;
}

int SBgeNC::fillHistograms(int file, int uni, double wei){
	double vertex_pos[3] = {0,0,0};
	double convlength_thresh = 5.0;
	double photon_thresh = 0.05;
	double overall_reco_eff = 0.8;
	double angle_resolution = 3.0*3.14159/180.0;
	double muon_eff = 0.8;
	double additional_eff_no_vis_vertex = 0.1/0.8;
	double pot_mod = detectors[file]->potmodifier*66.0*detectors[file]->proposal_modifier;


	for(double ee = 0.14; ee< 40; ee+=0.01){
//		std::cout<<"TEST: "<<ee<<" "<<detectors[file]->pion_track_length(ee)<<std::endl;

	}

	

//	exit(EXIT_FAILURE);

	std::vector<myphoton> gammas;
	std::vector<myphoton> muons;
	std::vector<myphoton> electrons;

		double Enu_true = *vmapD[file]["Enu"];//   vars_d.at(file).at(0);
		double El_true =  *vmapD[file]["El"];
		double weight = *vmapD[file]["weight"];

		int PDGnu = *vmapI[file]["PDGnu"];
		int CC = *vmapI[file]["CC"];
		int NC = *vmapI[file]["NC"];
		int Npi0dph = *vmapI[file]["Npi0dph"];
		int Nph = *vmapI[file]["Nph"];

		double *Eph =  vmapDA[file]["Eph"]->data; 
		double *Epi0dph =  vmapDA[file]["Epi0dph"]->data; 
		double *pl =  vmapDA[file]["pl"]->data;

		double (*pph)[3] = vmapDA[file]["pph"]->data3; 
		double (*ppi0dph)[3] =  vmapDA[file]["ppi0dph"]->data3; 

		/********************************************************************
		 *			Starting main cutflow
		 ******************************************************************/


		//	detectors[0]->random_pos(rangen,vertex_pos); //Assign a random position for vertex in detector 
		double osclen = detectors[file]->osc_length(rangen);


		//Is there a visible vertex and how much energy is there!
		double Enu_reco = 0;
		bool vis_vertex = false;
		double vertex_energy = 	this->vertexEnergy(file);


		if(CC) total_CC++;
		if(NC) total_NC++;

		/************************************************************************************************
		 *				CCbit, muon and electron! 
		 ***********************************************************************************************/	
		bool observable_lepton = false;

		if ((PDGnu == 14 && CC == 1) || (PDGnu== -14 && CC == 1)){
			double El_smear = detectors[file]->smear_energy(El_true, 0.15, rangen);

			myphoton temp_vec;
			detectors[file]->random_pos(rangen,vertex_pos); //Assign a random position for vertex in detector 

			double Lmu = detectors[file]->muon_track_length(El_true);

			temp_vec.lorentz.SetPxPyPzE(pl[0],pl[1],pl[2],El_smear);

			double endpos[3] = {0,0,0};
			double temp_3vec[3] = {pl[0],pl[1],pl[2]};
			double observable_L = 0;

			detectors[file]->get_endpoint(vertex_pos, Lmu, temp_3vec, endpos);
			//std::cout<<"Plep: "<<pl[0]<<" "<<pl[1]<<" "<<pl[2]<<std::endl;
			//std::cout<<"Vertex: "<<vertex_pos[0]<<" "<<vertex_pos[1]<<" "<<vertex_pos[2]<<std::endl;
			//std::cout<<"End: "<<endpos[0]<<" "<<endpos[1]<<" "<<endpos[2]<<std::endl;


			/*
			if( !detectors[file]->is_fully_contained(vertex_pos, endpos) ) {
				observable_lepton = true;
			}else if (Lmu > 10.0){
				observable_lepton = true;
			}
			observable_L = Lmu;
			*/


			
			if( detectors[file]->is_fully_contained(vertex_pos, endpos) ) {
				observable_L = Lmu;
				//std::cout<<"Containted"<<std::endl;

				if(observable_L > 50.0){
					observable_lepton = true;
				}

			}
			else {
				observable_L = detectors[file]->track_length_escape(vertex_pos,endpos);
				//std::cout<<"Exiting"<<std::endl;

				if(observable_L > 100.0){
					observable_lepton = true;
				}
			}
			

			temp_vec.convlength = observable_L;
			muons.push_back(temp_vec);

		}

		if(!observable_lepton){
			for(auto & mu: muons){
				vertex_energy += mu.lorentz.E();
			}
		}

		//Check if we actually have a "visibe vertex"
		if(vertex_energy > vertex_thresh)
		{
			vis_vertex = true;
		} else 
		{
			vis_vertex = false;
		}



		/************************************************************************************************
		 *				NC photon bit, based off yeon-jae's photon-flow 
		 * **********************************************************************************************/		
		int num_revertex = 1;

		//checks to see is it a NC event?
		if(!observable_lepton){

			// "2 photon showers reconstructed close to the pi0 invariant mass, associated with a reconstructed vertex inside the fiducial volume, and no electron shower or muon track associated with the vertex"
			//first go though individual photons

			//photon flow

			if (Nph!=0){
				for(int j=0; j<Nph;j++ ){

					myphoton temp_vec;
					temp_vec.lorentz.SetPxPyPzE(pph[j][0],pph[j][1],pph[j][2],Eph[j]);

					temp_vec.lorentz_t = temp_vec.lorentz;

					double temp_energy_smeared = detectors[file]->smear_energy(Eph[j], 0.15, rangen);

					temp_vec.lorentz.SetE(temp_energy_smeared);
					if (temp_energy_smeared > photon_thresh) {

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

					double temp_energy_smeared = detectors[file]->smear_energy(Epi0dph[j], 0.15, rangen);

					temp_vec.lorentz_t = temp_vec.lorentz;

					temp_vec.lorentz.SetE(temp_energy_smeared);

					if (temp_energy_smeared > photon_thresh) {
						gammas.push_back(temp_vec);
						gammas.back().isPion=1;
					}

				}

			}


			for(int iv=0; iv<num_revertex; iv++){

				double wei = weight*pot_mod/((double)num_revertex);

				detectors[file]->random_pos(rangen,vertex_pos); //Assign a random position for vertex in detector 

				//For every photon give it a conv pt and conv length
				for(auto & pho: gammas){
					double tmpcovpoint[3];
					pho.convlength = detectors[file]->photon_conversion_length(pho.lorentz.E(), rangen);
					double temp_3vec[3] = {pho.lorentz.Px(),pho.lorentz.Py(),pho.lorentz.Pz()};
					detectors[file]->get_endpoint(vertex_pos, pho.convlength, temp_3vec, tmpcovpoint);
					pho.convpoint = tmpcovpoint;
				}

				//is there 2 photons  in fiducial?
				if( gammas.size() == 2 &&  detectors[file]->is_fiducial(vertex_pos)  ){
					//Are both start points in fiducial?
					if( detectors[file]->is_fiducial(gammas.at(0).convpoint) && 	detectors[file]->is_fiducial(gammas.at(1).convpoint) ){
						//Do they both convert more than 5cm away and have energy > 60mev?
						if(gammas.at(0).convlength > convlength_thresh && gammas.at(1).convlength > convlength_thresh && gammas.at(0).lorentz.E() > photon_thresh && gammas.at(1).lorentz.E() > photon_thresh  ){

							wei = wei*overall_reco_eff;
							if(!vis_vertex){
								wei = wei*additional_eff_no_vis_vertex;	
							}

							//double invarM=(gammas.at(0).lorentz+gammas.at(1).lorentz).M();
							double theta_z_0 = acos( gammas.at(1).lorentz.Pz()/gammas.at(1).lorentz_t.Vect().Mag());
							double theta_z_1 = acos( gammas.at(1).lorentz.Pz()/gammas.at(1).lorentz_t.Vect().Mag());

							double theta = gammas.at(0).lorentz.Angle(gammas.at(1).lorentz.Vect() );
							double theta_smear = detectors[0]->smear_angle(theta, angle_resolution, rangen);
							//std::cout<<"TH "<<theta<<" "<<theta_smear<<" "<<theta*180.0/3.14159<<" "<<theta_smear*180.0/3.14159<<std::endl;
							double invarM= sqrt( 2.0*gammas.at(0).lorentz.E()*gammas.at(1).lorentz.E()*(1.0-cos(theta_smear) ));

							double recoE=vertex_energy + gammas.at(0).lorentz.E()+gammas.at(1).lorentz.E();



							std::string prenam = detectors[file]->name;
							prenam =  "nu_" + prenam + "_";

							if(muons.size() > 0){
								hist.at(map_hist[prenam+"nchadron_ccmuon"]).Fill(vertex_energy,wei*muon_eff);
								hist.at(map_hist[prenam+"ncinvar_ccmuon"]).Fill(invarM,wei*muon_eff);
							}else if(gammas.at(0).isPion && gammas.at(1).isPion){
								hist.at(map_hist[prenam+"nchadron_ncpi0"]).Fill(vertex_energy,wei);
								hist.at(map_hist[prenam+"ncinvar_ncpi0"]).Fill(invarM,wei);

								//************************ Stream to text file of pi0 ************************
								if(NC){
									outstream<<Enu_true<<" "<<El_true<<" "<<gammas.at(0).lorentz_t.E()<<" "<<gammas.at(0).lorentz.E()<<" "<<gammas.at(1).lorentz_t.E()<<" "<<gammas.at(1).lorentz.E()<<" "<<theta<<" "<<theta_smear<<" "<<gammas.at(0).lorentz.Px()<<" "<<gammas.at(0).lorentz.Py()<<" "<<gammas.at(0).lorentz.Pz()<<" "<<gammas.at(1).lorentz.Px()<<" "<<gammas.at(1).lorentz.Py()<<" "<<gammas.at(1).lorentz.Pz()<<" "<<invarM<<" "<<vertex_energy<<" "<<weight<<" "<<PDGnu<<std::endl;

								}

							}else{
								hist.at(map_hist[prenam+"nchadron_nc2gamma"]).Fill(vertex_energy,wei);
								hist.at(map_hist[prenam+"ncinvar_nc2gamma"]).Fill(invarM,wei);
							}

							if(vertex_energy < 0.5&& invarM < 0.3){
								if(muons.size() > 0){
									hist.at(map_hist[prenam+"ncreco_ccmuon"]).Fill(recoE,wei*muon_eff);
								}else if(gammas.at(0).isPion && gammas.at(1).isPion){
									hist.at(map_hist[prenam+"ncreco_ncpi0"]).Fill(recoE,wei);
								}else{
									hist.at(map_hist[prenam+"ncreco_nc2gamma"]).Fill(recoE,wei);
								}
							}
						}
					}
				}

			}//END REVERTEX

		}


		return 0;
}

int SBgeNC::tidyHistograms(){


	//hist.at(map_hist["nu_SBND_elike_misphoton"]).SetBinContent(1,10);
	std::cout<<"Total NC: "<<total_NC<<" Total CC: "<<total_CC<<" Total : "<<total_NC+total_CC<<std::endl;
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
	SBgeNC testGen(xml);
	std::cout<<"GENTEST: doMC"<<std::endl;
	testGen.doMC();

}



double SBgeNC::vertexEnergy(int file){
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
