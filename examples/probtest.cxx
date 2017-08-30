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
#include "prob.h"
#include "SBNprob.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;









/*************************************************************
 *************************************************************
 *		BEGIN example.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{

	double T34=0;
	double T24=0;

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
				T34 = strtod(optarg,NULL);
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}

	}

	double L_DUNE = 1300;
	double L_SBND = 0.1;
	double L_ICARUS = 0.6;

	std::vector<std::string> modes = {"3p0vac","3p0mat","3p1vac","3p1mat","3p1mat_nonc"};
	std::ofstream logstream;
	logstream.open ("LOG.dat");


	double t12 = 33.5;
	double t23 = 41.6;
	double t13 = 8.5;
	double t14 = 0;
	double t24 = 20;
	double t34 = 0;

	double d13=0;
	double d24=0;
	double d34=0;

	double m21 = 7.5e-5;
	double m31 = 2.52e-3;
	double m41 = 1.0;

	logstream<<"T12 "<<t12<<std::endl;
	logstream<<"T23 "<<t23<<std::endl;
	logstream<<"T13 "<<t13<<std::endl;
	logstream<<"T14 "<<t14<<std::endl;
	logstream<<"T24 "<<t24<<std::endl;
	logstream<<"T34 "<<t34<<std::endl;

	logstream<<"D13 "<<d24<<std::endl;
	logstream<<"D24 "<<d24<<std::endl;
	logstream<<"D34 "<<d24<<std::endl;
	
	
	logstream<<"m21 "<<m21<<std::endl;
	logstream<<"m31 "<<m31<<std::endl;
	logstream<<"m41 "<<m41<<std::endl;
	
	logstream.close();


	for(auto &mode: modes){

		std::cout<<"Beginning mode :"<<mode<<std::endl;
		std::ofstream dunestream,icarusstream,sbndstream;

		dunestream.open ("DUNE_"+mode+".dat");
		icarusstream.open ("ICARUS_"+mode+".dat");
		sbndstream.open ("SBND_"+mode+".dat");


		std::vector<double> angles = {t12, t23, t13, t14, t24, t34};
		std::vector<double> phases = {d13,d24,d34};
		std::vector<double> mass = {m21,m31,m41};

		if(mode == "3p0vac" || mode == "3p0mat"){
			mass[2]=0;
			angles[3]=0; angles[4]=0; angles[5]=0;
		}


		SBNprob myprob(4, angles, phases, mass);

		if(mode == "3p0vac" || mode == "3p1vac")
		{
			myprob.setMatterEffect(false);
		}else{
			myprob.setMatterEffect(true);
		}
		if(mode == "3p1mat_nonc"){
			myprob.setNCMatterEffect(false);
		}



		//plotProbabilityMatter(int a, int b, double Emin, double Emax, double L, double percen, double n);
		myprob.plotProbabilityMatter(2, 1, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu->e
		myprob.plotProbabilityMatter(2, 2, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu->mu
		myprob.plotProbabilityMatter(2, 3, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu -> tau
		myprob.plotProbabilityMatter(2, 4, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu -> s
		myprob.plotProbabilityMatter(1, 4, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  e ->s

		myprob.plotProbabilityMatter(2, 1, -4, 1, L_SBND, 0.15, 1500,&sbndstream); //  mu->e
		myprob.plotProbabilityMatter(2, 2, -4, 1, L_SBND, 0.15, 1500,&sbndstream); //  mu->mu
		myprob.plotProbabilityMatter(2, 3, -4, 1, L_SBND, 0.15, 1500,&sbndstream); //  mu -> tau
		myprob.plotProbabilityMatter(2, 4, -4, 1, L_SBND, 0.15, 1500,&sbndstream); //  mu -> s
		myprob.plotProbabilityMatter(1, 4, -4, 1, L_SBND, 0.15, 1500,&sbndstream); //  e ->s

		myprob.plotProbabilityMatter(2, 1, -4, 1, L_ICARUS, 0.15, 1500,&icarusstream); //  mu->e
		myprob.plotProbabilityMatter(2, 2, -4, 1, L_ICARUS, 0.15, 1500,&icarusstream); //  mu->mu
		myprob.plotProbabilityMatter(2, 3, -4, 1, L_ICARUS, 0.15, 1500,&icarusstream); //  mu -> tau
		myprob.plotProbabilityMatter(2, 4, -4, 1, L_ICARUS, 0.15, 1500,&icarusstream); //  mu -> s
		myprob.plotProbabilityMatter(1, 4, -4, 1, L_ICARUS, 0.15, 1500,&icarusstream); //  e ->s

		myprob.setAntiNeutrinoMode(true);
		//AntiNeutrino
		myprob.plotProbabilityMatter(-2, 1, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu->e
		myprob.plotProbabilityMatter(-2, 2, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu->mu
		myprob.plotProbabilityMatter(-2, 3, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu -> tau
		myprob.plotProbabilityMatter(-2, 4, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu -> s
		myprob.plotProbabilityMatter(-1, 4, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  e ->s


		dunestream.close();
		sbndstream.close();
		icarusstream.close();

	}//end mode loop

	/*
	   if(false){//this is giff mode!
	   myprob.t34 = T34*3.14159/180.0;
	   myprob.init();



	   std::cout<<"# t34: "<<T34<<std::endl;
//int plotProbabilityMatter(int a, int b, double Emin, double Emax, double L, double percen, double n);
//		myprob.plotProbabilityMatter(1, 0, -1, 2, 1300, 0.1, 1500);
//		myprob.plotProbabilityMatter(1, 1, -1, 2, 1300, 0.1, 1500);
//		myprob.plotProbabilityMatter(1, 2, -1, 2, 1300, 0.1, 1500);
//		myprob.plotProbabilityMatter(1, 3, -1, 2, 1300, 0.1, 1500);
}
*/
}


