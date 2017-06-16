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


	std::vector<double> angles = {33, 40, 8, 20, 20, 20};
	//std::vector<double> angles = {33,40, 8, 0,0,0};

	std::vector<double> phases = {0,0,0};
	std::vector<double> mass = {7.5e-5, 2.552e-3, 1.0};

	SBNprob myprob(4, angles, phases, mass);
	myprob.setMatterEffect(true);

	myprob.plotProbabilityMatter(1, 0, -3, 2, 1300, 0.1, 1500);
	myprob.plotProbabilityMatter(1, 1, -3, 2, 1300, 0.1, 1500);
	myprob.plotProbabilityMatter(1, 2, -3, 2, 1300, 0.1, 1500);
	myprob.plotProbabilityMatter(1, 3, -3, 2, 1300, 0.1, 1500);



	if(false){//this is giff mode!
		myprob.t34 = T34*3.14159/180.0;
		myprob.init();



		std::cout<<"# t34: "<<T34<<std::endl;
		//int plotProbabilityMatter(int a, int b, double Emin, double Emax, double L, double percen, double n);
		myprob.plotProbabilityMatter(1, 0, -1, 2, 1300, 0.1, 1500);
		myprob.plotProbabilityMatter(1, 1, -1, 2, 1300, 0.1, 1500);
		myprob.plotProbabilityMatter(1, 2, -1, 2, 1300, 0.1, 1500);
		myprob.plotProbabilityMatter(1, 3, -1, 2, 1300, 0.1, 1500);
	}
}


