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


complex_matrix M(3);

M.real(0,0) = 1.0;
M.real(1,1) = 2.0;
M.real(2,2) = 3.0;

M.real(0,1) = 0.1;
M.real(1,0) = 0.2;
M.real(0,2) = 0.3;
M.real(1,2) = 0.4;
M.real(2,0) = 0.5;
M.real(2,1) = 0.6;

M.imag(0,0) = -1.0;
M.imag(1,1) = -2.0;
M.imag(2,2) = -3.0;

M.imag(0,1) = -0.1;
M.imag(1,0) = -0.2;
M.imag(0,2) = -0.3;
M.imag(1,2) = -0.4;
M.imag(2,0) = -0.5;
M.imag(2,1) = -0.6;

//M.matrixExp();

//M.print();

M.mult(&M);
//M.print();




SBNprob myprob(4);

//double ans = myprob.probabilityMatterExact(2,1,5.0,1300);
//std::cout<<ans<<std::endl;

//return 0;

for(double ee=-2; ee<2; ee+=0.0033){
	double Ee = pow(10,ee);
	double ansA = myprob.probabilityMatterExact(1,0,Ee,1300);
	double ansD = myprob.probabilityMatterExact(1,1,Ee,1300);
	double ansT = myprob.probabilityMatterExact(1,2,Ee,1300);

std::cout<<Ee<<" "<<ansA<<" "<<ansD<<" "<<ansT<<std::endl;
}


}


