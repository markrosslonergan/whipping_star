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

class SBgeNicarus : public SBgeN{
	public:
	SBgeNicarus(std::string xmlname):SBgeN(xmlname){};
	bool eventSelection(int file);
	int fillHistograms(int file, int uni, double wei);
	int tidyHistograms();
	~SBgeNicarus(){};
};


bool SBgeNicarus::eventSelection(int file){

return true;
}

int SBgeNicarus::fillHistograms(int file, int uni, double wei){

	//std::cout<< *vmapD[file]["Enu"]<<" "<<vars_d.at(file).at(0)<<std::endl;	
	double en = *vmapD[file]["Enu"];//   vars_d.at(file).at(0);
	double p1 =  vmapDA[file]["Ep"]->data[0];

	hist.at( map_hist["nu_SBND_elike_dirt"] ).Fill(en);

return 0;
}

int SBgeNicarus::tidyHistograms(){
	hist.at(map_hist["nu_SBND_elike_misphoton"]).SetBinContent(1,10);
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


SBgeNicarus testGen(xml);
testGen.doMC();

}
