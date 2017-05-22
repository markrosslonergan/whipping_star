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
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovar.h"
#include "MCEventWeight.h"

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



std::cout<<"Strating covariance"<<std::endl;
SBNcovar my_covar(filename,xml);

my_covar.formCovarianceMatrix();

return 0;

	SBNspec bkg_spec("../data/precomp/SBN_bkg_all", xml);
	bkg_spec.compressVector();
	
	SBNchi test_chi(bkg_spec,my_covar.full_covariance);
	SBNspec sig_spec("../data/precomp/SBN_bkg_all", xml);


	int n = 75;
	double x[n], chi[n];
	for(int i=0; i<n; i++){
		double scale =0.25+i*0.025;		
		
		SBNspec loop_spec = sig_spec;
		loop_spec.Scale("elike_misphoton", scale);
		loop_spec.compressVector();

		x[i]=scale;
		chi[i]=test_chi.CalcChi(loop_spec);
		
		std::cout<<"scaling: "<<scale<<" "<<"Chi^2: "<<chi[i]<<std::endl;
	}
   	
	TGraph *gr  = new TGraph(n,x,chi);
	TFile * ff = new TFile("example_1.root","RECREATE");
	TCanvas *c1 = new TCanvas("c1","example_1");

	gr->SetTitle("Scaling of misidentified photon rate ");
	gr->GetXaxis()->SetTitle("Scale factor");
	gr->GetYaxis()->SetTitle("#chi^{2}");
	gr->Draw("ACP");

	c1->Write();
	ff->Close();





}
