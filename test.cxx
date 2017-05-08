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

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;


/*************************************************************
 *************************************************************
 *		BEGIN Main::sbnfit.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{


std::string xml = "default.xml";
int iarg = 0;
opterr=1;
int index; 
int test_mode=0;
/*************************************************************
 *************************************************************
 *		Command Line Argument Reading
 ************************************************************
 ************************************************************/

const struct option longopts[] = 
{
	{"xml", 		required_argument, 	0, 'x'},
	{"test",		required_argument,	0, 't'},
	{0,			no_argument, 		0,  0},
};


while(iarg != -1)
{
	iarg = getopt_long(argc,argv, "x:t:", longopts, &index);

	switch(iarg)
	{
		case 'x':
			xml = optarg;
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
//Test 1: just a signal, varying the normalisation, fit a landau, fit something more cpmplcated 3+1 sensitivity
//
//Test 2: landau fit, then inject a signal and other fit



double uboonepot = 0.5;


SBNspec bkgSpec("precomp/SBN_bkg_all", xml);
bkgSpec.Scale("uBooNE", uboonepot);
bkgSpec.compressVector();


SBNchi chi(bkgSpec);

/*SBNspec sig("precomp/SBN_bkg_all", "sbn.xml");

for(double x =0.5 ;x <= 1.5; x+=0.025){
	SBNspec loopSpec = sig;
	
	loopSpec.Scale("elike_misphoton", x);
	loopSpec.calcFullVector();
	loopSpec.compressVector();
	std::cout<<"scaling: "<<x<<" "<<"Chi^2: "<< chi.CalcChi(loopSpec)<<std::endl;
}

std::cout<<"Begining SBNosc study"<<std::endl;
*/



SBNspec bkgSpec2("precomp/SBN_bkg_all", xml);
bkgSpec2.Scale("uBooNE_elike_mismuon", 5);
bkgSpec2.Scale("uBooNE_mlike_intrinsic", 0.4);

bkgSpec2.compressVector();
SBNchi chi2(bkgSpec2);
SBNspec sigSpec2("precomp/SBN_bkg_all",xml);

/*
TF1 *fLan = new TF1("fLan","TMath::Landau(x,[0],[1],0)",0,5);


	for(int i=0; i<100; i++){
		double mpv = i*0.05;
		fLan->SetParameters(mpv, 1.3);
		
		SBNspec tempSpec = sigSpec2;
		tempSpec.Scale("uBooNE_elike_intrinsic", fLan);
		tempSpec.compressVector();

		std::cout<<mpv<<" "<<chi2.CalcChi(tempSpec)<<std::endl;

	}
*/

/*
std::cout<<"construct my fitter"<<std::endl;
SBNfit myfit(bkgSpec2, sigSpec2, 2);

std::vector<std::pair<std::string, int>> myin;
myin.push_back(std::make_pair("uBooNE_elike_mismuon",0) );
myin.push_back(std::make_pair("uBooNE_mlike_intrinsic",1) );

std::vector<double> init  {0.99,1.001};
std::vector<double> low  {0.1,0.1};
std::vector<double> up  {10,10};
std::vector<double> step  {0.02,0.02};
std::vector<std::string> nam  {"p1","p2"};


std::cout<<"initialize my fitter"<<std::endl;
myfit.initialize_norm(myin);

std::cout<<"set stuff"<<std::endl;
myfit.setInitialValues(init);
myfit.setNames(nam);
myfit.setLowerValues(low);
myfit.setUpperValues(up);
myfit.setStepSizes(step);

std::cout<<"minimize!"<<std::endl;
myfit.Minimize();

std::cout<<"Minimized! chi^2: "<<myfit.BFchi<<" p1: "<<myfit.BFparam[0]<<" p2: "<<myfit.BFparam[1]<<" #:"<<myfit.BFncalls<<std::endl;

*/

	SBNosc injectSig("precomp/SBN_bkg_all",xml);
	neutrinoModel signalModel(1,0.15,0.18);
	signalModel.numsterile = 1;
		
	injectSig.load_model(signalModel);
	injectSig.OscillateThis();
	//injectSig.poissonScale();
	injectSig.compressVector();

	//injectSig.writeOut("inject.root");

SBNosc testSig("precomp/SBN_bkg_all",xml);



SBNfit3p1  fit3p1(injectSig, testSig, 3);
fit3p1.setMethod("GSLMultiMin", "BFGS2");


std::vector<std::string> nam 	{"m4", "Ue4", "Um4"};
std::vector<double> init  	{1.0, 0.05,  0.05}; 
std::vector<double> low  	{0.1,  0,     0};
std::vector<double> up    	{10,   1,     1}; 
std::vector<double> step  	{0.2,  0.05, 0.05}; 
std::vector<int> fix 		{1,    0,      0};

fit3p1.setInitialValues(init);
fit3p1.setNames(nam);
fit3p1.setLowerValues(low);
fit3p1.setUpperValues(up);
fit3p1.setStepSizes(step);
fit3p1.setFixed(fix);

fit3p1.Minimize();

std::cout<<"Minimized! chi^2: "<<fit3p1.bf_chi<<" Ue1: "<<fit3p1.bf_params[1]<<" Um1: "<<fit3p1.bf_params[2]<<" #:"<<fit3p1.num_func_calls<<std::endl;



return 0;
				SBNosc oscSig("precomp/SBN_bkg_all","sbn.xml");
				oscSig.setAppMode();
				oscSig.Scale("uBooNE",uboonepot);		


	TCanvas *c1 = new TCanvas("c1","c1",600,400);
	TH2F *hcontz = new TH2F("hcontz","",100,-5,0, 100,-2,2);
	hcontz->GetXaxis()->SetTitle("#sin^{2} 2 #theta_{e #mu}");
	hcontz->GetYaxis()->SetTitle("#Delta m^{2}_{41} (eV^{2})");

	for(double m = -2.00; m <=2.04; m=m+0.04){
	 for(double sins2 = 0.0 ; sins2 >= -5; sins2 = sins2 - 0.05){

	 	double uei = 0.1;
		double umi = sqrt(pow(10,sins2))/(2*uei);
		
				double imn[3] = {sqrt(pow(10,m)),0,0};
				double iue[3] = {umi,0,0};
				double ium[3] = {uei,0,0};
				double iph[3] = {0,0,0};
				neutrinoModel signalModel(imn,iue,ium,iph);
				signalModel.numsterile = 1;
				
				oscSig.load_model(signalModel);
				std::vector<double> ans = oscSig.Oscillate();

				std::cout<<m<<" "<<sins2<<" "<<chi.CalcChi(ans)<<std::endl;
			
		hcontz->SetBinContent( 1+floor(-(-5.0-sins2)/0.05+0.00001) , 1+floor(-(-2.00-m)/0.04+0.00001), chi.CalcChi(ans));


	 }
	}
		  Double_t contours[1];
        	  contours[0] = 1.64;
		  hcontz->SetContour(1, contours);

				c1->cd();

				hcontz->Draw("CONT3");
				TFile * ff = new TFile("test.root","RECREATE");
				ff->cd();
				c1->Write();
				ff->Close();



}
