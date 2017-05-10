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

/*************************************************************
 *************************************************************
 *		Example 1:
 ************************************************************
 ************************************************************/
if(test_mode==1){


	SBNspec bkg_spec("../data/precomp/SBN_bkg_all", xml);
	bkg_spec.compressVector();
	SBNchi test_chi(bkg_spec);
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




return 0;

/*************************************************************
 *************************************************************
 *		Example 2:
 ************************************************************
 ************************************************************/

} else if (test_mode ==2){



SBNspec bkg_spec("../data/precomp/SBN_bkg_all", xml);
bkg_spec.Scale("uBooNE_elike_mismuon", 2);
bkg_spec.Scale("uBooNE_mlike_intrinsic", 0.5);
bkg_spec.compressVector();

SBNchi test_chi(bkg_spec);
SBNspec sig_spec("../data/precomp/SBN_bkg_all",xml);


	TF1 *fLan = new TF1("fLan","TMath::Landau(x,[0],[1],0)",0,5);

	int n = 100;
	double x[n], chi2[n];

	for(int i=0; i<n; i++){
		double mpv = i*0.05;
		fLan->SetParameters(mpv, 1.3);
		
		SBNspec loop_spec = sig_spec;
		loop_spec.Scale("uBooNE_elike_intrinsic", fLan);
		loop_spec.compressVector();

		x[i]=mpv;
		chi2[i]=test_chi.CalcChi(loop_spec);
	
		std::cout<<"MPV: "<<mpv<<" chi^2: "<<chi2[i]<<std::endl;

	}

	TGraph *gr2  = new TGraph(n,x,chi2);
	TFile  *ff2 = new TFile("example_2.root","RECREATE");
	TCanvas *c2 = new TCanvas("c2","example_2");

	gr2->SetTitle("Scaling intrinsic nu_e by energy dependant landau");
	gr2->GetXaxis()->SetTitle("MPV");
	gr2->GetYaxis()->SetTitle("#chi^{2}");
	gr2->Draw("ACP");

	c2->Write();
	ff2->Close();


return 0;
/*************************************************************
 *************************************************************
 *		Example 3: Sterile neutrino, 3+1 sensitivity
 ************************************************************
 ************************************************************/
}else if(test_mode==3){

	double uboonepot=2;
	SBNspec bkg_spec("../data/precomp/SBN_bkg_all", xml);
	bkg_spec.Scale("uBooNE", uboonepot);
	bkg_spec.compressVector();
	
	SBNchi test_chi(bkg_spec);
	
	SBNosc oscSig("../data/precomp/SBN_bkg_all",xml);
	oscSig.setAppMode();
	oscSig.Scale("uBooNE",uboonepot);		


	TCanvas *c1 = new TCanvas("c1","c1",600,400);
	TH2F *hcontz = new TH2F("hcontz","MicroBooNE 3+1 90\% C.L ",100,-5,0, 100,-2,2);
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

			double tchi=test_chi.CalcChi(ans); 
			std::cout<<"Dm^2: "<<m<<" sin^2 th: "<<sins2<<" chi^2: "<<tchi<<std::endl;

			hcontz->SetBinContent( 1+floor(-(-5.0-sins2)/0.05+0.00001) , 1+floor(-(-2.00-m)/0.04+0.00001), tchi);


	 }
	}
			  Double_t contours[1];
        		  contours[0] = 1.64;
			  hcontz->SetContour(1, contours);

				c1->cd();

				hcontz->Draw("CONT3");
				TFile * ff = new TFile("example_3.root","RECREATE");
				ff->cd();
				c1->Write();
				ff->Close();

return 0;
/*************************************************************
 *************************************************************
 *		Example 3: Sterile neutrino, 3+1 sensitivity
 ************************************************************
 ************************************************************/
} else if(test_mode ==4){


	SBNspec bkg_spec("../data/precomp/SBN_bkg_all", xml);
	bkg_spec.compressVector();
	
	SBNchi test_chi(bkg_spec);
	
	SBNspec sig_spec("../data/precomp/SBN_bkg_all",xml);
	sig_spec.compressVector();



	SBNfit myfit(bkg_spec, sig_spec, 2);

	std::vector<std::pair<std::string, int>> myin;
	myin.push_back(std::make_pair("uBooNE_elike_mismuon",0) );
	myin.push_back(std::make_pair("uBooNE_mlike_intrinsic",1) );

	std::vector<double> init  {0.99,1.001};
	std::vector<double> low  {0.1,0.1};
	std::vector<double> up  {10,10};
	std::vector<double> step  {0.02,0.02};
	std::vector<std::string> nam  {"p1","p2"};


	myfit.initialize_norm(myin);

	myfit.setInitialValues(init);
	myfit.setNames(nam);
	myfit.setLowerValues(low);
	myfit.setUpperValues(up);
	myfit.setStepSizes(step);

	std::cout<<"minimize!"<<std::endl;
	myfit.Minimize();

	std::cout<<"Minimized! chi^2: "<<myfit.bf_chi<<" p1: "<<myfit.bf_params[0]<<" p2: "<<myfit.bf_params[1]<<" #:"<<myfit.num_func_calls<<std::endl;

return 0;
/*************************************************************
 *************************************************************
 *		Example 5: Injected Signal 3+1
 ************************************************************
 ************************************************************/

} else if(test_mode==5){

	SBNosc inject_sig("../data/precomp/SBN_bkg_all",xml);
	
	double dmSq=1;
	double UE4=0.15;
	double UM4=0.18;
	neutrinoModel signalModel(dmSq,UE4,UM4);
	signalModel.numsterile = 1;
		
	inject_sig.load_model(signalModel);
	inject_sig.OscillateThis();
	inject_sig.poissonScale();
	inject_sig.compressVector();

	inject_sig.writeOut("example_5_injected_signal.root");


	SBNosc test_sig("../data/precomp/SBN_bkg_all",xml);



	SBNfit3p1  fit3p1(inject_sig, test_sig, 3);
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
}			

}
