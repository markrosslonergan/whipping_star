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

	//SBNspec bkg_spec("../../data/precomp/SBN_bkg_all", xml);
	//bkg_spec.Scale("uBooNE",0.5);
	SBNosc bkgSpec("../../data/precomp/SBN_bkg_all",xml);
	bkgSpec.Scale("uBooNE",0.5);


	SBNosc oscSig("../../data/precomp/SBN_bkg_all",xml);
	oscSig.Scale("uBooNE",0.5);
	oscSig.setBothMode();
//	oscSig.setWierdMode();
//	oscSig.setDisMode();
//	oscSig.setAppMode();
//	oscSig.setDisEMode();

	oscSig.writeOut("zenB.root");
	//0.46 0.15 0.13 0.77 0.13 0.14 5.56
	double imn[3] = {sqrt(0.46),sqrt(0.77),0};
	double iue[3] = {0.15,0.13,0};
	double ium[3] = {0.13,0.14,0};
	double iph[3] = {5.56,0,0};

	//construct a signalModel
	neutrinoModel signalModel(imn,iue,ium,iph);
	signalModel.numsterile = 2; //this isnt really necessary as it can tell from imn, but nice for reading

	oscSig.load_model(signalModel);

	oscSig.OscillateThis();
	oscSig.writeOut("zenS.root");
	oscSig.compressVector();

	std::vector<int> cols = {kBlack, kBlack, kBlue-7, kGreen-6, kRed-7,kOrange-3, kMagenta-3, kBlue-7, kRed-7};


	TLegend legE(0.6,0.35,0.875,0.875);
	TLegend legM(0.6,0.6,0.875,0.875);


	TCanvas *c = new TCanvas("","",2000,2000);
	c->Divide(2,2);
	c->cd(1);
	THStack * he 	   = new THStack("#nu_{e} like events ","#nu_{e} like events");
	int n=0;

	std::vector<TH1D> temp = bkgSpec.hist;
	std::cout<<"temp length "<<temp.size()<<std::endl;


	TH1D *sigE = (TH1D*)temp.at(0).Clone("sige");
	TH1D *sigM = (TH1D*)bkgSpec.hist.at(8).Clone("sigm");
	std::vector<double> errE;
	std::vector<double> errM;

	std::cout<<"Length comp: "<<oscSig.compVec.size()<<std::endl;
	for(int i=0; i<11; i++){
		std::cout<<"E: "<<i<<std::endl;
		sigE->SetBinContent(i+1,oscSig.compVec.at(i));
		errE.push_back(sqrt( oscSig.compVec.at(i)));
	}	
	for(int i=0; i<19; i++){
		std::cout<<"M: "<<i<<std::endl;
		sigM->SetBinContent(i+1,oscSig.compVec.at(i+11));
		errM.push_back(sqrt( oscSig.compVec.at(i+11)));
	}

	std::cout<<"Done filling"<<std::endl;	
	sigE->SetMarkerStyle(21);
	sigE->SetError(&errE[0]);
	sigM->SetError(&errM[0]);
	sigM->SetMarkerStyle(21);
	sigM->SetMarkerColor(kBlack);
	sigE->SetMarkerColor(kBlack);
	sigE->SetLineColor(kBlack);
	sigM->SetLineColor(kBlack);
	sigE->SetMarkerSize(2);
	sigM->SetMarkerSize(2);
	sigE->SetLineWidth(1);
	sigM->SetLineWidth(1);
	sigE->Scale(1,"width");
	sigM->Scale(1,"width");

	std::cout<<"temp"<<std::endl;
	for(auto &h: temp){
		h.Scale(1,"width,nosw2");
		h.SetLineColor(kBlack);
		h.SetFillColor(cols.at(n));
		n++;
	}

	std::cout<<"adding"<<std::endl;
	he->Add(&temp.at(2) );				
	he->Add(&temp.at(3) );				
	he->Add(&temp.at(4) );				
	he->Add(&temp.at(5) );				
	he->Add(&temp.at(6) );

	legE.AddEntry(&temp.at(2),"Intrinsic #nu_{e} CC", "f");
	legE.AddEntry(&temp.at(3),"Mis ID NC", "f");
	legE.AddEntry(&temp.at(4),"Mis ID #nu_{#mu} CC", "f");
	legE.AddEntry(&temp.at(5),"Dirt", "f");
	legE.AddEntry(&temp.at(6),"Cosmics", "f");
	


	he->Draw();	
	he->GetYaxis()->SetTitle("Events/GeV");
	he->GetXaxis()->SetTitle("E_{#nu} Reco [GeV]");
	he->GetYaxis()->SetTitleOffset(1.5);
	he->SetMaximum(500);
	std::cout<<"drawi sigE"<<std::endl;
	sigE->Draw("E1,same");

	std::string tl = "3+2 Best Fit";//#Delta m_{41}^{2}= 0.46 eV^{2}";//,  #Delta m_{51}^2 = 0.77 eV^{2}/\nU_{e4} = 0.15, U_{e5} = 0.13\nU_{#mu 4} = 0.13, U_{#mu 5} = 0.14";
	legE.AddEntry(sigE, tl.c_str(),"lep");
	legE.Draw();
	
	c->cd(2);
	THStack * hm 	   = new THStack("#nu_{#mu} like events","#nu_{#mu} like events");
	legM.AddEntry(&temp.at(7),"Intrinsic #nu_{#mu} CC", "f");
	legM.AddEntry(&temp.at(8),"Mis ID #pi^{#pm}", "f");


	std::cout<<"adding M7"<<std::endl;
	hm->Add(&temp.at(8) );				
	std::cout<<"adding M8"<<std::endl;
	hm->Add(&temp.at(7) );
	std::cout<<"madded"<<std::endl;
	std::cout<<"drawing stack M"<<std::endl;
	hm->Draw();
	
	hm->GetYaxis()->SetTitle("Events/GeV");
	hm->GetXaxis()->SetTitle("E_{#nu} Reco [GeV]");
	hm->GetYaxis()->SetTitleOffset(1.5);
	sigM->SetFillColor(kWhite);
	sigM->Draw("E1,same");
	legM.AddEntry(sigM, "3+2 Best Fit","lep");	
//	legM.AddEntry(sigM, "#Delta m_{41}^{2}= 0.46 eV^{2}, #Delta m_{51}^{2} = 0.77 eV^{2}","");	
	legM.Draw();

	c->cd(3);	
	TLegend leg(0.2,0.3,0.875,0.95);
	leg.AddEntry(sigM, "3+2 Best Fit Signal","lep");
	leg.AddEntry(sigM, "#Delta m_{41}^{2}= 0.46 eV^{2}, #Delta m_{51}^{2} = 0.77 eV^{2}  ","");
	leg.AddEntry(sigM, "U_{e4} = 0.15, U_{e5} = 0.13","");
	leg.AddEntry(sigM, "U_{#mu 4} = 0.13, U_{#mu 5} = 0.14","");
	leg.AddEntry(sigM, "#phi_{45} = 5.54 ", "");

	leg.Draw();
	c->SaveAs("zen.pdf");

	//Its the full osc...Thats why


}
