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
#include "TAxis.h"

#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

#include <fstream>

using namespace sbn;

int oscillateFast(double dm, double sinth, TH1D* input, TTree * tin){

	double Reco = 0;
	double True = 0;
	double L = 0;
	double weight = 0;

	tin->SetBranchAddress("Ereco",&Reco);
	tin->SetBranchAddress("Etruth",&True);
	tin->SetBranchAddress("L",&L);
	tin->SetBranchAddress("weight",&weight);

	int Nent = tin->GetEntries();

	for(int i=0; i<Nent; i++){
		tin->GetEntry(i);
		//std::cout<<" R: "<<Reco<<" T: "<<True<<" L: "<<L<<" weight: "<<weight<<std::endl;
		double osc = sinth*pow(sin(1.27*dm*L/True),2.0);
		input->Fill(Reco, osc*weight/((double)Nent) );
	}

	return 0;
}


int oscillate(double dm, double sinth, TH1D* input){
	std::ifstream in;
	in.open("/home/mark/work/uBooNE/request_of_MB_data/data_release_public/nu_only_mode/miniboone_numunuefullosc_ntuple.txt");

	int n=0;
	while (1) {
		//containing information on reconstructed neutrino energy, true neutrino energy, neutrino baseline, and event weight for each event
		n++;
		double Reco=0;
		double True=0; 
		double L = 0;
		double weight = 0;
		in >> Reco;
		in >> True;
		in >> L;
		in >> weight;	

		if (!in.good()) break;
		//std::cout<<n<<" R: "<<Reco/1000.0<<" T: "<<True/1000.0<<" L: "<<L/100000.0<<" weight: "<<weight<<std::endl;
		double osc = sinth*pow(sin( 1.27*dm*(L/100000.0)/(True/1000.0) ),2.0);
		//std::cout<<Reco/1000.0<<" "<<osc<<std::endl;
		if(osc>1){std::cout<<"ERROR prob<1"<<std::endl;exit(EXIT_FAILURE);}
		input->Fill(Reco/1000.0, osc*weight);
	}
	in.close();

	return 0;
}


/*

   class SBNfitMiniBooNE : public SBNfit {

   protected:
   public:
   TFile *f ;//=  new TFile("/home/mark/work/sbnfit_reduce/data/miniboone/fullosc_ntuple_nu.root");
   TTree * tup ;//= (TTree*)f->Get("fullosc_nu_mode;1");
   SBNspec data_spec;
   SBNspec sig_spec;

   SBNfitMiniBooNE(SBNspec inData, SBNspec inSig, int npa) : SBNfit(inData,inSig,npa), sig_spec(inSig), data_spec(inData) {
   f =  new TFile("/home/mark/work/sbnfit_reduce/data/miniboone/fullosc_ntuple_nu.root");
   tup = (TTree*)f->Get("fullosc_nu_mode;1");

   };

   double MinimizerCalcChi(const double * X){
   num_func_calls++;
   SBNspec tempspec = bk_spec; 

   double Dm = pow(10,X[0]);
   double sinsq2th = pow(10,X[1]);



   bkgSpec = ; 
   this->loadBkg();

   std::vector<double> ans = tempOsc.Oscillate();

   lastChi =this->CalcChi(ans);
   return lastChi;

   };

   };

*/
















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

	std::cout<<"Getting ntuple files"<<std::endl;
	TFile *f =  new TFile("/home/mark/work/sbnfit_reduce/data/miniboone/nu_mode/fullosc_ntuple_nu.root");
	TFile *fb =  new TFile("/home/mark/work/sbnfit_reduce/data/miniboone/nubar_mode/fullosc_ntuple_nubar.root");

	std::cout<<"Loading up TTrees"<<std::endl;
	TTree * tup = (TTree*)f->Get("fullosc_nu_mode;1");
	TTree * tupB = (TTree*)fb->Get("fullosc_nu_mode;1");


	TCanvas *c1 = new TCanvas("c1","c1",1200,400);
	c1->Divide(2,1);

	double numX = 100; 
	double numY = numX;
	double Xmin = -3.7;
	double Xmax = 0;
	double Ymin =-2;
	double Ymax = 2;
	//	Minimum chi is 16.4416409170099 at m: 0.506666666666665 and -2.61466666666666 (manual box minimization)
	//	Minimum chi is 16.341378430506 at m: 0.491999999999999 and -2.64640000000001
	//	Minimum chi is 16.3408463665959 at m: 0.493125 and -2.64587500000001
	if(false){//bf scan eh
		numX = 80;
		numY = 80;
		Xmax = -2.64;
		Xmin = -2.65;
		Ymin = 0.48;
		Ymax = 0.51;
	}






	TH2F *hcontz = new TH2F("hcontz","MicroBooNE 3+1 90\% C.L ",numX, Xmin,Xmax,numX,Ymin,Ymax  );

	hcontz->GetXaxis()->SetTitle("#sin^{2} 2 #theta_{e #mu}");
	hcontz->GetYaxis()->SetTitle("#Delta m^{2}_{41} (eV^{2})");


	std::cout<<"Loading up comb rootfiles"<<std::endl;
	//SBNspec bkg_spec("../../data/miniboone/nu_mode/uBooNE_nu_bkg", xml);
	//SBNspec obs_spec("../../data/miniboone/nu_mode/uBooNE_nu_sig", xml);
//		SBNspec bkg_spec("../../data/miniboone/nubar_mode/uBooNE_nubar_bkg", xml);
//		SBNspec obs_spec("../../data/miniboone/nubar_mode/uBooNE_nubar_sig", xml);
		SBNspec bkg_spec("../../data/miniboone/comb_mode/uBooNE_comb_bkg", xml);
		SBNspec obs_spec("../../data/miniboone/comb_mode/uBooNE_comb_sig", xml);

	std::cout<<"writing them out for posteritys"<<std::endl;
	bkg_spec.writeOut("bkg.root");
	obs_spec.writeOut("obs.root");

	TAxis *yaxis = hcontz->GetYaxis();
	TAxis *xaxis = hcontz->GetXaxis();


	double min = 99999;
	double minM = 0;
	double minS=0;
	//so varying over all Dela M and sin^2 2 theta

	double xstep = (Xmax-Xmin)/numX;
	double ystep = (Ymax-Ymin)/numY;
	obs_spec.compressVector();

	SBNspec prev_spec = bkg_spec;
	double prev_min = 1e10;

	for(int i = 0; i<5; i++){//Iteartion eh

		SBNchi * uchi;
		if(i==0){
			uchi = new SBNchi(obs_spec);
		}else{
			uchi = new SBNchi(prev_spec);
		}
		
		min = 999999;
		for(double m = Ymin; m <=Ymax; m=m+ystep){
			double tchi =0;
			for(double sins2 = Xmax ; sins2 >= Xmin; sins2 = sins2 - xstep){

				SBNspec sig_spec = bkg_spec;
				oscillateFast( pow(10,m)  , pow(10,sins2), &sig_spec.hist.at(0), tup );
			//	oscillateFast( pow(10,m)  , pow(10,sins2), &sig_spec.hist.at(0), tupB );
				oscillateFast( pow(10,m)  , pow(10,sins2), &sig_spec.hist.at(3), tupB );
				sig_spec.compressVector();
				tchi = uchi->CalcChi(&sig_spec,&obs_spec);
				//std::cout<<"Dm^2: "<<m<<" sin^2 th: "<<sins2<<" chi^2: "<<tchi<<std::endl;
				//and save wherever you like , this si just a quick hodge podge example
				hcontz->SetBinContent( 1+floor(-(Xmin-sins2)/xstep+0.00001) , 1+floor(-(Ymin-m)/ystep+0.00001), tchi);
				//hcontz->SetBinContent( xaxis->FindBin(sins2), yaxis->FindBin(m), tchi);

				if(tchi < min){min=tchi; minS = sins2; minM = m;}

			}

			std::cout<<"Dm^2: "<<m<<" chi^2: "<<tchi<<" "<<min<<std::endl;
		}
		SBNspec prev_spec = bkg_spec;
		//	oscillateFast( pow(10,minM)  , pow(10,minS), &prev_spec.hist.at(0), tupB );
			oscillateFast( pow(10,minM)  , pow(10,minS), &prev_spec.hist.at(0), tup );
			oscillateFast( pow(10,minM)  , pow(10,minS), &prev_spec.hist.at(3), tupB );
		prev_spec.compressVector();

		std::cout<<std::setprecision(15)<<"On Iteration: "<<i<<" and Minimum chi is "<<min<<" at m: "<<minM<<" and "<<minS<<std::endl;
		
		if(fabs(min-prev_min)<=0.2){
			std::cout<<"Converded on iteration : "<<i<<" prev_chi^2: "<<prev_min<<" Current_chi^2 :"<<min<<std::endl;
			break;
		}
		if(i==4){
			std::cout<<"Not Converded on iteration : "<<i<<" ending anyway: "<<min<<std::endl;
			break;
		}	

		prev_min = min;
		hcontz->Reset();

	}//end of iterations!


	//	min = 16.3408463665959;

	std::vector<double> contours = {2.3+min, 4.61+min, 5.99+min, 9.21+min};
	//std::vector<double> contours = {10,15,20,30,100};



	TFile * fout = new TFile("ans_comb.root","RECREATE");

	//	TFile * fout = new TFile("ans_nubar.root","RECREATE");
	//	TFile * fout = new TFile("ans_comb.root","RECREATE");
	fout->cd();
	hcontz->Write();
	fout->Close();


	TH2F *hcontz2 = (TH2F*)hcontz->Clone("colz");
	hcontz->SetContour(contours.size(), &contours[0]);
	c1->cd(1);
	hcontz->Draw("CONT1Z");
	c1->cd(2)->SetLogz();
	hcontz2->Draw("colz");
	c1->SaveAs("mb_comb.pdf");
	//	TFile * ff = new TFile("example_3.root","RECREATE");




	f->Close();
}

