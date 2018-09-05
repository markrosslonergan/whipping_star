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






	if(test_mode==1){

		//Load a spectra for background
		SBNspec bkg_spec("../../data/precomp/SBN_bkg_all", xml);
		std::cout<<"test_mode 1 || loaded"<<std::endl;
		//and compress down the subchannels into a channel only vector
		bkg_spec.compressVector();

		//initilize a SBNchi class with your background and covariance matrix as located in xml file
		SBNchi test_chi(bkg_spec);

		//Load up a signal, we will use the same background here
		SBNspec sig_spec("../../data/precomp/SBN_bkg_all", xml);

		//So say you think your photon-mis rate is wrong, and want to see the effects of scaling it..

		int n = 75;
		double x[n], chi[n];
		for(int i=0; i<n; i++)
		{

			double scale =0.25+i*0.025;		

			SBNspec loop_spec = sig_spec;

			//Scale the dummy spectra, scaling only histograms that match "elike_misphoton" if you used Scale.("elike") it would scale ALL elikle subchannels
			loop_spec.Scale("elike_misphoton", scale);
			loop_spec.compressVector();

			//Set these for plotting later
			x[i]=scale;
			//And calculate the chi^2
			chi[i]=test_chi.CalcChi(loop_spec);

			std::cout<<"scaling: "<<scale<<" "<<"Chi^2: "<<chi[i]<<std::endl;
		}

		//Dump your results
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
		 *	Example 2:	parameter dependant scaling
		 ************************************************************
		 ************************************************************/

	} else if (test_mode ==2){


		//Again load up a background)
		SBNspec bkg_spec("../../data/precomp/SBN_bkg_all", xml);
		//Pretend we have mis judged the muon-mis ID rate by 100% and over estimated the intrinsics by 100%
		//not exactly likely but bear with me
		bkg_spec.Scale("uBooNE_elike_mismuon", 2);
		bkg_spec.Scale("uBooNE_mlike_intrinsic", 0.5);
		bkg_spec.compressVector();

		SBNchi test_chi(bkg_spec);

		SBNspec sig_spec("../../data/precomp/SBN_bkg_all",xml);

		//And say that instead of noticing we have misjudged the muon rates, we think that our intrinis nu_e is wrong.
		//Some theorists suggests correcting the spectra by a landau function and we want to test this 

		//set up a function, can be anything
		TF1 *fLan = new TF1("fLan","TMath::Landau(x,[0],[1],0)",0,5);

		int n = 100;
		double x[n], chi2[n];

		for(int i=0; i<n; i++){
			double mpv = i*0.05;
			fLan->SetParameters(mpv, 1.3);

			//and scale the dummy specetra by the function
			//this is a HIST level so its calling the fucntion at Bin centers!
			SBNspec loop_spec = sig_spec;
			loop_spec.Scale("uBooNE_elike_intrinsic", fLan);
			loop_spec.compressVector();

			x[i]=mpv;
			chi2[i]=test_chi.CalcChi(loop_spec);

			std::cout<<"MPV: "<<mpv<<" chi^2: "<<chi2[i]<<std::endl;

		}
		//	and print out
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

		//Time for a more model dependant example
		//using just SBNspec does not give you a huge amount of power. The idea is to make your own model dependant classes from this
		//such as here, SBNosc, a SBNspec precoded to do oscillation phyiscs (at SBL)

		//Load up your background, uboone was wierdly scaled in the rootfiles so fix here
		double uboonepot=0.5;
		SBNspec bkg_spec("../../data/precomp/SBN_bkg_all", xml);
		bkg_spec.Scale("uBooNE", uboonepot);
		bkg_spec.compressVector();

		//again create a SBNchi from this spectra
		SBNchi test_chi(bkg_spec);

		//Now create a oscillation spectra, constructed the same.
		SBNosc oscSig("../../data/precomp/SBN_bkg_all",xml);

		//Say we want to just look at apperance mode (as its easiest to plot 2d!)
		oscSig.Scale("uBooNE",uboonepot);		
		oscSig.setAppMode();

		//Want to contour plot sensitivity eventually so som standard root
		TCanvas *c1 = new TCanvas("c1","c1",600,400);
		TH2F *hcontz = new TH2F("hcontz","MicroBooNE 3+1 90\% C.L ",100,-5,0, 100,-2,2);
		hcontz->GetXaxis()->SetTitle("#sin^{2} 2 #theta_{e #mu}");
		hcontz->GetYaxis()->SetTitle("#Delta m^{2}_{41} (eV^{2})");

		//so varying over all Dela M and sin^2 2 theta
		for(double m = -2.00; m <=2.04; m=m+0.04){
			for(double sins2 = 0.0 ; sins2 >= -5; sins2 = sins2 - 0.05){

				//always work in proper UPMNS elements!!
				double uei = 0.1;
				double umi = sqrt(pow(10,sins2))/(2*uei);

				//This is where you can set up 3+N
				double imn[3] = {sqrt(pow(10,m)),0,0};
				double iue[3] = {umi,0,0};
				double ium[3] = {uei,0,0};
				double iph[3] = {0,0,0};

				//construct a signalModel
				neutrinoModel signalModel(imn,iue,ium,iph);
				signalModel.numsterile = 1; //this isnt really necessary as it can tell from imn, but nice for reading

				//And load thus model into our spectra. At this point its comuted all the necessary mass-splittins and which frequencies they are
				oscSig.load_model(signalModel);

				//And apply this oscillaion! Adding to it the bkgSpec that it was initilised with.
				std::vector<double> ans = oscSig.Oscillate();

				//Then calculate a chu^2
				double tchi=test_chi.CalcChi(ans); 

				std::cout<<"Dm^2: "<<m<<" sin^2 th: "<<sins2<<" chi^2: "<<tchi<<std::endl;
				//and save wherever you like , this si just a quick hodge podge example
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
		 *		Example 3: SBNfit
		 ************************************************************
		 ************************************************************/
	} else if(test_mode ==4){

		//SBNfit is the class that does some of the fitting and links to roots minimization schemes
		//It is like SBNspec, it contains basic functionality but you probably have to extend it with our own model dependant classes. provided that the 
		// your class superseeds the virtual MinimizeCalcChi  that SBNfit has.
		//
		// It inherits from SBNchi, loading the correlation matrix in the same manner

		//set up usual signal and background
		SBNspec bkg_spec("../../data/precomp/SBN_bkg_all", xml);
		bkg_spec.compressVector();

		SBNchi test_chi(bkg_spec);

		SBNspec sig_spec("../../data/precomp/SBN_bkg_all",xml);
		sig_spec.compressVector();


		//So someone asks, lets see what normalisation shifts on elike mis-ided muon and mulike intrinics looks like, maybe we can match observations better?


		//Set up a SBNfit with background and signal, and number of paramaters
		SBNfit myfit(bkg_spec, sig_spec, 2);

		//SBNfit contains a simple routine for scaling any number of historgrams by a parameter, and minimizing over them
		//so make a vector of pairs, saying which hists to scale with which parameter
		std::vector<std::pair<std::string, int>> myin;
		myin.push_back(std::make_pair("uBooNE_elike_mismuon",0) );
		myin.push_back(std::make_pair("uBooNE_mlike_intrinsic",1) );


		myfit.initialize_norm(myin);

		//We initilise parametrs and give them names
		std::vector<double> init  {0.99,1.001};
		std::vector<double> low  {0.1,0.1};
		std::vector<double> up  {10,10};
		std::vector<double> step  {0.02,0.02};
		std::vector<std::string> nam  {"p1","p2"};

		myfit.setInitialValues(init);
		myfit.setNames(nam);
		myfit.setLowerValues(low);
		myfit.setUpperValues(up);
		myfit.setStepSizes(step);

		//and minimize, obviously as we are comparing it to itself, we hope that it find that no scaling is best!

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

		//Ok a slightly more complicated example, SBNfit on its own cant do much. So here we show a model specific example for 3+N oscillations


		//setup a SBNosc for a signal
		SBNosc inject_sig("../../data/precomp/SBN_bkg_all",xml);

		//we want to simulate  a single eV sterile
		double dmSq=1;
		double UE4=0.15;
		double UM4=0.18;
		neutrinoModel signalModel(dmSq,UE4,UM4);
		signalModel.numsterile = 1;

		//we oscillate it
		inject_sig.load_model(signalModel);
		inject_sig.OscillateThis();
		//and give it some poisson noise
		inject_sig.poissonScale();
		inject_sig.compressVector();

		//here you can see what it looks like!
		inject_sig.writeOut("example_5_injected_signal.root");


		//We then create a test point
		SBNosc test_sig("../../data/precomp/SBN_bkg_all",xml);


		//and a SBNfit3p1 (SBNfit3pN also exists and obviously can be used for n=1 also)

		SBNfit3p1  fit3p1(inject_sig, test_sig, 3);

		//Tell it what methods and algorithms we want to use
		fit3p1.setMethod("GSLMultiMin", "BFGS2");


		std::vector<std::string> nam 	{"m4", "Ue4", "Um4"};
		std::vector<double> init  	{1.0, 0.05,  0.05}; 
		std::vector<double> low  	{0.1,  0,     0};
		std::vector<double> up    	{10,   1,     1}; 
		std::vector<double> step  	{0.2,  0.05, 0.05}; 
		//As the oscillation SBNosc reads in the precomputed spectra, it is disrete in DeltaM4 (in steps of 0.04). the GSL algorithms are VERY bad at discrete minimization so 
		//here we fix the m4 variable and only minimize over the mixing elements
		std::vector<int> fix 		{1,    0,      0};
		//If you want to minimze over discrete variables, minut or GSL simulated Annealing is more successful but requires model dependant tweaking, or exhaustive search.

		fit3p1.setInitialValues(init);
		fit3p1.setNames(nam);
		fit3p1.setLowerValues(low);
		fit3p1.setUpperValues(up);
		fit3p1.setStepSizes(step);
		fit3p1.setFixed(fix);

		//minimze!
		fit3p1.Minimize();

		//As we added noise, the minimim may not be exactly the same as inputted, as expected but gives an idea of resolution
		std::cout<<"Minimized! chi^2: "<<fit3p1.bf_chi<<" Ue1: "<<fit3p1.bf_params[1]<<" Um1: "<<fit3p1.bf_params[2]<<" #:"<<fit3p1.num_func_calls<<std::endl;



		return 0;
	}else if(test_mode == 99)
	{



		//OK everything below here a bit more random and unnecessary







		//Time for a more model dependant example
		//using just SBNspec does not give you a huge amount of power. The idea is to make your own model dependant classes from this
		//such as here, SBNosc, a SBNspec precoded to do oscillation phyiscs (at SBL)

		//Load up your background, uboone was wierdly scaled in the rootfiles so fix here
		double uboonepot=2;
		SBNspec bkg_spec("../../data/precomp/SBN_bkg_all", xml);
		bkg_spec.Scale("uBooNE", uboonepot);
		bkg_spec.compressVector();

		//again create a SBNchi from this spectra
		SBNchi test_chi(bkg_spec);

		//Now create a oscillation spectra, constructed the same.
		SBNosc oscSig("../../data/precomp/SBN_bkg_all",xml);

		//Say we want to just look at apperance mode (as its easiest to plot 2d!)
		oscSig.Scale("uBooNE",uboonepot);		

		//Want to contour plot sensitivity eventually so som standard root
		TCanvas *c1 = new TCanvas("c1","c1",600,400);
		TH2F *hcontz = new TH2F("hcontz","MicroBooNE 3+1 90\% C.L ",100,-5,0, 100,-2,2);
		hcontz->GetXaxis()->SetTitle("#sin^{2} 2 #theta_{e #mu}");
		hcontz->GetYaxis()->SetTitle("#Delta m^{2}_{41} (eV^{2})");


		double th_14 = 20 ; 
		double th_24 = 10 ;
		double Ue4=sin(th_14*3.14/180.0);
		double Um4=cos(th_14*3.14/180.0)*sin(th_24*3.14/180.0);

		for(double im = -1; im<=1; im+=0.02){
			double m4 = pow(10,im);

			SBNosc aSig=oscSig;	
			SBNosc dSig=oscSig;	
			SBNosc bSig=oscSig;	

			neutrinoModel signalModel(m4,Ue4, Um4);
			signalModel.numsterile = 1; //this isnt really necessary as it can tell from imn, but nice for reading

			//And load thus model into our spectra. At this point its comuted all the necessary mass-splittins and which frequencies they are
			aSig.load_model(signalModel);
			bSig.load_model(signalModel);
			dSig.load_model(signalModel);

			aSig.setAppMode();
			dSig.setDisMode();
			bSig.setBothMode();
			//And apply this oscillaion! Adding to it the bkgSpec that it was initilised with.
			aSig.OscillateThis();
			bSig.OscillateThis();
			dSig.OscillateThis();

			//Then calculate a chu^2
			double tchi1=test_chi.CalcChi(aSig); 
			double tchi2=test_chi.CalcChi(bSig); 
			double tchi3=test_chi.CalcChi(dSig); 

			std::cout<<"th_14: "<<th_14<<" th_24: "<<th_24<<" Ue4: "<<Ue4<<" Um4: "<<Um4<<" Deltam: "<<pow(m4,2)<<" log(Dm): "<<2*im<<" Chi^2 App: "<<tchi1<<" Chi^2 Dis: "<<tchi3<<" Chi^2 Both: "<<tchi2<<std::endl;
		}
		return 0;
	}else if(test_mode==911){

		std::cout<<"1"<<std::endl;
		SBNspec bkg_spec1("../../data/precomp/SBN_bkg_all", xml);
		SBNspec bkg_spec2("../../data/precomp/SBN_bkg_all", xml);
		SBNspec bkg_spec05("../../data/precomp/SBN_bkg_all", xml);

		std::cout<<"2"<<std::endl;
		bkg_spec2.Scale("uBooNE",2);
		bkg_spec05.Scale("uBooNE",0.5);

		std::cout<<"3"<<std::endl;
		SBNchi test_chi1(bkg_spec1);
		SBNchi test_chi2(bkg_spec2);
		SBNchi test_chi05(bkg_spec05);

		std::cout<<"4"<<std::endl;
		//bkg_spec.Scale("uBooNE", uboonepot);

		std::cout<<"5"<<std::endl;
		SBNosc oscSig("../../data/precomp/SBN_bkg_all",xml);
		oscSig.setAppMode();



		//Its the full osc...Thats why


		//Want to contour plot sensitivity eventually so som standard root
		TCanvas *c1 = new TCanvas("c1","c1",600,400);
		TH2F *hcontz1 = new TH2F("hcontz1","MicroBooNE 3+1 90\% C.L ",100,-5,0, 100,-2,2);
		TH2F *hcontz2 = new TH2F("hcontz2","MicroBooNE 3+1 90\% C.L ",100,-5,0, 100,-2,2);
		TH2F *hcontz05 = new TH2F("hcontz05","MicroBooNE 3+1 90\% C.L ",100,-5,0, 100,-2,2);
		hcontz1->GetXaxis()->SetTitle("#sin^{2} 2 #theta_{e #mu}");
		hcontz05->GetYaxis()->SetTitle("#Delta m^{2}_{41} (eV^{2})");

		std::cout<<"7"<<std::endl;
		//so varying over all Dela M and sin^2 2 theta
		for(double m = -2.00; m <=2.04; m=m+0.04){
			for(double sins2 = 0.0 ; sins2 >= -5; sins2 = sins2 - 0.05){

				//always work in proper UPMNS elements!!
				double uei = 0.1;
				double umi = sqrt(pow(10,sins2))/(2*uei);

				//This is where you can set up 3+N
				double imn[3] = {sqrt(pow(10,m)),0,0};
				double iue[3] = {umi,0,0};
				double ium[3] = {uei,0,0};
				double iph[3] = {0,0,0};

				//construct a signalModel
				neutrinoModel signalModel(imn,iue,ium,iph);
				signalModel.numsterile = 1; //this isnt really necessary as it can tell from imn, but nice for reading

				//And load thus model into our spectra. At this point its comuted all the necessary mass-splittins and which frequencies they are
				oscSig.load_model(signalModel);

				//And apply this oscillaion! Adding to it the bkgSpec that it was initilised with.
				std::vector<double> ans1 = oscSig.Oscillate();
				std::vector<double> ans2 = oscSig.Oscillate(2.0);
				std::vector<double> ans05 = oscSig.Oscillate(0.5);


				//Then calculate a chu^2
				double tchi1=test_chi1.CalcChi(ans1); 
				double tchi2=test_chi2.CalcChi(ans2); 
				double tchi05=test_chi05.CalcChi(ans05); 

				//for(int i=0; i<ans1.size(); i++){
				//	std::cout<<"#@# "<<ans1[i]<<" "<<ans2[i]<<" "<<ans05[i]<<" "<<test_chi1.bkgSpec.compVec[i]<<" "<<test_chi2.bkgSpec.compVec[i]<<" "<<test_chi05.bkgSpec.compVec[i]<<std::endl;
				//}

				std::cout<<"Dm^2: "<<m<<" sin^2 th: "<<sins2<<" chi2: "<<tchi2<<" chi1: "<<tchi1<<" chi05: "<<tchi05<<std::endl;
				//and save wherever you like , this si just a quick hodge podge example
				hcontz1->SetBinContent( 1+floor(-(-5.0-sins2)/0.05+0.00001) , 1+floor(-(-2.00-m)/0.04+0.00001), tchi1);
				hcontz2->SetBinContent( 1+floor(-(-5.0-sins2)/0.05+0.00001) , 1+floor(-(-2.00-m)/0.04+0.00001), tchi2);
				hcontz05->SetBinContent( 1+floor(-(-5.0-sins2)/0.05+0.00001) , 1+floor(-(-2.00-m)/0.04+0.00001), tchi05);


			}
		}
		Double_t contours[1];
		contours[0] = 1.64;
		hcontz1->SetContour(1, contours);
		hcontz2->SetContour(1, contours);
		hcontz05->SetContour(1, contours);

		hcontz1->SetLineColor(kBlack);
		hcontz2->SetLineColor(kGreen);
		hcontz05->SetLineColor(kRed);

		c1->cd();

		hcontz1->Draw("CONT3");
		hcontz2->Draw("CONT3 SAME");
		hcontz05->Draw("CONT3 SAME");

		TFile * ff = new TFile("example_911.root","RECREATE");
		ff->cd();
		c1->Write();
		ff->Close();

		return 0;
}



	/*************************************************************
	 *************************************************************
	 *		Example 1: 	Simple loading and scaling
	 ************************************************************
	 ************************************************************/
	if(test_mode==10){
		SBNspec bkg_spec("../../build/examples/SBN_CV", xml);
		SBNspec sig_spec("../../build/examples/SBN_LEE_siggal", xml);

		sig_spec.hist[1]=bkg_spec.hist[1];	

		double kal_sig = 186.4;
		double kal_in = 362.8;

		//bkg_spec.Scale("elike", kal_in/bkg_spec.hist[0].GetSumOfWeights() );
		//sig_spec.Scale("elike", kal_in/bkg_spec.hist[0].GetSumOfWeights() );


		bkg_spec.compressVector();
		sig_spec.compressVector();




		std::cout<<"Full: "<<sig_spec.fullVec.size()<<" comp "<<sig_spec.compVec.size()<<std::endl;

		SBNchi test_chi(bkg_spec);

		std::cout<<"mat "<<test_chi.vMcI[0].size()<<" "<<test_chi.vMcI[0][0]<<std::endl;


		int n = 100;
		double x[n], chi[n];
		for(int i=0; i<n; i++){

			double scale =0+i*0.02;		

			SBNspec loop_spec = sig_spec;


			std::vector<double>excess;

			for(int k=1; k<=bkg_spec.hist[0].GetNbinsX(); k++){
				excess.push_back(bkg_spec.hist[0].GetBinContent(k)-sig_spec.hist[0].GetBinContent(k)  );
			}


			double sc=kal_in/kal_sig * bkg_spec.hist[0].GetSumOfWeights()/sig_spec.hist[0].GetSumOfWeights();

			for(int k=1; k<=bkg_spec.hist[0].GetNbinsX(); k++){
				loop_spec.hist[0].SetBinContent(k, bkg_spec.hist[0].GetBinContent(k)+ excess.at(k-1)*scale*sc);
			}


			loop_spec.compressVector();

			//Set these for plotting later
			x[i]=scale;
			//And calculate the chi^2
			chi[i]=test_chi.CalcChi(loop_spec);

			std::cout<<scale<<" "<<chi[i]<<std::endl;
		}

		//Dump your results
		TGraph *gr  = new TGraph(n,x,chi);
		TFile * ff = new TFile("lee_scale.root","RECREATE");
		TCanvas *c1 = new TCanvas("c1","Scale LEE ");

		gr->SetTitle("chi^2 of LEE off CV ");
		gr->GetXaxis()->SetTitle("Scale factor");
		gr->GetYaxis()->SetTitle("#chi^{2}");
		gr->Draw("ACP");

		c1->Write();
		ff->Close();



		return 0;


	}






}
