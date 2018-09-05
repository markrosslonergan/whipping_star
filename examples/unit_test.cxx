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
#include "TStyle.h"
#include "TLine.h"

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

TH2D * TMatrix2TH2(TMatrixT<double> * in){
	int n = in->GetNrows();
	TH2D * tmp = new TH2D("","",n,0,n,n,0,n);
	for(int i = 0; i<n; i++){
		for(int j = 0; j<n; j++){
			tmp->SetBinContent(i+1,j+1, (*in)(i,j));
		}
	}

	return tmp;
}


int TH2grid (TH2* h)
{
	int nx    = h->GetNbinsX(); 
	int ny    = h->GetNbinsY();
	double x1 = h->GetXaxis()->GetXmin();
	double x2 = h->GetXaxis()->GetXmax();
	double y1 = h->GetYaxis()->GetXmin();
	double y2 = h->GetYaxis()->GetXmax();
	double x=x1;
	double y=y1;
	TLine l; 
	for (int i=0; i<=nx; i++) {
		l.SetLineWidth(1);
		l.DrawLine(x,y1,x,y2);   
		x = x+h->GetXaxis()->GetBinWidth(i);
	}
	for (int j=0; j<=ny; j++) {
		l.SetLineWidth(1);
		l.DrawLine(x1,y,x2,y); 
		y = y+h->GetYaxis()->GetBinWidth(j);
	}
}

/*************************************************************
 *************************************************************
 *		BEGIN unit_tests.cxx
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
	std::cout<<"Open new TFile to output results."<<std::endl;
	TFile   *file = new TFile("unit_test.root","RECREATE");
	gStyle->SetOptStat(kFALSE);


	/*************************************************************
	 *************************************************************
	 *		Unit Test 1:
	 ************************************************************
	 ************************************************************/
	std::cout<<"#UNIT_TEST_1:: Begininig unit test 1."<<std::endl;
	if(test_mode == 1 || test_mode == 0){
		TCanvas *c1 = new TCanvas("unit_test_1","unit_test_1",2000,1000);
		c1->Divide(2,1);
		TPad * p1 = (TPad*)c1->cd(1);
		//c1->SetGrid(1,1);

		std::cout<<"#UNIT_TEST_1:: Constructing SBNspec for passed xml, 'unit.xml'"<<std::endl;
		SBNspec unit_spec(xml);
		unit_spec.compressVector();

		std::cout<<"#UNIT_TEST_1:: Creating a toy flat covariance matrix, 0.1 in each entry for nu"<<std::endl;
		TMatrixT<double> flat_covar(unit_spec.num_bins_total, unit_spec.num_bins_total);

		flat_covar = 0.1;

		std::cout<<"#UNIT_TEST_1:: Drawing said toy covariance matrix to the canvas"<<std::endl;
		TH2D * tmp = TMatrix2TH2(&flat_covar);
		//p1->SetPhi(0.0001);
		//p1->SetTheta(90);
		tmp->Draw("colz");
		TH2grid(tmp);	

		std::cout<<"#UNIT_TEST_1:: Constructing a SBNchi object from SBNspec and toy covariance"<<std::endl;
		SBNchi unit_chi(unit_spec, flat_covar);

		TPad* p2 = (TPad*)c1->cd(2);
		//p2->SetPhi(0.0001);
		//p2->SetTheta(90);

		TMatrixT<double> flat_covar_compressed(unit_spec.num_bins_total_compressed, unit_spec.num_bins_total_compressed); 
		unit_chi.collapse_layer3(flat_covar, flat_covar_compressed);


		//TMatrix2TH2(&flat_covar_compressed)->Draw("lego2");
		TH2D * tmp2 = TMatrix2TH2(&flat_covar_compressed);
		tmp2->Draw("colz");
		TH2grid(tmp2);
		//TMatrix2TH2(&flat_covar_compressed)->Draw("box same");

		c1->Write();
		c1->SaveAs( "unit_1.pdf","pdf");

		std::cout<<"#UNIT_TEST_1:: End of unit test 1."<<std::endl;
	}

	/*************************************************************
	 *************************************************************
	 *		Unit Test 2:
	 ************************************************************
	 ************************************************************/
	std::cout<<"#UNIT_TEST_2:: Begininig unit test 2."<<std::endl;

	if (test_mode == 2 || test_mode ==0){
		std::cout<<"#UNIT_TEST_2:: Constructing SBNspec for passed xml."<<std::endl;
		SBNspec unit_spec(xml);
		//unit_spec.setAsGaussian(1,0.2, 10000);	
		unit_spec.setAsFlat(100);	
		unit_spec.compressVector();
		unit_spec.writeOut("test2.root");

		SBNspec unit_spec_inf(xml);
		unit_spec_inf.setAsFlat(1000000000);
		unit_spec_inf.compressVector();



		std::cout<<"#UNIT_TEST_2:: Creating a toy flat covariance matrix, 0.1 in each entry for nu"<<std::endl;
		TMatrixT<double> flat_covar(unit_spec.num_bins_total, unit_spec.num_bins_total);
		flat_covar = 0.0;

		SBNchi unit_chi_stat_only(unit_spec, flat_covar);

		flat_covar = 0.1*0.1;
		SBNchi unit_chi1(unit_spec_inf, flat_covar);

		flat_covar = 0.05*0.05;
		SBNchi unit_chi2(unit_spec_inf, flat_covar);

		flat_covar = 0.2*0.2;
		SBNchi unit_chi3(unit_spec_inf, flat_covar);


		int n = 50;
		double x[n], chi1[n], chi2[n],chi3[n], chi_stat_only[n];
		for(int i=0; i<n; i++)
		{

			SBNspec loop_spec = unit_spec;
			SBNspec loop_spec_inf = unit_spec_inf;

			double scale =0.5+i*0.02;		

			loop_spec.Scale("elike", scale);
			loop_spec.compressVector();


			loop_spec_inf.Scale("elike", scale);
			loop_spec_inf.compressVector();

			//Set these for plotting later
			x[i]=scale;
			//And calculate the chi^2
			chi_stat_only[i]=unit_chi_stat_only.CalcChi(loop_spec)/20.0;
			chi1[i]=unit_chi1.CalcChi(loop_spec_inf);
			chi2[i]=unit_chi2.CalcChi(loop_spec_inf);
			chi3[i]=unit_chi3.CalcChi(loop_spec_inf);

			std::cout<<"UNIT_TEST_2:: Scaling: "<<scale<<" "<<"Chi^2: "<<chi1[i]<<" "<<chi2[i]<<" "<<chi3[i]<<" Chi^2 Stat Only: "<<chi_stat_only[i]<<std::endl;
		}

		//Dump your results
		TGraph *gr1  = new TGraph(n,x,chi1);
		TGraph *gr2  = new TGraph(n,x,chi2);
		TGraph *gr3  = new TGraph(n,x,chi3);
		TGraph *gr_stat_only  = new TGraph(n,x,chi_stat_only);

		//Plotting here	
		TCanvas *c = new TCanvas("unit_test_2","unit_test_2",2000,1000);
		c->Divide(2,1);

		TPad * p1 = (TPad*)c->cd(1);
		double ymax = 2;
		gr_stat_only->SetTitle("Scaling of flat 100 event distribution (Stat Only, 20 dof)");
		gr_stat_only->GetXaxis()->SetTitle("Scale Factor");
		gr_stat_only->GetYaxis()->SetTitleOffset(1.3);
		gr_stat_only->GetYaxis()->SetTitle("#chi^{2} / ndof");
		gr_stat_only->GetHistogram()->SetMaximum(ymax);
		gr_stat_only->GetXaxis()->SetLimits(0.8,1.2);
		gr_stat_only->SetLineColor(kBlue-7); 
		gr_stat_only->SetLineWidth(3);
		TLine lv; lv.SetLineColor(kRed-6); lv.SetLineWidth(3);
		TLine lh; lh.SetLineColor(kRed-6); lh.SetLineWidth(3);
		gr_stat_only->Draw("ACP");
		lv.DrawLine(0.9,0,0.9,1.0);lv.DrawLine(1.1,0,1.1,1.0);
		lh.DrawLine(0.8,1,1.2,1);


		TPad * p2 = (TPad*)c->cd(2);
		gr1->SetTitle("Scaling of flat event distribution (Sys only, 1 dof). Fully correlated, #sigma^{2} = 0.1");
		gr1->GetXaxis()->SetTitle("Scale Factor");
		gr1->GetYaxis()->SetTitleOffset(1.3);
		gr1->GetYaxis()->SetTitle("#chi^{2} / ndof");
		gr1->GetXaxis()->SetLimits(0.8,1.2);
		gr1->GetHistogram()->SetMaximum(ymax);

		gr1->SetLineWidth(3);gr1->SetLineColor(kGreen-6);
		gr2->SetLineWidth(3);gr2->SetLineColor(kRed-7);
		gr3->SetLineWidth(3);gr3->SetLineColor(kGreen-6);

		gr1->Draw("ACP");
		//	gr2->Draw("same CP");
		//	gr3->Draw("same CP");

		lv.DrawLine(0.9,0,0.9,1.0);lv.DrawLine(1.1,0,1.1,1.0);
		lh.DrawLine(0.8,1,1.2,1);


		file->cd();
		c->Write();




	}// end of unit_test 2 


	/*************************************************************
	 *************************************************************
	 *		Unit Test 3:
	 ************************************************************
	 ************************************************************/
	std::cout<<"#UNIT_TEST_2:: Begininig unit test 2."<<std::endl;

	if (test_mode == 3 || test_mode ==0){
		std::cout<<"#UNIT_TEST_3:: Constructing SBNspec for passed xml."<<std::endl;
		SBNspec unit_spec(xml);
		unit_spec.setAsGaussian(1,0.5, 10000000);	
		unit_spec.compressVector();


		std::cout<<"#UNIT_TEST_3:: Creating a toy flat covariance matrix, 0.1 in each entry for nu"<<std::endl;
		TMatrixT<double> flat_covar(unit_spec.num_bins_total, unit_spec.num_bins_total);
		flat_covar = 0.0;

		SBNchi unit_chi_stat_only(unit_spec, flat_covar);


		TH1D * h_chi_stat_only= new TH1D("h_chi_stat_only","",20,0,2);

		int n = 100;
		double x[n], chi1[n], chi2[n],chi3[n], chi_stat_only[n];
		for(int i=0; i<n; i++)
		{
			double ex = (double)i*7.0/((double)n);
			int fac = pow(10,ex);

			SBNspec loop_spec = unit_spec;
			loop_spec.setAsGaussian(1,0.5,fac);
			loop_spec.NormAll(unit_spec.hist.at(0).GetSumOfWeights());
			loop_spec.compressVector();


			std::cout<<"UNIT_TEST_3:: sum of w: "<<unit_spec.hist.at(0).GetSumOfWeights()<<" "<<loop_spec.hist.at(0).GetSumOfWeights()<<std::endl;

			x[i]=fac;
			chi_stat_only[i]=unit_chi_stat_only.CalcChi(loop_spec);

			std::cout<<"UNIT_TEST_3:: n "<<fac<<" Chi^2 Stat Only: "<<chi_stat_only[i]<<std::endl;
			h_chi_stat_only->Fill(chi_stat_only[i]);
		}

		//Plotting here	
		TCanvas *c = new TCanvas("unit_test_3","unit_test_3",2000,1000);
		c->Divide(2,1);

		TPad * p1 = (TPad*)c->cd(1);
		TGraph *gr_stat_only  = new TGraph(n,x,chi_stat_only);
		gr_stat_only->Print();
			p1->SetLogy();
			p1->SetLogx();
		double ymax = 2;
		gr_stat_only->SetTitle("Spectral fit to a gaussian ");
		gr_stat_only->GetXaxis()->SetTitle("Number of generated events in test gaussian");
		gr_stat_only->GetYaxis()->SetTitleOffset(1.3);
		gr_stat_only->GetYaxis()->SetTitle("#chi^{2} / ndof");
		//	gr_stat_only->GetHistogram()->SetMaximum(ymax);
		gr_stat_only->GetXaxis()->SetLimits(10,1000000);
		gr_stat_only->SetLineColor(kBlue-7); 
		gr_stat_only->SetLineWidth(3);
		TLine lv; lv.SetLineColor(kRed-6); lv.SetLineWidth(3);
		TLine lh; lh.SetLineColor(kRed-6); lh.SetLineWidth(3);
		gr_stat_only->Draw("AlP");
		//	lv.DrawLine(0.9,0,0.9,1.0);lv.DrawLine(1.1,0,1.1,1.0);
		//	lh.DrawLine(0.8,1,1.2,1);


		TPad * p2 = (TPad*)c->cd(2);


		file->cd();
		c->Write();




	}





	file->Close();
}
