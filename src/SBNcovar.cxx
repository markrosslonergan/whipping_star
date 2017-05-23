#include "SBNcovar.h"
//#include "MCEventWeight.h"

using namespace sbn;

SBNcovar::SBNcovar(std::string rootfile, std::string xmlname) : SBNconfig(xmlname) {
	std::cout<<"inside SBNcovar"<<std::endl;	
	gROOT->ProcessLine("#include <map>");
	gROOT->ProcessLine("#include <vector>");
	gROOT->ProcessLine("#include <string>");
	gSystem->Load("/uboone/app/users/markrl/sbnfit/whipping_star/src/mdict_h.so");
	gStyle->SetOptStat(0);

	tolerence_positivesemi = 1e-10;
	is_small_negative_eigenvalue = false;
	bool DEBUG_BOOL = false;

	//Initialise all the things
	//for every multisim, create a SBNspec
	for(int m=0; m<num_multisim; m++){
		SBNspec temp_spec(xmlname,m);
		multi_hists.push_back(temp_spec);
	}
	SBNspec tm(xmlname,-1);
	spec_CV = tm;

	TFile *f = new TFile(rootfile.c_str());
	TTree * full_mc=  (TTree*)f->Get(multisim_name.c_str());
	std::cout<<"Got tree"<<std::endl;
	
	std::vector<double> weights;

	std::vector<int> vars_i (branch_names_int.size(),0); 
	std::vector<double> vars_d (branch_names_double.size(),0);
	std::map<std::string, std::vector<double> > *fWeight = 0;


	TBranch *bweight=0;
	full_mc->SetBranchAddress("mcweight", &fWeight, &bweight);		


	std::cout<<"Branch set"<<std::endl;

	for(int i=0; i< branch_names_int.size(); i++){
				std::cout<<"Setting int: "<<branch_names_int[i]<<std::endl;
				full_mc->SetBranchAddress(branch_names_int[i].c_str(), &(vars_i[i]) );

		}
	for(int i=0; i< branch_names_double.size(); i++){
				std::cout<<"Setting double: "<<branch_names_double[i]<<std::endl;
				full_mc->SetBranchAddress(branch_names_double[i].c_str(), &(vars_d[i]) );

		}


		 int nentries = full_mc->GetEntries(); // read the number of entries

		 //	nentries = 50000;
		 // So for every entry..
		 std::cout<<"Starting loop over entries: "<<nentries<<std::endl;

		 //for (int i = 50000; i <  60000 ; i++) {
		 for (int i = 0; i <  nentries ; i++) {
			int mm=0;
			weights.clear(); 
			 double glob_corr=1;
			 
			if(i%1000==0) std::cout<<"On entry :"<<i<<" over entries: "<<nentries<<std::endl;

		 	if(DEBUG_BOOL) std::cout<<"On entry :"<<i<<" over entries: "<<nentries<<std::endl;
			 full_mc->GetEntry(i);
		
			//if(vars_i[2] != 12 || vars_i[3] != 11) continue;
			if(vars_i[0] != 1001 ) continue;
	
			 for(int h=0;h<branch_names_double.size();h++){
				if(DEBUG_BOOL)std::cout<<h<<" "<<branch_names_double[h].c_str()<<" "<<vars_d[h]<<std::endl;
			 }
	
			 for(int h=0;h<branch_names_int.size();h++){
				if(DEBUG_BOOL)std::cout<<h<<" "<<branch_names_int[h].c_str()<<" "<<vars_i[h]<<std::endl;
			 }

			
			for(std::map<std::string, std::vector<double> >::iterator  it = fWeight->begin(); it != fWeight->end(); ++it) {
				if(DEBUG_BOOL)std::cout << it->first <<" "<<it->second.size()<<std::endl;
				
				if(it->first == "bnbcorrection_FluxHist"){
					if(DEBUG_BOOL)std::cout<<"Got the bnb flux correction"<<std::endl;	
					glob_corr = it->second.at(0);		
					continue;
				}
				for(auto &j: it->second){

					if(i==1529 || i == 1629)continue;
					if(DEBUG_BOOL)std::cout<<i<<" Weight: "<<j<<" d: "<<mm<<std::endl;
						double wei = j;
						if(wei==0) {
							wei= 0.0000001;
							if(DEBUG_BOOL)std::cout<<"ERROR: Weight 0 on event: "<<i<<" Weight: "<<j<<" multisim: "<<mm<<" "<<it->first<<" size: "<<it->second.size()<<std::endl;
						}
						weights.push_back(wei*glob_corr);
						mm++;
				}
				
			}


			 // loop over every multisim we have in this entry
			 // the first must be 1, i.e the CV or mean case.
			double num_sim = weights.size();	 

			//if(num_sim != multi_hists.size()){
			//	std::cout<<"ERROR: numsim!= nulti_hist"<<num_sim<<" "<<multi_hists.size()<<std::endl;
				//exit(EXIT_FAILURE);
			//}

			if(DEBUG_BOOL)std::cout<<"Num Sim "<<num_sim<<std::endl;
			for(int m=0; m<num_sim; m++){
				//std::cout<<"On Uni :"<<m<<" of "<<num_sim<<std::endl;
				//Then in this sim, fill every histogram! 
				for(auto & h: multi_hists[m].hist){	

				if(vars_d[0]!=vars_d[0] || weights[m]!= weights[m]){
				std::cout<<"ERROR: vars_d[0]: "<<vars_d[0]<<" weights: "<<weights[m]<<" m: "<<m<<std::endl;
				//exit(EXIT_FAILURE);
				//continue;

				}


				h.Fill(vars_d[0],weights[m]);
				}
			 }
			
				for(auto & h: spec_CV.hist){	
					h.Fill(vars_d[0],glob_corr);
				}
			
			
		 } //end of entry loop

		f->Close();

	}



	int SBNcovar::formCovarianceMatrix(){

		full_covariance.ResizeTo(num_bins_total, num_bins_total);
		frac_covariance.ResizeTo(num_bins_total, num_bins_total);
		full_correlation.ResizeTo(num_bins_total, num_bins_total);

		TH2D * h2 = new TH2D("Frac Cov","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);
		TH2D * h3 = new TH2D("Corr","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);
		TH2D * h4 = new TH2D("Full Cov","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);
		

		std::cout<<"Starting formCovariance: we have "<<multi_hists.size()<<" histograms "<<std::endl;
		for(int i=0; i< multi_hists.size(); i++){
			multi_hists[i].calcFullVector();
		}
	
		spec_CV.calcFullVector();	
		std::vector<double> CV = spec_CV.fullVec;


		for(int i=0; i<num_bins_total; i++){
		  for(int j=0; j<num_bins_total; j++){
			
			full_covariance(i,j)=0;

			//Remember that m=0 is CV!
			for(int m=1; m < multi_hists.size(); m++){
				full_covariance(i,j) += (CV[i]-multi_hists[m].fullVec[i])*(CV[j]-multi_hists[m].fullVec[j]);
			

				if(full_covariance(i,j)!=full_covariance(i,j)){
					
				std::cout<<"ERROR: nan : at (i,j):  "<<i<<" "<<j<<" fullcov: "<<full_covariance(i,j)<<" multi hist sise "<<multi_hists.size()<<" CV: "<<CV[i]<<" "<<CV[j]<<" multihisg "<<multi_hists[m].fullVec[i]<<" "<<multi_hists[m].fullVec[j]<<" on dim : "<<m<<std::endl;
					
				}



			}
			full_covariance(i,j) = full_covariance(i,j)/( multi_hists.size()-1.0);
			
		  }
		}
		std::cout<<"Final Matrix"<<std::endl;
	
	
		for(int i=0; i<num_bins_total; i++){
		  for(int j=0; j<num_bins_total; j++){
			std::cout<<i<<" "<<j<<" "<<full_covariance(i,j)<<std::endl;

frac_covariance(i,j) = full_covariance(i,j)/(spec_CV.fullVec[i]*spec_CV.fullVec[j]) ;
full_correlation(i,j)= full_covariance(i,j)/(sqrt(full_covariance(i,i))*sqrt(full_covariance(j,j)));

			h2->SetBinContent(i+1,j+1,frac_covariance(i,j));
			h3->SetBinContent(i+1,j+1,full_correlation(i,j));
			h4->SetBinContent(i+1,j+1,full_covariance(i,j));

		  }
		}

		/************************************************************
		 *		Quality Testing Suite			    *
		 * *********************************************************/

		TFile *ftest=new TFile("qtest.root","RECREATE");
		ftest->cd();
		TCanvas *c1 =  new TCanvas();
		c1->cd();
		//full_covariance.Draw();
		h2->SetTitle("Fractional Covariance Matrix (sys only)");
		h2->GetYaxis()->SetTitle("E_{#nu}^{truth}");
		h2->GetXaxis()->SetTitle("E_{#nu}^{truth}");
		h2->Draw("COLZ");
		c1->Write();

		TCanvas *c2 =  new TCanvas();
		c2->cd();

		h3->SetTitle("Correlation Matrix (sys only)");
		h3->GetYaxis()->SetTitle("E_{#nu}^{truth}");
		h3->GetXaxis()->SetTitle("E_{#nu}^{truth}");
		h3->Draw("COLZ");
		c2->Write();	

		TCanvas *c3 =  new TCanvas();
		c3->cd();

		h4->SetTitle("Covariance Matrix (sys only)");
		h4->GetYaxis()->SetTitle("E_{#nu}^{truth}");
		h4->GetXaxis()->SetTitle("E_{#nu}^{truth}");
		h4->Draw("COLZ");
		c3->Write();
	

		full_covariance.Write();
		ftest->Close();


		spec_CV.writeOut("CV.root");




	

		if(full_covariance.IsSymmetric()){
			std::cout<<"Generated covariance matrix is symmetric"<<std::endl;
		}else{
			std::cerr<<"ERROR: SBNcovar::formCovarianceMatrix, result is not symmetric!"<<std::endl;
			exit(EXIT_FAILURE);
		}


		//if a matrix is (a) real and (b) symmetric (checked above) then to prove positive semi-definite, we just need to check eigenvalues and >=0;
		TMatrixDEigen eigen (full_covariance);
		TVectorD eigen_values = eigen.GetEigenValuesRe();


		for(int i=0; i< eigen_values.GetNoElements(); i++){
			if(eigen_values(i)<0){
				is_small_negative_eigenvalue = true;
				if(fabs(eigen_values(i))> tolerence_positivesemi ){
					std::cerr<<"ERROR: SBNcovar::formCovarianceMatrix, contains (at least one)  negative eigenvalue: "<<eigen_values(i)<<std::endl;
					exit(EXIT_FAILURE);
				}
			}
		}


		if(is_small_negative_eigenvalue){	
			std::cout<<"Generated covariance matrix is (allmost) positive semi-definite. It did contain small negative values of absolute value <= :"<<tolerence_positivesemi<<std::endl;
		}else{
			std::cout<<"Generated covariance matrix is also positive semi-definite."<<std::endl;
		}


/*
		for (int i=1;i<=11 ;i++){
		       	h2->GetYaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
		       	h2->GetXaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
		       	h4->GetYaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
		       	h4->GetXaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
		       	h3->GetYaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
		       	h3->GetXaxis()->SetBinLabel(i,std::to_string( bin_edges[0][i-1] ).c_str());
	}
*/

	return 0;
	}

