#include "SBNcovar.h"
//#include "MCEventWeight.h"

using namespace sbn;

SBNcovar::SBNcovar(std::string rootfile, std::string xmlname) : SBNconfig(xmlname) {
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


	std::cout<<" multi_hist_size: "<<multi_hists.size()<<" internal size: "<<multi_hists.at(0).hist.size()<<"  deeper: "<<multi_hists.at(0).hist.at(0).GetBinContent(1)<<std::endl;

	SBNspec tm(xmlname,-1);
	spec_CV = tm;

	//Load all files as per xml
	std::vector<TFile *> files;	
	std::vector<TTree *> trees;	
	
	for(auto &fn: multisim_file){
		files.push_back(new TFile(fn.c_str()));
	}

	for(int i=0; i<multisim_name.size(); i++){
		trees.push_back((TTree*)files.at(i)->Get(multisim_name[i].c_str()) );
	}

	

	std::vector<int> vars_i (branch_names_int.size(),0); 
	std::vector<double> vars_d (branch_names_double.size(),0);
	std::map<std::string, std::vector<double> > *fWeight = 0;

	std::vector<int> vars_i2 (branch_names_int.size(),0); 
	std::vector<double> vars_d2 (branch_names_double.size(),0);
	std::map<std::string, std::vector<double> > *fWeight2 = 0;


	TBranch *bweight=0;
	TBranch *bweight2 =0;
	full_mc->SetBranchAddress("mcweight", &fWeight, &bweight);		
	full_mc2->SetBranchAddress("mcweight", &fWeight2, &bweight2);		

	TFile *fS = new TFile("/uboone/app/users/mastbaum/leerw/fittree_lee_rw.root");
	TTree * tnS =  (TTree*)fS->Get("leerw;1");


	float lee=0;
	tnS->SetBranchAddress("leerw",&lee);
	SBNspec spec_sig(xmlname,-1);








	std::cout<<"Branch set"<<std::endl;

	for(int i=0; i< branch_names_int.size(); i++){
		std::cout<<"Setting int: "<<branch_names_int[i]<<std::endl;
		full_mc->SetBranchAddress(branch_names_int[i].c_str(), &(vars_i[i]) );
		full_mc2->SetBranchAddress(branch_names_int[i].c_str(), &(vars_i2[i]) );
	}
	for(int i=0; i< branch_names_double.size(); i++){
		std::cout<<"Setting double: "<<branch_names_double[i]<<std::endl;
		full_mc->SetBranchAddress(branch_names_double[i].c_str(), &(vars_d[i]) );
		full_mc2->SetBranchAddress(branch_names_double[i].c_str(), &(vars_d2[i]) );
	}


	int nentries = full_mc->GetEntries(); // read the number of entries
	int nentries2 = full_mc2->GetEntries(); // read the number of entries

	std::cout<<"Starting loop over entries: numu#: "<<nentries<<" numu#: "<<full_mc2->GetEntries()<<std::endl;
	
	std::cout<<"Init facts: multihistsize: "<<multi_hists.size()<<std::endl;
	std::cout<<"start with mu"<<std::endl;

	


	for (int i2 = 0; i2 <  nentries2 ; i2++) {
		bool isValid2 =true;
		int mm2 =0;
		std::vector<double> weights2;
		double glob_corr2 = 1;

		if(i2%10000==0)std::cout<<"On entry :"<<i2<<" over entries numu: "<<nentries2<<std::endl;
		full_mc2->GetEntry(i2);

		if(vars_i2[0] == 1001){

			//Loop over all things in mcweight map
			for(std::map<std::string, std::vector<double> >::iterator  it = fWeight2->begin(); it != fWeight2->end(); ++it) {
				
				if(!isValid2) break;

				if(it->first == "bnbcorrection_FluxHist"){
					glob_corr2 = it->second.at(0);		

					if(std::isinf(glob_corr2)){
						isValid2=false;
						break;
					}
					
					continue;
				}

				for(auto &j: it->second){

				//	std::cout<<i<<" Weight: "<<j<<" d: "<<mm<<std::endl;
					double wei = j*glob_corr2;

					if(std::isnan(wei) || (wei!= wei) || std::isinf(glob_corr2) ){
						std::cout<<"ERROR: "<<" weight has a value of: "<<j<<" and glob_corr: "<<glob_corr2<<" this wei: "<<wei<<". So I am skipping this event #: "<<i2<<std::endl;
						isValid2 = false;
						break;
					}

					weights2.push_back(wei);
					mm2++;

				}

			}//end of iterator


			if(isValid2){
				// loop over every multisim we have in this entry
				// the first must be 1, i.e the CV or mean case.
				double num_sim = weights2.size();	 
				
				if(num_sim > multi_hists.size()){
					std::cout<<"ERROR: numu loop numsim > nulti_hist"<<num_sim<<" "<<multi_hists.size()<<std::endl;
					num_sim = multi_hists.size();
				}

				for(int m=0; m<num_sim; m++){
						
						multi_hists[m].hist.at(1).Fill(vars_d2[0],weights2[m]);

						//important check. failure mode	
						if(weights2[m]!=weights2[m] || isinf(weights2[m]) ){
							std::cout<<"ERROR: weight has a value of: "<<weights2[m]<<". So I am killing all. on Dim: "<<m<<" energy"<<vars_d2[0]<<" global_corr is "<<glob_corr2<<std::endl;
							exit(EXIT_FAILURE);
						}



				}

					spec_CV.hist.at(1).Fill(vars_d2[0],glob_corr2);


			}//end of validity if



		}//end of mu like

	} //end of entry loop



	std::cout<<"Starting nue"<<std::endl;
	for (int i = 0; i < nentries ; i++) {
		bool isValid=true;
		int mm=0;
		std::vector<double> weights;
		double glob_corr = 1;

		if(i%10000==0)std::cout<<"On entry nue :"<<i<<" over entries: "<<nentries<<std::endl;

		full_mc->GetEntry(i);
		tnS->GetEntry(i);

		if(vars_i[0] == 1001){

			for(std::map<std::string, std::vector<double> >::iterator  it = fWeight->begin(); it != fWeight->end(); ++it) {

				if(!isValid) break;
				if(it->first == "bnbcorrection_FluxHist"){
					glob_corr = it->second.at(0);		
					
					if(std::isinf(glob_corr)){
						isValid=false;
						break;
					}


					continue;
				}

				for(auto &j: it->second){

					
					double wei = j*glob_corr;

					if(std::isnan(wei) || (wei!= wei) || std::isinf(glob_corr) ){
						std::cout<<"ERROR: "<<" weight has a value of: "<<j<<" and glob_corr: "<<glob_corr<<" this wei: "<<wei<<". So I am skipping this event #: "<<i<<std::endl;
						isValid = false;
						break;
					}


					weights.push_back(wei);
					mm++;


				}

			}//end of iterator
			if(isValid){
				// loop over every multisim we have in this entry
				// the first must be 1, i.e the CV or mean case.
				double num_sim = weights.size();	 
				if(num_sim > multi_hists.size()){
					std::cout<<"ERROR: nue loop: numsim> multi_hists: numsim: "<<num_sim<<" multihi"<<multi_hists.size()<<std::endl;
					num_sim = multi_hists.size();
				}

				for(int m=0; m<num_sim; m++){
						//std::cout<<"ATTN: "<<m<<" "<<multi_hists[m].hist.size()<<" nummultisim_xml: "<<num_multisim<<" where as weightts: "<<weights.size()<<" "<<num_sim<<std::endl;
						multi_hists[m].hist.at(0).Fill(vars_d[0],weights[m]);

						if(weights[m]!=weights[m] || std::isinf(weights[m])){
							std::cout<<"ERROR: weight has a value of: "<<weights[m]<<". So I am killing all. on Dim: "<<m<<" energy"<<vars_d[0]<<" global_corr is "<<glob_corr<<std::endl;
							exit(EXIT_FAILURE);
						}


				}
					spec_CV.hist[0].Fill(vars_d[0],glob_corr);

					spec_sig.hist[0].Fill(vars_d[0],glob_corr*(1+lee));

			}//end of validity if



		}//end of nue file



	} //end of entry loop

	spec_sig.calcFullVector();
	spec_sig.writeOut("leesig.root");

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

		for(auto &h: multi_hists){
			h.calcFullVector();
		}

		spec_CV.calcFullVector();	
		std::vector<double> CV = spec_CV.fullVec;

		for(auto k=0;k<spec_CV.fullVec.size();k++){
			std::cout<<" "<<spec_CV.fullVec[k]<<" "<<multi_hists.at(0).fullVec[k]<<std::endl;
		}

		for(int i=0; i<num_bins_total; i++){
			for(int j=0; j<num_bins_total; j++){

				full_covariance(i,j)=0;

				for(int m=0; m < multi_hists.size(); m++){
					full_covariance(i,j) += (CV[i]-multi_hists.at(m).fullVec.at(i))*(CV[j]-multi_hists.at(m).fullVec.at(j));


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


		frac_covariance.Write();
		full_covariance.Write();
		full_correlation.Write();
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

