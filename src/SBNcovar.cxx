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
		multi_sbnspec.push_back(temp_spec);
	}


	std::cout<<" multi_hist_size: "<<multi_sbnspec.size()<<" internal size: "<<multi_sbnspec.at(0).hist.size()<<"  deeper: "<<multi_sbnspec.at(0).hist.at(0).GetBinContent(1)<<std::endl;

	SBNspec tm(xmlname,-1);
	spec_CV = tm;

	//Load all files as per xml
	std::vector<TFile *> files;	
	std::vector<TTree *> trees;	

	int Nfiles = multisim_file.size();

	for(auto &fn: multisim_file){
		files.push_back(new TFile(fn.c_str()));
	}

	for(int i=0; i<multisim_name.size(); i++){
		trees.push_back((TTree*)files.at(i)->Get(multisim_name[i].c_str()) );
	}

	std::vector<int> nentries;
	for(auto &t: trees){
		nentries.push_back(t->GetEntries());
	}

	
	
	vars_i(Nfiles   , std::vector<int>(branch_names_int.size(),0));
	vars_d(Nfiles   , std::vector<double>(branch_names_double.size(),0.0));

	std::vector< std::map<std::string, std::vector<double> > * > fWeights(multisim_name.size(),0);

	std::vector< TBranch *> bweight(Nfiles,0);

	for(int i=0; i< Nfiles; i++){
		trees.at(i)->SetBranchAddress("mcweight", &(fWeights.at(i)), &(bweight.at(i)) );

		for(auto &bni: branch_names_int){
			trees.at(i)->SetBranchAddress(bni.c_str(), &(vars_i.at(i)));
		}
		for(auto &bnd: branch_names_double){
			trees.at(i)->SetBranchAddress(bnd.c_str(), &(vars_d.at(i)));
		}
	}



	//signal crap todo
	TFile *fS = new TFile("/uboone/app/users/mastbaum/leerw/fittree_lee_rw.root");
	TTree * tnS =  (TTree*)fS->Get("leerw;1");
	float lee=0;
	tnS->SetBranchAddress("leerw",&lee);
	SBNspec spec_sig(xmlname,-1);


	for(int j=0;j<Nfiles;j++){
		for(int i=0; i< nentries.at(j); i++){
			bool is_valid = true;
			int n_uni = 0;
			std::vector<double> weights = 1;
			double global_weight = 1;

			if(i%2500==0)std::cout<<"Event: "<<j<<" of "<<nentries[j]<<" from File: "<<multisim_file[j]<<std::endl;

			trees.at(j)->GetEntry(i);

			//here we put the selection criteria, for example nuance interaction 1001 == CCQE, virtual bool
			if( this->eventSelection(j) ){

				//Loop over all things in mcweight map
				for(std::map<std::string, std::vector<double> >::iterator  it = fWeight.at(j)->begin(); it != fWeight.at(j)->end(); ++it) {

					if(!is_valid) break;

					if(it->first == "bnbcorrection_FluxHist"){
						global_weight = it->second.at(0);		

						if(std::isinf(global_weight)){
							is_valid=false;
							break;
						}

						continue;
					}



					for(auto &wei: it->second){

						//TODO: here have a SBNconfig variable where we can either have "all" for all variables
						// or a list of strings containing the factors


						if(std::isnan(wei) || std::isinf(wei) || std::isinf(global_weight) || std::isnan(global_weight)  ){
							std::cout<<"ERROR: "<<" weight has a value of: "<<wei<<" and global_weight: "<<global_weight;
							std::cout<<". So I am skipping this event #: "<<i<<std::endl;
							is_valid = false;
							break;
						}

						weights.at(j).push_back(wei*global_weight);
						n_uni++;

					}

				}//end of weight (it->second) iterator


				if(is_valid){
					// loop over every multisim we have in this entry
					double num_sim = weights.at(j).size();	 

					if(num_sim > multi_sbnspec.size()){
						std::cout<<"ERROR: numu loop numsim > nulti_hist"<<num_sim<<" "<<multi_sbnspec.size()<<std::endl;
						num_sim = multi_sbnspec.size();
					}

					for(int m=0; m<num_sim; m++){
						//This is the part where we will every histogram in this Universe
						this->fillHistograms(j, m, weights.at(j).at(m) );

						//important check. failure mode	
						if(weights.at(j)[m]!=weights.at(j)[m] || isinf(weights.at(j)[m]) ){
							std::cout<<"ERROR: weight has a value of: "<<weights.at(k)[m]<<". So I am killing all. on Dim: "<<m<<" energy"<<vars_d.at(j)[0]<<" global_eright is "<<global_weight<<std::endl;
							exit(EXIT_FAILURE);
						}

					}

					//blarg, how will I treat this spectrum
					spec_CV.hist.at(j).Fill(vars_d.at(j)[0],global_weight);


				}//end of valid if-check
			}//end of CCQE..interaction type check.. OH no this should be in fillHistograms
		} //end of entry loop
	}//end of file loop 


	/***************************************************************
	 *		Now some clean-up and Writing
	 * ************************************************************/

	spec_sig.calcFullVector();
	spec_sig.writeOut("LEE_siggal.root");

	for(auto &f: files){
		f->Close();
	}

}
	/***************************************************************
	 *		Some virtual functions for selection and histogram filling
	 * ************************************************************/


bool eventSelection(int which_file){
     //from here have access to vars_i  and vars_d  to make a selection
	
     bool ans = false;
     if(vars_i.at(which_file)[0] == 1001){
	     ans = true;
     }

     return ans;
}

int fillHistograms(int file, int uni, double wei){
	//Fill the histograms
	multi_sbnspec.at(uni).hist.at(file).Fill( vars_d.at(file).at(0),wei);
return 0;
}

	/***************************************************************
	 *		And then form a covariance matrix (well 3)
	 * ************************************************************/


int SBNcovar::formCovarianceMatrix(){

	full_covariance.ResizeTo(num_bins_total, num_bins_total);
	frac_covariance.ResizeTo(num_bins_total, num_bins_total);
	full_correlation.ResizeTo(num_bins_total, num_bins_total);

	//prepare three TH2D for plotting 
	TH2D * h2 = new TH2D("Frac Cov","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);
	TH2D * h3 = new TH2D("Corr","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);
	TH2D * h4 = new TH2D("Full Cov","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);


	for(auto &h: multi_sbnspec){
		h.calcFullVector();
	}


	spec_CV.calcFullVector();	
	std::vector<double> CV = spec_CV.fullVec;

	//for(auto k=0;k<spec_CV.fullVec.size();k++){
	//	std::cout<<" "<<spec_CV.fullVec[k]<<" "<<multi_sbnspec.at(0).fullVec[k]<<std::endl;
	//}

	for(int i=0; i<num_bins_total; i++){
		for(int j=0; j<num_bins_total; j++){

			full_covariance(i,j)=0;

			for(int m=0; m < multi_sbnspec.size(); m++){
				full_covariance(i,j) += (CV[i]-multi_sbnspec.at(m).fullVec.at(i))*(CV[j]-multi_sbnspec.at(m).fullVec.at(j));


				if(full_covariance(i,j)!=full_covariance(i,j)){
					std::cout<<"ERROR: nan : at (i,j):  "<<i<<" "<<j<<" fullcov: "<<full_covariance(i,j)<<" multi hist sise "<<multi_sbnspec.size()<<" CV: "<<CV[i]<<" "<<CV[j]<<" multihisg "<<multi_sbnspec[m].fullVec[i]<<" "<<multi_sbnspec[m].fullVec[j]<<" on dim : "<<m<<std::endl;
				}



			}
			full_covariance(i,j) = full_covariance(i,j)/( multi_sbnspec.size()-1.0);



		}
	}


	for(int i=0; i<num_bins_total; i++){
		for(int j=0; j<num_bins_total; j++){

			frac_covariance(i,j) = full_covariance(i,j)/(spec_CV.fullVec[i]*spec_CV.fullVec[j]) ;
			full_correlation(i,j)= full_covariance(i,j)/(sqrt(full_covariance(i,i))*sqrt(full_covariance(j,j)));

			h2->SetBinContent(i+1,j+1,frac_covariance(i,j));
			h3->SetBinContent(i+1,j+1,full_correlation(i,j));
			h4->SetBinContent(i+1,j+1,full_covariance(i,j));

		}
	}

	/************************************************************
	 *			Saving to file				    *
	 * *********************************************************/

	TFile *ftest=new TFile("covariance_matrix.root","RECREATE");
	ftest->cd();
	TCanvas *c1 =  new TCanvas("Fractional Covariance Matrix");
	c1->cd();
	
	h2->SetTitle("Fractional Covariance Matrix (sys only)");
	h2->GetYaxis()->SetTitle("E_{#nu}^{truth}");
	h2->GetXaxis()->SetTitle("E_{#nu}^{truth}");
	h2->Draw("COLZ");
	c1->Write();

	TCanvas *c2 =  new TCanvas("Correlation Matrix");
	c2->cd();

	h3->SetTitle("Correlation Matrix (sys only)");
	h3->GetYaxis()->SetTitle("E_{#nu}^{truth}");
	h3->GetXaxis()->SetTitle("E_{#nu}^{truth}");
	h3->Draw("COLZ");
	c2->Write();	

	TCanvas *c3 =  new TCanvas("Covariance Matrix");
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


	/************************************************************
	 *		Quality Testing Suite			    *
	 * *********************************************************/


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

