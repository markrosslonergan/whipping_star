#include "SBNcovar.h"
using namespace sbn;

SBNcovar::SBNcovar(std::string rootfile, std::string xmlname) : SBNconfig(xmlname) {

	tolerence_positivesemi = 1e-10;
	is_small_negative_eigenvalue = false;


	//Initialise all the things
	//for every multisim, create a SBNspec
	for(int m=0; m<num_multisim; m++){
		SBNspec temp_spec(xmlname,m);
		multi_hists.push_back(temp_spec);
	}


	char namei[200];
	sprintf(namei,"%s.root",rootfile.c_str());	
	TFile *f = new TFile(namei);

	TTree * full_mc=  (TTree*)f->Get(multisim_name.c_str());

	std::vector<double> *weights=0 ;
	std::vector<double> vars (branch_names.size(),0); 


	TBranch *bweight=0;
	full_mc->SetBranchAddress("weights", &weights, &bweight);		

	for(int i=0; i< branch_names.size(); i++){
		if(branch_names[i]!="weights"){

			full_mc->SetBranchAddress(branch_names[i].c_str(), &(vars[i]) );
		}
	}



	 int nentries = full_mc->GetEntries(); // read the number of entries


	 // So for every entry..
	 for (int i = 0; i < nentries; i++) {
		 full_mc->GetEntry(i);

		 // loop over every multisim we have in this entry
		 // the first must be 1, i.e the CV or mean case.
		 //
		 for(int m=0; m<num_multisim; m++){
			
			if(m==0 && weights->at(m)!=1){
				std::cerr<<"ERROR: Im just going to arbitrarily say the first multisim MUST have weight 1. ie the CV"<<std::endl;
				exit(EXIT_FAILURE);
			}
			
			//Then in this sim, fill every histogram! 
			for(auto & h: multi_hists[m].hist){	
				h.Fill(vars[0],weights->at(m));
			}
		 }
	 }



	f->Close();

}



int SBNcovar::formCovarianceMatrix(){

	full_covariance.ResizeTo(num_bins_total, num_bins_total );
	TH2D * h2 = new TH2D("test","",num_bins_total,1,num_bins_total, num_bins_total,1,num_bins_total);

	for(int i=0; i< num_multisim;i++){
		multi_hists[i].calcFullVector();
	}
	
	std::vector<double> CV = multi_hists[0].fullVec;

	for(int i=0; i<num_bins_total; i++){
	  for(int j=0; j<num_bins_total; j++){
		
		full_covariance(i,j)=0;

		//Remember that m=0 is CV!
		for(int m=1; m < num_multisim; m++){
			full_covariance(i,j) += (CV[i]-multi_hists[m].fullVec[i])*(CV[j]-multi_hists[m].fullVec[j]);
		}
		full_covariance(i,j) = full_covariance(i,j)/( (double)num_multisim-1.0);
		h2->SetBinContent(i+1,j+1,full_covariance(i,j));

	  }
	}

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



	TFile *ftest=new TFile("qtest.root","RECREATE");
	ftest->cd();
	TCanvas *c1 =  new TCanvas();
	c1->cd();
	//full_covariance.Draw();
	h2->Draw("COLZ");
	c1->Write();
	ftest->Close();






return 0;
}

