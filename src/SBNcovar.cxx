#include "SBNcovar.h"
//#include "MCEventWeight.h"

using namespace sbn;

SBNcovar::SBNcovar(std::string rootfile, std::string xmlname) : SBNconfig(xmlname) {
	std::cout<<"inside SBNcovar"<<std::endl;	
	gROOT->ProcessLine("#include <map>");
	gROOT->ProcessLine("#include <vector>");
	gROOT->ProcessLine("#include <string>");
	gSystem->Load("/home/mark/work/sbnfit_reduce/src/mdict_h.so");


	tolerence_positivesemi = 1e-10;
	is_small_negative_eigenvalue = false;


	//Initialise all the things
	//for every multisim, create a SBNspec
	for(int m=0; m<num_multisim; m++){
		SBNspec temp_spec(xmlname,m);
		multi_hists.push_back(temp_spec);
	}


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


		 // So for every entry..
		 std::cout<<"Starting loop over entries: "<<nentries<<std::endl;
		 for (int i = 0; i < nentries; i++) {
			weights.clear(); 
			 double glob_corr=1;
			 
		 	 std::cout<<"On entry :"<<i<<" over entries: "<<nentries<<std::endl;
			 full_mc->GetEntry(i);
		 	 std::cout<<"Gotten entry :"<<i<<" over entries: "<<nentries<<std::endl;
			
			 for(int h=0;h<branch_names_double.size();h++){
				std::cout<<h<<" "<<branch_names_double[h].c_str()<<" "<<vars_d[h]<<std::endl;
			 }
	
			 for(int h=0;h<branch_names_int.size();h++){
				std::cout<<h<<" "<<branch_names_int[h].c_str()<<" "<<vars_i[h]<<std::endl;
			 }

			
			std::cout<<"Gonna loop over weight maps now:"<<std::endl;			
			for(std::map<std::string, std::vector<double> >::iterator  it = fWeight->begin(); it != fWeight->end(); ++it) {
				std::cout << it->first <<" "<<it->second.size()<<std::endl;
				
				if(it->first == "bnbcorrection_FluxHist"){
					std::cout<<"Got the bnb flux correction"<<std::endl;	
					glob_corr = it->second.at(0);		
					continue;
				}

				std::cout<<" starting weight loop:"<<std::endl;
				for(auto &j: it->second){
					//	std::cout<<"Weight: "<<j<<std::endl;
						weights.push_back(j*glob_corr);
				}
				
			}


			 // loop over every multisim we have in this entry
			 // the first must be 1, i.e the CV or mean case.
			double num_sim = weights.size();	 

			std::cout<<"Num Sim "<<num_sim<<std::endl;
			for(int m=0; m<num_sim; m++){
				//std::cout<<"On Uni :"<<m<<" of "<<num_sim<<std::endl;
				//Then in this sim, fill every histogram! 
				for(auto & h: multi_hists[m].hist){	
					h.Fill(vars_i[0],weights[m]);
				}
			 }
		 } //end of entry loop

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

