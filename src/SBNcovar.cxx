#include "SBNcovar.h"
//#include "MCEventWeight.h"

using namespace sbn;

SBNcovar::SBNcovar(std::string xmlname) : SBNconfig(xmlname) {
	gROOT->ProcessLine("#include <map>");
	gROOT->ProcessLine("#include <vector>");
	gROOT->ProcessLine("#include <string>");
	gSystem->Load("/uboone/app/users/markrl/sbnfit/whipping_star/src/mdict_h.so");
	gStyle->SetOptStat(0);

	universes_used = 0;
	tolerence_positivesemi = 1e-10;
	is_small_negative_eigenvalue = false;
	bool DEBUG_BOOL = false;


	//std::cout<<"parm size: "<<parameter_names.size()<<std::endl;
	///for(int i=0; i< parameter_names.size(); i++){
	//	std::cout<<i<<" 1down"<<parameter_names.at(i).size()<<std::endl;
	//	for(auto &nam: parameter_names.at(i)){
	///		std::cout<<nam<<std::endl;
	//	}
	//	}

	//Initialise all the things
	//for every multisim, create a SBNspec
	for(int m=0; m<num_multisim[0]; m++){
		SBNspec temp_spec(xmlname,m);
		multi_sbnspec.push_back(temp_spec);
	}

	
	

	//std::cout<<" multi_hist_size: "<<multi_sbnspec.size()<<" internal size: "<<multi_sbnspec.at(0).hist.size()<<"  deeper: "<<multi_sbnspec.at(0).hist.at(0).GetBinContent(1)<<std::endl;

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
		//nentries.push_back(1000);
		nentries.push_back(t->GetEntries());
	}


	//signal crap todo
	//
	TFile *fS = new TFile("/uboone/app/users/mastbaum/leerw/leerw.root");
	TTree * tnS =  (TTree*)fS->Get("leerw;1");
	float lee=0;
	tnS->SetBranchAddress("leerw",&lee);
	SBNspec spec_sig(xmlname,-1);


	vars_i= std::vector<std::vector<int>>(Nfiles   , std::vector<int>(branch_names_int.at(0).size(),0));
	vars_d=std::vector<std::vector<double>>(Nfiles   , std::vector<double>(branch_names_double.at(0).size(),0.0));

	std::vector< std::map<std::string, std::vector<double> > * > fWeights(multisim_name.size(),0);

	std::vector< TBranch *> bweight(Nfiles,0);

	for(int i=0; i< Nfiles; i++){
		trees.at(i)->SetBranchAddress("mcweight", &(fWeights.at(i)), &(bweight.at(i)) );

		for(auto &bfni: branch_names_int){
			for(int k=0; k< bfni.size();k++){
				trees.at(i)->SetBranchAddress(bfni[k].c_str(), &(vars_i.at(i).at(k)));
			}
		}
		for(auto &bfnd: branch_names_double){
			for(int k=0; k< bfnd.size();k++){
				trees.at(i)->SetBranchAddress(bfnd[k].c_str(), &(vars_d.at(i).at(k)));
			}
		}
	}




	for(int j=0;j<Nfiles;j++){
		double pot_factor = pot.at(j)/(pot_scaling.at(j) * (double)nentries.at(j));

		for(int i=0; i< nentries.at(j); i++){
			bool is_valid = true;
			int n_uni = 0;
			std::vector<double> weights;
			double global_weight = 1;

			if(i%2500==0)std::cout<<"Event: "<<i<<" of "<<nentries[j]<<" from File: "<<multisim_file[j]<<" POT factor: "<<pot_factor<<std::endl;

			trees.at(j)->GetEntry(i);
			if(j==0){tnS->GetEntry(i);}
			//here we put the selection criteria, for example nuance interaction 1001 == CCQE, virtual bool
			if( this->eventSelection(j) ){

				//Loop over all things in mcweight map
				for(std::map<std::string, std::vector<double> >::iterator  it = fWeights.at(j)->begin(); it != fWeights.at(j)->end(); ++it) {
					bool isThis=false;
					
					for(auto nam: parameter_names.at(j)){
						if(nam == it->first || nam == "ALL"){
							isThis=true;				
							break;	
						}
					}

					if(!isThis) continue;
					if(!is_valid) break;

					if(it->first == "bnbcorrection_FluxHist"){
						global_weight = it->second.at(0);		

						//Skip events that fail (should they be included in POT?)
						if(std::isinf(global_weight)){
							is_valid=false;
							break;
						}
						global_weight = global_weight*pot_factor;
						
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

						weights.push_back(wei*global_weight);
						n_uni++;

					}

				}//end of weight (it->second) iterator


				if(is_valid){
					// loop over every multisim we have in this entry
					double num_sim = weights.size();	 

					if(num_sim > multi_sbnspec.size()){
						//std::cout<<"ERROR: numu loop numsim > nulti_hist"<<num_sim<<" "<<multi_sbnspec.size()<<std::endl;
						num_sim = multi_sbnspec.size();
					}

					for(int m=0; m<num_sim; m++){
						//This is the part where we will every histogram in this Universe
						this->fillHistograms(j, m, weights.at(m) );

						universes_used++;
						//important check. failure mode	
						if(weights[m]!=weights[m] || isinf(weights[m]) ){
							std::cout<<"ERROR: weight has a value of: "<<weights.at(m)<<". So I am killing all. on Dim: "<<m<<" energy"<<vars_d.at(j)[0]<<" global_eright is "<<global_weight<<std::endl;
							exit(EXIT_FAILURE);
						}

					}

					//blarg, how will I treat this spectrum
					spec_CV.hist.at(j).Fill(vars_d.at(j)[0],global_weight);
					
					if(j==0){
						spec_sig.hist.at(j).Fill(vars_d.at(j)[0],global_weight*(1+lee));
					}


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


bool SBNcovar::eventSelection(int which_file){
	//from here have access to vars_i  and vars_d  to make a selection

	bool ans = false;
	if(vars_i.at(which_file)[0] == 1001){
		ans = true;
	}

	return ans;
}

int SBNcovar::fillHistograms(int file, int uni, double wei){
	//Fill the histograms
	//
/*	TRandom3 rangen(0);
	double sigma = 0;
	if(file==0){
		sigma=;
	}else if(file==1){
		sigma=;
	}	
	double en = rangen.Gaus( vars_d.at(file), sigma );
*/	
	double en = vars_d.at(file)[0];

	multi_sbnspec.at(uni).hist.at(file).Fill(en, wei);
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
			full_covariance(i,j) = full_covariance(i,j)/( universes_used-1.0);



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
	std::string nn = "covariance_matrix_" + parameter_names[0][0]+".root";

	TFile *ftest=new TFile(nn.c_str(),"RECREATE");
	ftest->cd();


	TCanvas *cspline =  new TCanvas("Splines");
	int num_hists = multi_sbnspec.at(0).hist.size();
	cspline->Divide(num_hists,1);

	for(int h=0; h<spec_CV.hist.size(); h++){
		cspline->cd(h+1);
		spec_CV.hist.at(h).SetLineColor(kBlack);
		spec_CV.hist.at(h).SetLineWidth(4);
		spec_CV.hist.at(h).SetTitle( (fullnames[h]+" : "+parameter_names[0][0]).c_str() );
		spec_CV.hist.at(h).Draw("L SAME");
	}

	for(int m=0; m< multi_sbnspec.size(); m++){
		for(int h=0; h<multi_sbnspec.at(m).hist.size(); h++){
			cspline->cd(h+1);
			TRandom3 * rangen = new TRandom3(0);
			multi_sbnspec.at(m).hist.at(h).SetLineColor(rangen->Uniform(400,900));
	
			multi_sbnspec.at(m).hist.at(h).Draw("L SAME");
		}
	}

	for(int h=0; h<spec_CV.hist.size(); h++){
		cspline->cd(h+1);
		spec_CV.hist.at(h).SetLineWidth(4);
		spec_CV.hist.at(h).Draw("L SAME");
	}


	cspline->Write();
	std::string ppsp = "splines_"+parameter_names[0][0]+".pdf";

	cspline->SaveAs(ppsp.c_str());







	//matricies
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

	std::string pp = "Fractional Covarariance and Correlation: "+parameter_names[0][0];
	TCanvas *cboth = new TCanvas(pp.c_str());
	cboth->SetCanvasSize(1200,600);

	h2->SetTitle( ("Fractional Covariance: "+ parameter_names[0][0]).c_str());
	h3->SetTitle( ("Correlation: "+ parameter_names[0][0]).c_str());
	cboth->Divide(2,1);
	cboth->SetFixedAspectRatio();
	cboth->cd(1);

	h2->Draw("COLZ");
	cboth->cd(2);
	h3->Draw("COLZ");
	cboth->Write();
	std::string ppdf = "covar_plots_"+parameter_names[0][0]+".pdf";
	cboth->SaveAs(ppdf.c_str());


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

