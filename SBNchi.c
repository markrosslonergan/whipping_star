#include "SBNchi.h"
using namespace sbn;


void help(std::string in){

	std::cout<<in.c_str()<<std::endl;

}


/***********************************************
 *		Constructors
 * ********************************************/

SBNchi::SBNchi(SBNspec in) : SBNconfig(in.xmlname), bkgSpec(in){
		lastChi = -9999999;
		load_bkg(bkgSpec);
}


SBNchi::SBNchi(SBNspec in, std::string newxmlname) : SBNconfig(newxmlname), bkgSpec(in){
		
		if(fullnames.size() !=in.fullnames.size()){
			std::cerr<<"ERROR: SBNchi::SBNchi | Selected covariance matrix and background spectrum are different sizes!"<<std::endl;
			exit(EXIT_FAILURE);
		}else{		
		 	for(int i=0; i< fullnames.size(); i++){
				if(fullnames[i]!=in.fullnames[i]){
					std::cerr<<"ERROR: SBNchi::SBNchi | Spectrum and Covariance matrix have different (or different order) subchannels!"<<std::endl;
					exit(EXIT_FAILURE);
				}
			}
		}		


		lastChi = -9999999;

		load_bkg(bkgSpec);
}


/***********************************************
 *		Rest for now
 * ********************************************/

int SBNchi::load_bkg(SBNspec inSpec){

		TMatrixT <double> McI(num_bins_total_compressed,num_bins_total_compressed);
		// Fill systematics from pre-computed files
		TMatrixT <double> Msys(num_bins_total,num_bins_total);

		Msys = sys_fill_direct();



		// systematics per scaled event
		for(int i =0; i<Msys.GetNcols(); i++)
		{
			for(int j =0; j<Msys.GetNrows(); j++)
			{
				Msys(i,j)=Msys(i,j)*inSpec.fullVec[i]*inSpec.fullVec[j];
			}
		}
			
		// Fill stats from the back ground vector
		TMatrixT <double> Mstat(num_bins_total,num_bins_total);
		stats_fill(Mstat, inSpec.fullVec);


		//And then define the total covariance matrix in all its glory
		TMatrixT <double > Mtotal(num_bins_total,num_bins_total);

		Mtotal = Mstat+Msys;
		
		// Now contract back the larger antimatrix
		TMatrixT<double > Mctotal(num_bins_total_compressed,num_bins_total_compressed);

		collapse_layer3(Mtotal, Mctotal);
			
		vMc = to_vector(Mctotal);
		// just to hold determinant
		double invdet=0; 

		// Bit o inverting, root tmatrix seems perfectly and sufficiently fast for this, even with anti_mode
		McI = Mctotal.Invert(&invdet);

		// There is currently a bug, somehow a memory leak perhaps. converting the TMatrix to a vector of vectors fixes it for now. 
		vMcI = to_vector(McI);

return 0;

}


double SBNchi::CalcChi(SBNspec sigSpec){
		double tchi = 0;	

		if(sigSpec.compVec.size()==0){
			std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
			sigSpec.compressVector();
		}

		for(int i =0; i<num_bins_total_compressed; i++){
				for(int j =0; j<num_bins_total_compressed; j++){
						tchi += (bkgSpec.compVec[i]-sigSpec.compVec[i])*vMcI[i][j]*(bkgSpec.compVec[j]-sigSpec.compVec[j] );
				}
		}

		lastChi = tchi;
		return tchi;
}


double SBNchi::CalcChi(std::vector<double> sigVec){
		double tchi = 0;	

		if(sigVec.size() != num_bins_total_compressed ){
			std::cerr<<"ERROR: SBNchi::CalcChi(std::vector<double>) ~ your inputed vector does not have correct dimensions"<<std::endl;
			std::cerr<<"sigVec.size(): "<<sigVec.size()<<" num_bins_total_compressed"<<num_bins_total_compressed<<std::endl;
			exit(EXIT_FAILURE);
		}
	

		for(int i =0; i<num_bins_total_compressed; i++){
				for(int j =0; j<num_bins_total_compressed; j++){
						tchi += (bkgSpec.compVec[i]-sigVec[i])*vMcI[i][j]*(bkgSpec.compVec[j]-sigVec[j] );
				}
		}

		lastChi = tchi;
		return tchi;
}




void SBNchi::fake_fill(TMatrixT <double> &M){

	//Fills a square matrix of dim matrix_size with random numbers for now.
	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();
	
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,std::time(0));

	int matrix_size=M.GetNrows();

	if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}
   	 
	for(int i=0; i<matrix_size; i++){
		for (int j = i;j<matrix_size;j++){
				M(i,j)=gsl_rng_uniform(r);
				M(j,i)=M(i,j);
		}
	
	}

	gsl_rng_free(r);
 return ;
}


std::vector<std::vector<double >> SBNchi::to_vector(TMatrixT <double > Min)
{
	int dimension =  Min.GetNrows();

	std::vector<std::vector<double >>  ans(dimension, std::vector<double>(dimension));

	for(int i = 0; i< dimension; i++){
		for(int k = 0; k< dimension; k++){
			ans[i][k]=Min(i,k);
			if(ans[i][k]==-0){
				ans[i][k]=0;
			}
		}	
	}
return ans;


}


void SBNchi::stats_fill(TMatrixT <double> &M, std::vector<double> diag){
	int matrix_size = M.GetNrows();
	

	if(matrix_size != diag.size()){std::cout<<"#ERROR: stats_fill, matrix not equal to diagonal"<<std::endl;}
	if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}


	M.Zero();
	
	for(int i=0; i<matrix_size; i++)
	{
		M(i,i) = diag[i];	

	}



 return ;
}



TMatrixT<double> SBNchi::sys_fill_direct(){
		return sys_fill_direct(correlation_matrix_rootfile, correlation_matrix_name);
}


TMatrixT<double > SBNchi::sys_fill_direct(std::string rootname, std::string matname){


		TMatrixT<double> temp2(num_bins_total,num_bins_total);
		TFile *fm= new TFile(rootname.c_str());
		TMatrixT<float> * temp = (TMatrixT <float>* )fm->Get(matname.c_str());
	

		std::vector<std::vector<double>> mcont;

		for(int p:used_bins){
			std::vector<double> tvec;
			for(int u:used_bins){
				tvec.push_back( (*temp)(p,u) );
			}					
			mcont.push_back(tvec);
		}	

		


		for(int i =0; i<num_bins_total; i++)
		{
			for(int j =0; j<num_bins_total; j++)
			{
				temp2(i,j)=mcont[i][j];
			}
		}
		delete temp;


		fm->Close();
		delete fm;
		return temp2;


}


//This is the powerhouse, takes each detector matrix filled with num_channels channels of num_subchannels[i] subchannels, and collapses it.
void SBNchi::collapse_layer1(TMatrixT <double> & M, TMatrixT <double> & Mc){
	bool debug = false;
	if(debug)	std::cout<<"Starting:M "<<M.GetNcols()<<" "<<M.GetNrows()<<" "<<115<<std::endl;
	if(debug)	std::cout<<"Starting:Mc "<<Mc.GetNcols()<<" "<<Mc.GetNrows()<<" "<<30<<std::endl;


		
		std::vector<std::vector<TMatrixT<double>>> Summed(num_channels, std::vector<TMatrixT<double>>(num_channels) );	//Initialise a matrix of matricies, to ZERO.
		for(int ic = 0; ic < num_channels; ic++){ 
			for(int jc =0; jc < num_channels; jc++){
			Summed[ic][jc].ResizeTo(num_bins[jc],num_bins[ic]) ;// This is CORRECT, do not switch (ie Summed[0][1] = size (num_bins[1], num_bins[0])
			Summed[ic][jc] = 0.0;
			}
		}
		

		int mrow = 0.0;
		int mcol = 0.0;

		for(int ic = 0; ic < num_channels; ic++){ 	  //Loop over all rows
			for(int jc =0; jc < num_channels; jc++){ //Loop over all columns
				
							if(debug)std::cout<<"Diagonal! : "<<ic<<" "<<jc<<" mcol is: "<<mcol<<" mrow is: "<<mrow<<std::endl;				
					
				for(int m=0; m < num_subchannels[ic]; m++){
					for(int n=0; n< num_subchannels[jc]; n++){ //For each big block, loop over all subchannels summing together
						Summed[ic][jc] +=  M.GetSub(mrow+n*num_bins[jc] ,mrow + n*num_bins[jc]+num_bins[jc]-1, mcol + m*num_bins[ic], mcol+ m*num_bins[ic]+num_bins[ic]-1 );
					}
				}	
				

				mrow += num_subchannels[jc]*num_bins[jc];//As we work our way left in columns, add on that many bins
			}//end of column loop

			mrow = 0; // as we end this row, reset row count, but jump down 1 column
			mcol += num_subchannels[ic]*num_bins[ic];
		}//end of row loop

	

///********************************* And put them back together! ************************//
	Mc.Zero(); 
	mrow = 0;
	mcol = 0;

	//Repeat again for Contracted matrix
	for(int ic = 0; ic < num_channels; ic++){ 	  
		for(int jc =0; jc < num_channels; jc++){ 
	
				Mc.SetSub(mrow,mcol,Summed[ic][jc]);	
				mrow += num_bins[jc];
			}

		mrow = 0;
		mcol +=num_bins[ic];
	}

return;
}



//This is the detector layer, Take a given mode and run over each detector V detector sub matrix
void SBNchi::collapse_layer2(TMatrixT <double> & M, TMatrixT <double> & Mc){
		Mc.Zero();
		int nrow = num_bins_detector_block;// N_e_bins*N_e_spectra+N_m_bins*N_m_spectra;
		int crow = num_bins_detector_block_compressed; //N_e_bins+N_m_bins;
	
		for(int m =0; m< num_detectors; m++){
			for(int n =0; n< num_detectors; n++){
				TMatrixT<double> imat(nrow,nrow);
				TMatrixT<double> imatc(crow,crow);
				
				imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);
				collapse_layer1(imat,imatc);
				Mc.SetSub(n*crow,m*crow,imatc);
			}
		}

return;
}
	
//This is the Mode layer, Take a given full matrix and runs over each Mode V Mode sub matrix
void SBNchi::collapse_layer3(TMatrixT <double> & M, TMatrixT <double> & Mc){
		Mc.Zero();
		int nrow = num_bins_mode_block;// (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
		int crow=  num_bins_mode_block_compressed;// (N_e_bins+N_m_bins)*N_dets;
	
		for(int m =0; m< num_modes ; m++){
			for(int n =0; n< num_modes; n++){
				
				TMatrixT<double> imat(nrow,nrow);
				TMatrixT<double> imatc(crow,crow);
				
				imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);
				
				collapse_layer2(imat,imatc);
				Mc.SetSub(n*crow,m*crow,imatc);

			}
		}

return;
}


