#include "SBNchi.h"

void help(std::string in){

	std::cout<<in.c_str()<<std::endl;

}


/***********************************************
 *		Constructors
 * ********************************************/

SBNchi::SBNchi(SBNspec in) : SBNconfig(in.xmlname), bkgSpec(in){
		help("Starting chi constructor");	

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

		TMatrixT <double> McI(TallComp,TallComp);
		// Fill systematics from pre-computed files
		TMatrixT <double> Msys(Tall,Tall);

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
		TMatrixT <double> Mstat(Tall,Tall);
		stats_fill(Mstat, inSpec.fullVec);


		//And then define the total covariance matrix in all its glory
		TMatrixT <double > Mtotal(Tall,Tall);

		Mtotal = Mstat+Msys;
		
		// Now contract back the larger antimatrix
		TMatrixT<double > Mctotal(TallComp,TallComp);

		collapse_layer3(Mtotal, Mctotal);
			
		vMc = to_vector(Mctotal);
		// just to hold determinant
		double invdet=0; 

		// Bit o inverting, root tmatrix seems perfectly and sufficiently fast for this, even with anti_mode
		McI = Mctotal.Invert(&invdet);

		// There is currently a bug, somehow a memory leak perhaps. converting the TMatrix to a vector of vectors fixes it for now. 
		vMcI = to_vector(McI);

return 1;

}


double SBNchi::calc_chi(SBNspec sigSpec){
		double tchi = 0;	

		for(int i =0; i<TallComp; i++){
				for(int j =0; j<TallComp; j++){
						tchi += (bkgSpec.compVec[i]-sigSpec.compVec[i])*vMcI[i][j]*(bkgSpec.compVec[j]-sigSpec.compVec[j] );
				}
		}

		lastChi = tchi;
		return tchi;
}


double SBNchi::calc_chi(std::vector<double> sigVec){
		double tchi = 0;	

		if(sigVec.size() != TallComp ){
			std::cerr<<"ERROR: SBNchi::calc_chi(std::vector<double>) ~ your inputed vector does not have correct dimensions"<<std::endl;
			exit(EXIT_FAILURE);
		}
	

		for(int i =0; i<TallComp; i++){
				for(int j =0; j<TallComp; j++){
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
		return sys_fill_direct(CorrMatRoot, CorrMatName);
}


TMatrixT<double > SBNchi::sys_fill_direct(std::string rootname, std::string matname){


		TMatrixT<double> temp2(Tall,Tall);
		TFile *fm= new TFile(rootname.c_str());
		TMatrixT<float> * temp = (TMatrixT <float>* )fm->Get(matname.c_str());
	

		std::vector<std::vector<double>> mcont;

		for(int p:useBins){
			std::vector<double> tvec;
			for(int u:useBins){
				tvec.push_back( (*temp)(p,u) );
			}					
			mcont.push_back(tvec);
		}	

		


		for(int i =0; i<Tall; i++)
		{
			for(int j =0; j<Tall; j++)
			{
				temp2(i,j)=mcont[i][j];
			}
		}
		delete temp;


		fm->Close();
		delete fm;
		return temp2;


}


//This is the powerhouse, takes each detector matrix filled with Nchan channels of Chan[i] subchannels, and collapses it.
void SBNchi::collapse_layer1(TMatrixT <double> & M, TMatrixT <double> & Mc){
	bool debug = false;
	if(debug)	std::cout<<"Starting:M "<<M.GetNcols()<<" "<<M.GetNrows()<<" "<<115<<std::endl;
	if(debug)	std::cout<<"Starting:Mc "<<Mc.GetNcols()<<" "<<Mc.GetNrows()<<" "<<30<<std::endl;


		
		std::vector<std::vector<TMatrixT<double>>> Summed(Nchan, std::vector<TMatrixT<double>>(Nchan) );	//Initialise a matrix of matricies, to ZERO.
		for(int ic = 0; ic < Nchan; ic++){ 
			for(int jc =0; jc < Nchan; jc++){
			Summed[ic][jc].ResizeTo(Bins[jc],Bins[ic]) ;// This is CORRECT, do not switch (ie Summed[0][1] = size (Bins[1], Bins[0])
			Summed[ic][jc] = 0.0;
			}
		}
		

		int mrow = 0.0;
		int mcol = 0.0;

		for(int ic = 0; ic < Nchan; ic++){ 	  //Loop over all rows
			for(int jc =0; jc < Nchan; jc++){ //Loop over all columns
				
							if(debug)std::cout<<"Diagonal! : "<<ic<<" "<<jc<<" mcol is: "<<mcol<<" mrow is: "<<mrow<<std::endl;				
					
				for(int m=0; m < Chan[ic]; m++){
					for(int n=0; n< Chan[jc]; n++){ //For each big block, loop over all subchannels summing together
						Summed[ic][jc] +=  M.GetSub(mrow+n*Bins[jc] ,mrow + n*Bins[jc]+Bins[jc]-1, mcol + m*Bins[ic], mcol+ m*Bins[ic]+Bins[ic]-1 );
					}
				}	
				

				mrow += Chan[jc]*Bins[jc];//As we work our way left in columns, add on that many bins
			}//end of column loop

			mrow = 0; // as we end this row, reset row count, but jump down 1 column
			mcol += Chan[ic]*Bins[ic];
		}//end of row loop

	

///********************************* And put them back together! ************************//
	Mc.Zero(); 
	mrow = 0;
	mcol = 0;

	//Repeat again for Contracted matrix
	for(int ic = 0; ic < Nchan; ic++){ 	  
		for(int jc =0; jc < Nchan; jc++){ 
	
				Mc.SetSub(mrow,mcol,Summed[ic][jc]);	
				mrow += Bins[jc];
			}

		mrow = 0;
		mcol +=Bins[ic];
	}

return;
}



//This is the detector layer, Take a given mode and run over each detector V detector sub matrix
void SBNchi::collapse_layer2(TMatrixT <double> & M, TMatrixT <double> & Mc){
		Mc.Zero();
		int nrow = Tdet;// N_e_bins*N_e_spectra+N_m_bins*N_m_spectra;
		int crow = TdetComp; //N_e_bins+N_m_bins;
	
		for(int m =0; m< Ndet; m++){
			for(int n =0; n< Ndet; n++){
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
		int nrow = Tmode;// (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
		int crow=  TmodeComp;// (N_e_bins+N_m_bins)*N_dets;
	
		for(int m =0; m< Nmode ; m++){
			for(int n =0; n< Nmode; n++){
				
				TMatrixT<double> imat(nrow,nrow);
				TMatrixT<double> imatc(crow,crow);
				
				imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);
				
				collapse_layer2(imat,imatc);
				Mc.SetSub(n*crow,m*crow,imatc);

			}
		}

return;
}


