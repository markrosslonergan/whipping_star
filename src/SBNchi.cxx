#include "SBNchi.h"
using namespace sbn;


/***********************************************
 *		Constructors
 * ********************************************/

SBNchi::SBNchi(SBNspec in) : SBNchi(in,false){}

SBNchi::SBNchi(SBNspec in, std::string newxmlname) : SBNconfig(newxmlname), bkgSpec(in){
	stat_only = false;

	matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
	MfracCov.ResizeTo(num_bins_total, num_bins_total);
	
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

	MfracCov = sys_fill_direct();
	lastChi = -9999999;
	bkgSpec.compressVector();
	Msys.ResizeTo(num_bins_total,num_bins_total);
	Msys.Zero();
	Msys=MfracCov;

	this->reload_core_spec(&in);
}


SBNchi::SBNchi(SBNspec in, TMatrixT<double> Msysin) : SBNconfig(in.xmlname), bkgSpec(in){
	lastChi = -9999999;
	stat_only= false;
	
	matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
	Msys.ResizeTo(Msysin.GetNrows(), Msysin.GetNcols());
	MfracCov.ResizeTo(Msysin.GetNrows(), Msysin.GetNcols());
	
	MfracCov = Msysin;
	Msys.Zero();

	this->reload_core_spec(&bkgSpec);

}

SBNchi::SBNchi(SBNspec in, bool is_stat_only): SBNconfig(in.xmlname), bkgSpec(in), stat_only(is_stat_only){
	lastChi = -9999999;
	
	matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
	Msys.ResizeTo(num_bins_total, num_bins_total);
	MfracCov.ResizeTo(num_bins_total, num_bins_total);

	
	if(is_stat_only){
		MfracCov.Zero();	
		Msys.Zero();

	}else{
		MfracCov = sys_fill_direct();
		Msys.Zero();
	}
	

	this->reload_core_spec(&bkgSpec);

}



/***********************************************
 *		Rest for now
 * ********************************************/

int SBNchi::reload_core_spec(SBNspec *bkgin){
	if(isVerbose)std::cout<<"SBNchi::reload_core_spec || Begininning to reload core spec! First set new core spec"<<std::endl;
	bkgSpec = *bkgin;
	bkgSpec.compressVector();
	
	if(isVerbose)std::cout<<"SBNchi::reload_core_spec || Clear all previous chi^2 data"<<std::endl;
	lastChi_vec.clear();
	lastChi_vec.resize(num_bins_total_compressed, std::vector<double>( num_bins_total_compressed,0) );
	
	//Reset Msys to fractional
	if(isVerbose)std::cout<<"SBNchi::reload_core_spec, Reseting Msys to MfracCov"<<std::endl;
	Msys = MfracCov;


	if(Msys.IsSymmetric()){
		if(isVerbose)std::cout<<"At Start of reload_core_spec Inputted Msys covariance matrix is symmetric"<<std::endl;
	}else{
		std::cerr<<"ERROR: SBNchi::reload_core_spec, Msys input is not symmetric!"<<std::endl;
		//exit(EXIT_FAILURE);
	}


	if(Msys.GetNcols()!=num_bins_total ){
		std::cerr<<"ERROR: trying to pass a matrix to SBNchi that isnt the right size"<<std::endl;
		std::cerr<<"ERROR: num_bins_total: "<<num_bins_total<<" and matrix is: "<<Msys.GetNcols()<<std::endl;
		exit(EXIT_FAILURE);
	}


	if(isVerbose)std::cout<<"SBNchi::reload_core_spec || Go from fracCovariance to fullCovariance. Msys.GetNcols(): "<<Msys.GetNcols()<<" Msys.GetNrows(): "<<Msys.GetNrows()<<" core->fullvec.size(): "<<bkgSpec.fullVec.size()<<std::endl;

	// systematics per scaled event
	for(int i =0; i<Msys.GetNcols(); i++)
	{
		for(int j =0; j<Msys.GetNrows(); j++)
		{

			Msys(i,j)=Msys(i,j)*bkgSpec.fullVec.at(i)*bkgSpec.fullVec.at(j);

		}
	}

	if(isVerbose)std::cout<<"SBNchi::reload_core_spec || Filling stats into cov matrix"<<std::endl;
	// Fill stats from the back ground vector
	TMatrixT <double> Mstat(num_bins_total, num_bins_total);
	stats_fill(Mstat, bkgSpec.fullVec);

	

	if(Mstat.IsSymmetric()){
		if(isVerbose)std::cout<<"Stat matrix is symmetric (it is just diagonal core)"<<std::endl;
	}else{
		std::cerr<<"ERROR: SBNchi::formCovarianceMatrix, stats  is not symmetric!"<<std::endl;
		exit(EXIT_FAILURE);
	}


	//And then define the total covariance matrix in all its glory
	TMatrixT <double> Mtotal(num_bins_total,num_bins_total);
	Mtotal.Zero();

	if(stat_only){
		if(isVerbose)std::cout<<"SBNchi::SBNchi(SBNspec,TMatrixD) || Using stats only in covariance matrix"<<std::endl;
		Mtotal = Mstat;
	}else{
		if(isVerbose)std::cout<<"SBNchi::SBNchi(SBNspec, TMatrixD) || Using stats+sys in covariance matrix"<<std::endl;
		Mtotal = Mstat + Msys;
	}

	if(isVerbose)std::cout<<"Mstat: "<<Mstat.GetNrows()<<" x "<<Mstat.GetNcols()<<std::endl;
	if(isVerbose)std::cout<<"Msys: "<<Msys.GetNrows()<<" x "<<Msys.GetNcols()<<std::endl;
	if(isVerbose)std::cout<<"Mtotal: "<<Mtotal.GetNrows()<<" x "<<Mtotal.GetNcols()<<std::endl;

	if(Mtotal.IsSymmetric() ){
	if(isVerbose)	std::cout<<"Total Mstat +Msys is symmetric"<<std::endl;
	}else{
		
		double tol = 1e-13;
		double biggest_deviation = 0;
		int bi =0;
		int bj=0;

		if(isVerbose)std::cerr<<"WARNING: SBNchi::SBNchi(SBNspec, TMatrixD), stats + sys result appears to be not symmetric!"<<std::endl;
		for(int i=0; i<Mtotal.GetNrows(); i++){
			for(int j=0; j<Mtotal.GetNcols(); j++){
				double dev = fabs(Mtotal(i,j)-Mtotal(j,i));
				if(dev>biggest_deviation){
					biggest_deviation = 2*dev/(fabs(Mtotal(i,j))+fabs(Mtotal(j,i)));
					bi=i;
					bj=j;
				}	
			}
		}	

	if(isVerbose)	std::cerr<<"WARNING: SBNchi::SBNchi(SBNspec, TMatrixD) Biggest Relative Deviation from symmetry is i:"<<bi<<" j: "<<bj<<" of order "<<biggest_deviation<<" M(j,i)"<<Mtotal(bj,bi)<<" M(i,j)"<<Mtotal(bi,bj)<<std::endl;
	
		if(biggest_deviation >tol){

			std::cerr<<"ERROR: SBNchi::SBNchi(SBNspec, TMatrixD) Thats too unsymettric, killing."<<std::endl;

			exit(EXIT_FAILURE);
		}else{

		if(isVerbose)	std::cerr<<"WARNING: SBNchi::SBNchi(SBNspec, TMatrixD) Thats within tolderence. Continuing."<<std::endl;
		}
	}

	TMatrixT<double > Mctotal(num_bins_total_compressed,num_bins_total_compressed);
	collapse_layer3(Mtotal, Mctotal);

	matrix_collapsed = Mctotal;

	for(int j = 0; j< num_bins_total_compressed; j++){
		Mctotal(j,j) = Mctotal(j,j)+1e-14;
	}

	vMc = to_vector(Mctotal);
	double invdet=0;

	TMatrixD McI(num_bins_total_compressed,num_bins_total_compressed);
	McI.Zero();	
	
	std::cout<<"About to do a SVD decomposition"<<std::endl;	
	TDecompSVD svd(Mctotal);
  	 if (!svd.Decompose()) {
  	    std::cout << "Decomposition failed, matrix not symettric? More likely singular, maybe try to rebin" << std::endl;
	    std::cout<<"ERROR: SBNchi::SBNchi(SBNspec, TMatrixD) The matrix to invert failed a SVD decomp!"<<std::endl;
	    std::cout<<"CONDITION: "<<svd.Condition()<<std::endl;
	    for(int j = 0; j< num_bins_total_compressed; j++){
		std::cout<<j<<"/"<<num_bins_total_compressed<<" "<<Mctotal(j,j)<<" "<<bkgSpec.compVec.at(j)<<std::endl;
	    }
	    exit(EXIT_FAILURE);

  	 } else {
  	    McI = svd.Invert();
   	}
	std::cout<<"Inverted!"<<std::endl;
	vMcI = to_vector(McI);

	
	//std::cout<<"Going to invert the stats+sys matrix now"<<std::endl;
	//McI = Mctotal.Invert(&invdet);
	
	if( !McI.IsValid()){
			std::cerr<<"ERROR: SBNchi::SBNchi(SBNspec, TMatrixD) The inverted matrix isnt valid"<<std::endl;
			exit(EXIT_FAILURE);

	}
	


	// test for validity
	bool is_small_negative_eigenvalue = false;
	double tolerence_positivesemi = 1e-5; 


	//if a matrix is (a) real and (b) symmetric (checked above) then to prove positive semi-definite, we just need to check eigenvalues and >=0;
	/*TMatrixDEigen eigen (Mtotal);
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
	*/
	if(is_small_negative_eigenvalue){	
	if(isVerbose)	std::cout<<"Generated covariance matrix is (allmost) positive semi-definite. It did contain small negative values of absolute value <= :"<<tolerence_positivesemi<<std::endl;
	}else{
	if(isVerbose)	std::cout<<"Generated covariance matrix is also positive semi-definite."<<std::endl;
	}

	bkgSpec.compressVector();


	return 0;
}



int SBNchi::load_bkg(){
	lastChi_vec.clear();
	lastChi_vec.resize(num_bins_total_compressed, std::vector<double>( num_bins_total_compressed,0) );

	bkgSpec.compressVector();
	bkgSpec.calcFullVector();

	TMatrixT <double> McI(num_bins_total_compressed,num_bins_total_compressed);
	// Fill systematics from pre-computed files
	TMatrixT <double> Msys(num_bins_total,num_bins_total);

	Msys = sys_fill_direct();

	MfracCov = Msys;

	// systematics per scaled event
	for(int i =0; i<Msys.GetNcols(); i++)
	{
		for(int j =0; j<Msys.GetNrows(); j++)
		{
			std::cout<<"#: "<<Msys(i,j)<<" "<<bkgSpec.fullVec[i]<<" "<<bkgSpec.fullVec[j]<<std::endl;
			Msys(i,j)=Msys(i,j)*bkgSpec.fullVec[i]*bkgSpec.fullVec[j];
		}
	}

	// Fill stats from the back ground vector
	TMatrixT <double> Mstat(num_bins_total,num_bins_total);
	stats_fill(Mstat, bkgSpec.fullVec);


	//And then define the total covariance matrix in all its glory
	TMatrixT <double > Mtotal(num_bins_total,num_bins_total);

	if(stat_only){
		std::cout<<"SBNchi::load_bkg() || Using stats only in covariance matrix"<<std::endl;
		Mtotal = Mstat;
	}else{
		//std::cout<<"SBNchi::load_bkg() || Using stats+sys in covariance matrix"<<std::endl;
		Mtotal = Mstat + Msys;
	}
	// Now contract back the larger antimatrix
	TMatrixT<double > Mctotal(num_bins_total_compressed,num_bins_total_compressed);

	collapse_layer3(Mtotal, Mctotal);

	matrix_collapsed = Mctotal;
	vMc = to_vector(Mctotal);
	// just to hold determinant
	double invdet=0; 

	// Bit o inverting, root tmatrix seems perfectly and sufficiently fast for this, even with anti_mode
	McI = Mctotal.Invert(&invdet);

	// There is currently a bug, somehow a memory leak perhaps. converting the TMatrix to a vector of vectors fixes it for now. 
	vMcI = to_vector(McI);

	return 0;

}


TMatrixT<double> * SBNchi::getCompressedMatrix(){
	TMatrixT<double> * tmp = new TMatrixT<double>(num_bins_total_compressed,num_bins_total_compressed);
	for(int i=0; i<num_bins_total_compressed;i++){
		for(int j=0; j<num_bins_total_compressed;j++){
			(*tmp)(i,j) = vMc.at(i).at(j);
		}
	}

	return tmp;
}

double SBNchi::CalcChiLog(SBNspec *sigSpec){
	double tchi = 0;	

	if(sigSpec->compVec.size()==0){
		std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
		sigSpec->compressVector();
	}

	for(int i =0; i<num_bins_total_compressed; i++){
		for(int j =0; j<num_bins_total_compressed; j++){
			lastChi_vec.at(i).at(j) =(bkgSpec.compVec[i]-sigSpec->compVec[i])*vMcI[i][j]*(bkgSpec.compVec[j]-sigSpec->compVec[j] ); 
			tchi += lastChi_vec.at(i).at(j);
		}
	}

	double absDetM = log(fabs(matrix_collapsed.Determinant())); 

	lastChi = tchi+absDetM;
	return tchi+absDetM;
}

//Use this one
double SBNchi::CalcChi(SBNspec *sigSpec){
	double tchi = 0;	

	if(sigSpec->compVec.size()==0){
	if(isVerbose)	std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
		sigSpec->compressVector();
	}

	int k=0;

	/* THIS IS SIMPLIER CHI
	for(int i =0; i<num_bins_total_compressed; i++){
		double ank = pow(bkgSpec.compVec.at(i)-sigSpec->compVec.at(i),2.0)*vMcI.at(i).at(i);
		tchi+=ank;
		//std::cout<<"AGHR: "<<i<<"\t\t"<<ank<<"\t\t"<<tchi<<"\t\t"<<bkgSpec.compVec.at(i)<<"\t\t"<<sigSpec->compVec.at(i)<<"\t\t"<<vMcI.at(i).at(i)<<"\t\t"<<1.0/vMc.at(i).at(i)<<" "<<vMc.at(i).at(i)<<" "<<MfracCov(i,i)<<std::endl;
				
	}
	return tchi;
	*/
	
	for(int i =0; i<num_bins_total_compressed; i++){
		for(int j =0; j<num_bins_total_compressed; j++){
			k++;
		//	if(i!=j && vMcI.at(i).at(j) != 0){
		//		std::cout<<"ITS NOT SYM"<<std::endl;
		//	} 
			if(i==j && vMcI.at(i).at(j)<0){
				std::cout<<"ERROR: SBNchi::CalcChi || diagonal of inverse covariance is negative!"<<std::endl;
			}


			lastChi_vec.at(i).at(j) =(bkgSpec.compVec.at(i)-sigSpec->compVec.at(i))*vMcI.at(i).at(j)*(bkgSpec.compVec.at(j)-sigSpec->compVec.at(j) ); 
			tchi += lastChi_vec.at(i).at(j);
			 //std::cout<<"AGHR: "<<k<<" "<<i<<" "<<j<<" "<<tchi<<" m: "<<vMc.at(i).at(j)<<" minv: "<<vMcI.at(i).at(j)<<" bi "<<bkgSpec.compVec.at(i)<<" si "<<sigSpec->compVec.at(i)<<" bj "<<bkgSpec.compVec.at(j)<<" sj "<<sigSpec->compVec.at(j)<<std::endl;
		}
	}

	lastChi = tchi;
	return tchi;
}


//Obsoloete version
double SBNchi::CalcChi(SBNspec sigSpec){
	double tchi = 0;	

	if(sigSpec.compVec.size()==0){
		std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
		sigSpec.compressVector();
	}


	for(int i =0; i<num_bins_total_compressed; i++){
		for(int j =0; j<num_bins_total_compressed; j++){
			lastChi_vec.at(i).at(j) =(bkgSpec.compVec[i]-sigSpec.compVec[i])*vMcI[i][j]*(bkgSpec.compVec[j]-sigSpec.compVec[j] ); 
			tchi += lastChi_vec.at(i).at(j);
		}
	}

	lastChi = tchi;
	return tchi;
}


double SBNchi::CalcChi(std::vector<double> sigVec){
	double tchi = 0;	

	if(sigVec.size() != num_bins_total_compressed ){
		std::cerr<<"ERROR: SBNchi::CalcChi(std::vector<double>) ~ your inputed vector does not have correct dimensions"<<std::endl;
		std::cerr<<"sigVec.size(): "<<sigVec.size()<<" num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;
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

double SBNchi::CalcChi(SBNspec *sigSpec, SBNspec *obsSpec){
	double tchi=0;
	if(sigSpec->compVec.size()==0){
		std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
		sigSpec->compressVector();
	}

	if(obsSpec->compVec.size()==0){
		std::cout<<"WARNING: SBNchi::CalcChi, inputted obsSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
		obsSpec->compressVector();
	}


	for(int i =0; i<num_bins_total_compressed; i++){
		for(int j =0; j<num_bins_total_compressed; j++){
			tchi += (obsSpec->compVec[i]-sigSpec->compVec[i])*vMcI[i][j]*(obsSpec->compVec[j]-sigSpec->compVec[j] );
		}
	}

	lastChi = tchi;
	return tchi;

}




void SBNchi::fake_fill(TMatrixT <double> &M){

	TRandom3 *rangen = new TRandom3(0);

	//Fills a square matrix of dim matrix_size with random numbers for now.



	int matrix_size=M.GetNrows();

	if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}

	for(int i=0; i<matrix_size; i++){
		for (int j = i;j<matrix_size;j++){
			M(i,j)=rangen->Uniform(0,1);
			M(j,i)=M(i,j);
		}

	}

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

		//This NEEDS to be removed soon
		//This was just for wierd MiniBooNE run
		//if(i>=11 && i< 30) continue;
		//if(i>=41) continue;
			M(i,i) = diag.at(i);	

	}



	return ;
}



TMatrixT<double> SBNchi::sys_fill_direct(){
	return sys_fill_direct(correlation_matrix_rootfile, correlation_matrix_name);
}


TMatrixT<double > SBNchi::sys_fill_direct(std::string rootname, std::string matname){

	std::cout<<"SBNchi::sys_fill_direct || filling from "<<rootname<<std::endl;

	TMatrixT<double> temp2(num_bins_total,num_bins_total);
	TFile *fm= new TFile(rootname.c_str());

	TMatrixT<float> * temp = (TMatrixT <float>* )fm->Get(matname.c_str());
	//TMatrixT<double> * temp = (TMatrixT <double>* )fm->Get(matname.c_str());


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

	std::cout<<"SBNchi::sys_fill_direct || loaded with dim : "<<temp2.GetNcols()<<" "<<temp2.GetNrows()<<std::endl;

	fm->Close();
	delete fm;

	if(temp2.IsSymmetric()){
		if(isVerbose)std::cout<<"Inputted fracCov covariance matrix is symmetric"<<std::endl;
	}else{
		std::cerr<<"ERROR: SBNchi::sys_fill_direct, Msys input is not symmetric!"<<std::endl;
		//exit(EXIT_FAILURE);
	}




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

	for(int ic = 0; ic < num_channels; ic++){ 	 //Loop over all rows
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

TH2D SBNchi::getChiogram(){
	TH2D tmp("","",num_bins_total_compressed,0, num_bins_total_compressed ,num_bins_total_compressed,0, num_bins_total_compressed);

	for(int i =0; i<num_bins_total_compressed; i++){
		for(int j =0; j<num_bins_total_compressed; j++){

			tmp.SetBinContent(i+1, j+1, lastChi_vec.at(i).at(j));
		}

	}


	return tmp;


}


//This one varies the input comparative spectrum, and as sucn has  only to calculate the Msys once
TH1D SBNchi::toyMC_varyInput(SBNspec *specin, int num_MC){
	double center = this->CalcChi(specin);
	int nlower=0;
	
	TRandom3 *rangen = new TRandom3(0);

	TH1D ans("","",100,0,50);
	//So save the core one that we will sample for
	ans.GetXaxis()->SetCanExtend(kTRUE);
	isVerbose = false;
	for(int i=0; i < num_MC;i++){

		SBNspec tmp = *specin;
		tmp.poissonScale(rangen);
		tmp.compressVector(); //this line important isnt it!
		//tmp.printFullVec();

		double thischi = this->CalcChi(&tmp);
		ans.Fill(thischi);
		if(thischi<=center)nlower++;
			
		if(i%100==0) std::cout<<"SBNchi::toyMC_varyInput(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
	}
	std::cout<<"pval: "<<nlower/(double)num_MC<<std::endl;

	isVerbose = true;
	return ans;


}


//This one varies the input comparative spectrum, and as sucn has  only to calculate the Msys once
std::vector<double> SBNchi::toyMC_varyInput_getpval(SBNspec *specin, int num_MC, std::vector<double> chival){
	std::vector<int> nlower(chival.size(),0);
	

	TRandom3 *rangen = new TRandom3(0);

	TH1D ans("","",100,0,50);
	//So save the core one that we will sample for
	ans.GetXaxis()->SetCanExtend(kTRUE);
	isVerbose = false;
	for(int i=0; i < num_MC;i++){

		SBNspec tmp = *specin;
		tmp.poissonScale(rangen);
		tmp.compressVector(); //this line important isnt it!
		//tmp.printFullVec();

		double thischi = this->CalcChi(&tmp);
		ans.Fill(thischi);
		
	for(int j=0; j< chival.size(); j++){
			if(thischi>=chival.at(j)) nlower.at(j)++;
		}
			
		if(i%1000==0) std::cout<<"SBNchi::toyMC_varyInput(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
	}
	std::vector<double> pval;
	for(auto n: nlower){
		pval.push_back(n/(double)num_MC);

	}

	isVerbose = true;
	return pval;


}




//This one varies the core spectrum, and as sucn has to recalculate the Msys each stem
TH1D SBNchi::toyMC_varyCore(SBNspec *specin, int num_MC){
	double center = this->CalcChi(specin);
	int nlower=0;
	
	TRandom3 *rangen = new TRandom3(0);

	TH1D ans("MCans","MCans",100,center-100,center+200);
	//So save the core one that we will sample for
	SBNspec core  = bkgSpec;
	
	isVerbose = false;
	for(int i=0; i<num_MC;i++){

		SBNspec tmp = core;
		tmp.poissonScale(rangen);
		this->reload_core_spec(&tmp);
		double thischi = this->CalcChi(specin);
		ans.Fill(thischi);
		if(thischi<=center)nlower++;
			
		if(i%1000==0) std::cout<<"SBNchi::toyMC_varyCore(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
	}
	std::cout<<"pval: "<<nlower/(double)num_MC<<std::endl;

	isVerbose = true;
	return ans;


}


