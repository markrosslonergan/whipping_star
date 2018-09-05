#include "SBNfit.h"
using namespace sbn;



SBNfit::SBNfit(SBNspec inBk, SBNspec inSg, TMatrixD mat, int npar) : SBNchi(inBk, mat), sigOsc(inSg), num_params(npar) {

	for(int i =0; i< num_params; i++){
		f_is_fixed.push_back(0);
		f_param_names.push_back("");
		f_initial_values.push_back(0.5);
		f_upper_values.push_back(1);
		f_lower_values.push_back(0);
		f_step_sizes.push_back(0.01);
	}

	f_minimizer_mode ="GSLMultiMin"; //"GSLSimAn"
	f_minimizer_algo= "BFGS2";

	num_func_calls = 0;
}


SBNfit::SBNfit(SBNspec inBk, SBNspec inSg, int npar) : SBNchi(inBk), sigOsc(inSg), num_params(npar) {


	for(int i =0; i< num_params; i++){
		f_is_fixed.push_back(0);
		f_param_names.push_back("");
		f_initial_values.push_back(0.5);
		f_upper_values.push_back(1);
		f_lower_values.push_back(0);
		f_step_sizes.push_back(0.01);
	}

	f_minimizer_mode ="GSLMultiMin"; //"GSLSimAn"
	f_minimizer_algo= "BFGS2";

	num_func_calls = 0;
}

int SBNfit::load_signal(SBNspec inSg){
	sigOsc = inSg;			
	return 0;
}

int SBNfit::initialize_norm(std::vector< std::pair< std::string, int > > vecIn  ){
	vec_scales = vecIn;

	return 0;
}

double SBNfit::MinimizerCalcChi(const double * X){
	num_func_calls++;
	fOsc = sigOsc; 

	for(auto& v: vec_scales){
		fOsc.Scale(v.first, X[v.second] );
	}	

	fOsc.compressVector();	

	lastChi =this->CalcChi(fOsc);

	return lastChi;

}
	
double SBNfit::Minimize(){	
	num_func_calls=0;
   
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(f_minimizer_mode, f_minimizer_algo);

 	min->SetMaxIterations(20000);  // for GSL
	//min->SetTolerance(200000); //times 4 for normal
	min->SetPrintLevel(1);
	min->SetPrecision(0.001);//times 4 for normal


        ROOT::Math::Functor f( this, &SBNfit::MinimizerCalcChi, num_params); 

   	min->SetFunction(f);

   	for(int i=0;i<num_params;i++){
		if(f_is_fixed[i]){
	   	min->SetFixedVariable(i,f_param_names[i],f_initial_values[i]);
	} else {
   		min->SetLimitedVariable(i,f_param_names[i],f_initial_values[i], f_step_sizes[i], f_lower_values[i],f_upper_values[i]);
	}

   	}

  	min->Minimize(); 
   
  	const double *xs = min->X();

	  bf_chi= MinimizerCalcChi(xs);;
	  bf_params = xs;
	  return bf_chi;

}

/****************************************************
 ***		Some initial setup things
 * *************************************************/

int SBNfit::setMethod(std::string mode, std::string algo){
	f_minimizer_mode =mode; //"GSLSimAn"
	f_minimizer_algo= algo;


return 0;
}


int SBNfit::setInitialValues(std::vector<double> inv){
	for(int i = 0; i< num_params; i++){
		f_initial_values[i]=inv[i];
	}
	return 0;
}
int SBNfit::setInitialValues(double in){
	for(int i = 0; i< num_params; i++){
		f_initial_values[i]=in;
	}
	return 0;
}



int SBNfit::setUpperValues(std::vector<double>  inv){
	for(int i = 0; i< num_params; i++){
		f_upper_values[i]=inv[i];
	}
	return 0;
}
int SBNfit::setUpperValues(double in){
	for(int i = 0; i< num_params; i++){
		f_upper_values[i]=in;
	}
	return 0;
}



int SBNfit::setLowerValues(std::vector<double>  inv){
	for(int i = 0; i< num_params; i++){
		f_lower_values[i]=inv[i];
	}
	return 0;
}
int SBNfit::setLowerValues(double in){
	for(int i = 0; i< num_params; i++){
		f_lower_values[i]=in;
	}
	return 0;
}

int SBNfit::setStepSizes(std::vector<double>  inv){
	for(int i = 0; i< num_params; i++){
		f_step_sizes[i]=inv[i];
	}
	return 0;
}
int SBNfit::setStepSizes(double in){
	for(int i = 0; i< num_params; i++){
		f_step_sizes[i]=in;
	}
	return 0;
}



int SBNfit::setFixed(std::vector<int>  inv){
	for(int i = 0; i< num_params; i++){
		f_is_fixed[i]=inv[i];
	}
	return 0;
}
int SBNfit::setFixed(int in){
	for(int i = 0; i< num_params; i++){
		f_is_fixed[i]=in;
	}
	return 0;
}


int SBNfit::setNames(std::vector<std::string> inv){
	for(int i = 0; i< num_params; i++){
		f_param_names[i]=inv[i];
	}
	return 0;
}

