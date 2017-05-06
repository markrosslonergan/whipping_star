#include "SBNfit.h"
using namespace SBNFIT;

SBNfit::SBNfit(SBNspec inBk, SBNspec inSg, int npar) : SBNchi(inBk), sigOsc(inSg), Npar(npar) {


	for(int i =0; i< Npar; i++){
		fFixed.push_back(0);
		fNames.push_back("");
		fInitialValues.push_back(0.5);
		fUpperValues.push_back(1);
		fLowerValues.push_back(0);
		fStepSizes.push_back(0.01);
	}

	fMinMode ="GSLMultiMin"; //"GSLSimAn"
	fMinAlgo= "BFGS2";

	BFncalls = 0;
}

int SBNfit::load_signal(SBNspec inSg){
	sigOsc = inSg;			
	return 0;
}

int SBNfit::initialize_norm(std::vector< std::pair< std::string, int > > vecIn  ){
	vecScales = vecIn;

	return 0;
}

double SBNfit::minim_calc_chi(const double * X){
	BFncalls++;
	fOsc = sigOsc; 

	for(auto& v: vecScales){
		fOsc.Scale(v.first, X[v.second] );
	}	

	fOsc.compressVector();	

	double ans =this->calc_chi(fOsc);

	lastChi = ans;
	std::cout<<X[0]<<" "<<X[1]<<" "<<lastChi<<std::endl;
	return ans;

}

int SBNfit::init_minim(std::string mode, std::string algo){
//	min = new ROOT::Math::GSLMinimizer(ROOT::Math::kConjugateFR);	
//	ROOT::Math::Minimizer* min = new ROOT::Math::GSLMinimizer(ROOT::Math::kVectorBFGS2);	
//	min = new ROOT::Math::GSLSimAnMinimizer();
  // 	min->SetMaxIterations(250);  // for GSL
//	min->SetTolerance(0.001); //times 4 for normal
//	min->SetPrintLevel(0);
//	min->SetPrecision(0.0001);//times 4 for normal
	fMinMode =mode; //"GSLSimAn"
	fMinAlgo= algo;



return 0;
}
	
double SBNfit::minimize(){	
	BFncalls=0;
	//ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
//   ROOT::Math::Minimizer* min = 
//   //          ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
//   //   ROOT::Math::Minimizer* min = 
//   //          ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
//   //   ROOT::Math::Minimizer* min = 
//   //          ROOT::Math::Factory::CreateMinimizer("Minuit2", "Scan");
//   //   ROOT::Math::Minimizer* min = 
//   //          ROOT::Math::Factory::CreateMinimizer("Minuit2", "Fumili");
//   //   ROOT::Math::Minimizer* min = 
//   //          ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "ConjugateFR");
//   //   ROOT::Math::Minimizer* min = 
//   //          ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "ConjugatePR");
//   //   ROOT::Math::Minimizer* min = 
//   //          ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS");
//   ROOT::Math::Minimizer* min =   ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS2");
   //   ROOT::Math::Minimizer* min = 
//   //          ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "SteepestDescent");
//   //   ROOT::Math::Minimizer* min = 
//   //          ROOT::Math::Factory::CreateMinimizer("GSLMultiFit", "");
   ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(fMinMode,fMinAlgo);


 	min->SetMaxIterations(20000);  // for GSL
	//min->SetTolerance(200000); //times 4 for normal
	min->SetPrintLevel(1);
	min->SetPrecision(0.001);//times 4 for normal



	std::cout<<"define functor"<<std::endl;
        ROOT::Math::Functor f( this, &SBNfit::minim_calc_chi, Npar); 

	std::cout<<"set functor"<<std::endl;
   	min->SetFunction(f);

	std::cout<<"set variables"<<std::endl;
   	for(int i=0;i<Npar;i++){
		if(fFixed[i]){
	   	min->SetFixedVariable(i,fNames[i],fInitialValues[i]);
	} else {
   		min->SetLimitedVariable(i,fNames[i],fInitialValues[i], fStepSizes[i], fLowerValues[i],fUpperValues[i]);
	}

   	}

	std::cout<<"internam min"<<std::endl;
  min->Minimize(); 
   
  const double *xs = min->X();
  double valAns = minim_calc_chi(xs);
	// or min->minValue();

  BFchi= valAns;
  BFparam = xs;
  return valAns;

}



int SBNfit::setInitialValues(std::vector<double> inv){
	for(int i = 0; i< Npar; i++){
		fInitialValues[i]=inv[i];
	}
	return 0;
}
int SBNfit::setInitialValues(double in){
	for(int i = 0; i< Npar; i++){
		fInitialValues[i]=in;
	}
	return 0;
}



int SBNfit::setUpperValues(std::vector<double>  inv){
	for(int i = 0; i< Npar; i++){
		fUpperValues[i]=inv[i];
	}
	return 0;
}
int SBNfit::setUpperValues(double in){
	for(int i = 0; i< Npar; i++){
		fUpperValues[i]=in;
	}
	return 0;
}



int SBNfit::setLowerValues(std::vector<double>  inv){
	for(int i = 0; i< Npar; i++){
		fLowerValues[i]=inv[i];
	}
	return 0;
}
int SBNfit::setLowerValues(double in){
	for(int i = 0; i< Npar; i++){
		fLowerValues[i]=in;
	}
	return 0;
}

int SBNfit::setStepSizes(std::vector<double>  inv){
	for(int i = 0; i< Npar; i++){
		fStepSizes[i]=inv[i];
	}
	return 0;
}
int SBNfit::setStepSizes(double in){
	for(int i = 0; i< Npar; i++){
		fStepSizes[i]=in;
	}
	return 0;
}



int SBNfit::setFixed(std::vector<int>  inv){
	for(int i = 0; i< Npar; i++){
		fFixed[i]=inv[i];
	}
	return 0;
}int SBNfit::setFixed(int in){
	for(int i = 0; i< Npar; i++){
		fFixed[i]=in;
	}
	return 0;
}


int SBNfit::setNames(std::vector<std::string> inv){
	for(int i = 0; i< Npar; i++){
		fNames[i]=inv[i];
	}
	return 0;
}

