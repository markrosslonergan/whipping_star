#include "SBNfit3pN.h"
using namespace sbn;

/****************************************************************
 *		Generic 3+N 
 * *************************************************************/


SBNfit3pN::SBNfit3pN(SBNosc inBk, SBNosc inSg, int npa) : SBNfit(inBk,inSg,npa), sigOsc(inSg) {

}


double SBNfit3pN::MinimizerCalcChi(const double * X){
	num_func_calls++;
	SBNosc tempOsc = sigOsc; 

		
		
		double imn[3] = {X[0],X[1],X[2]};
		double iue[3] = {X[3],X[4],X[5]};
		double ium[3] = {X[6],X[7],X[8]};
		double iph[3] = {X[9],X[10],X[11]};
		neutrinoModel signalModel(imn,iue,ium,iph);
					
		tempOsc.load_model(signalModel);
		std::vector<double> ans = tempOsc.Oscillate();
	

		lastChi =this->CalcChi(ans);
	return lastChi;

}


/****************************************************************
 *			3+1	Only 
 * *************************************************************/


SBNfit3p1::SBNfit3p1(SBNosc inBk, SBNosc inSg, int npa) : SBNfit(inBk,inSg,npa), sigOsc(inSg) {

}

double SBNfit3p1::MinimizerCalcChi(const double * X){
	num_func_calls++;
	SBNosc tempOsc = sigOsc; 
	
	neutrinoModel signalModel(X[0],X[1],X[2]);
	signalModel.numsterile=1;		

	tempOsc.load_model(signalModel);

	std::vector<double> ans = tempOsc.Oscillate();
	
	lastChi =this->CalcChi(ans);

	return lastChi;
}

