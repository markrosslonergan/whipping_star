#include "SBNfit3pN.h"

SBNfit3pN::SBNfit3pN(SBNosc inBk, SBNosc inSg, int npa) : SBNfit(inBk,inSg,npa), sigOsc(inSg) {

}


double SBNfit3pN::minim_calc_chi(const double * X){
	BFncalls++;
	SBNosc tempOsc = sigOsc; 

		
		double imn[3] = {X[0],X[1],X[2]};
		double iue[3] = {X[3],X[4],X[5]};
		double ium[3] = {X[6],X[7],X[8]};
		double iph[3] = {X[9],X[10],X[11]};
		neutrinoModel signalModel(imn,iue,ium,iph);
		signalModel.numsterile = 1;
					
		tempOsc.load_model(signalModel);
		std::vector<double> ans = tempOsc.oscillate();
	

		double ch =this->calc_chi(ans);

	lastChi = ch;
	std::cout<<X[0]<<" "<<X[3]<<" "<<X[6]<<" "<<lastChi<<std::endl;
	return ch;

}

