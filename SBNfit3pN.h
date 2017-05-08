#ifndef SBNFIT3PN_H_
#define SBNFIT3PN_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <utility>

#include "SBNosc.h"
#include "SBNfit.h"

#include <TH1D.h>
#include <string>
#include <TF1.h>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "Math/GSLMinimizer.h"
namespace sbn{

class SBNfit3pN : public SBNfit {
	
	protected:
	public:
	
	SBNosc sigOsc;
	double MinimizerCalcChi(const double * X);
	
	SBNfit3pN(SBNosc,SBNosc,int);

};


class SBNfit3p1 : public SBNfit {
	
	protected:
	public:
	
	SBNosc sigOsc;
	double MinimizerCalcChi(const double * X);
	
	SBNfit3p1(SBNosc,SBNosc,int);

};




};
#endif
