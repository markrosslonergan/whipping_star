#ifndef SBNFIT3pN_H_
#define SBNFIT3pN_H_

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
namespace SBNFIT{

class SBNfit3pN : public SBNfit {
	
	protected:
	public:
	
	SBNosc sigOsc;
	double minim_calc_chi(const double * X);
	
	SBNfit3pN(SBNosc,SBNosc,int);

};

};
#endif
