#include "SBNprob.h"
using namespace sbn;

SBNprob::SBNprob(int dim,std::vector<double> angles, std::vector<double> phases, std::vector<double> mass ) :hamiltonian(dim), hamil_kin(dim), potential(dim), UtVU(dim), U(dim), Uconj(dim),U_H0_Ut(dim){
	dimension=dim;
	Nneutrino = dim;
	degree = 3.14159/180.0;
	conversion_parameter = 5.06842*1e9;//this is frim km to inv ev
	
	this->setParameters(angles,phases,mass);
	rho = 2.8;

	useMatterEffect = true;
	useNCMatterEffect = true;
	useAntiNeutrino = false;
	this->init();
}
int SBNprob::setParameters(std::vector<double> angles, std::vector<double> phases, std::vector<double> mass){

	t12=angles.at(0)*degree;
	t23=angles.at(1)*degree;
	t13=angles.at(2)*degree;
	
	t14=angles.at(3)*degree;
	t24=angles.at(4)*degree;
	t34=angles.at(5)*degree;

	d13=phases.at(0)*degree;
	d24=phases.at(1)*degree;
	d34=phases.at(2)*degree;

	Dm21=mass.at(0);
	Dm31=mass.at(1);
	Dm41=mass.at(2);

	Ms1 = fabs(Dm21);
	Ms2 = fabs(Dm21)+Dm21;
	Ms3 = fabs(Dm21)+Dm31;
	Ms4 = Dm41; 

	this->init();

		return 0;
}

int SBNprob::init(){

	//eV2GeV = 1e-9;	
	
	double Ye=0.4957;
	//formula from KOPP theisis
	Vcc= 7.56e-14*rho*Ye;// want potential to be in this unitless way for now;


	Vnc= -0.5*Vcc;

	if(!useNCMatterEffect){
		Vnc =0.0;
	}

	complex_matrix R12(Nneutrino);
	complex_matrix R13(Nneutrino);
	complex_matrix R23(Nneutrino);
	complex_matrix R14(Nneutrino);
	complex_matrix R24(Nneutrino);
	complex_matrix R34(Nneutrino);

	R12.setRotation(1,2,t12);
	R13.setComplexRotation(1,3,t13,d13);
	R23.setRotation(2,3,t23);

	R14.setRotation(1,4,t14);
	R24.setComplexRotation(2,4,t24,d24);
	R34.setComplexRotation(3,4,t34,d34);


	U.setIdentity();
	U.mult(&R34);
	U.mult(&R24);
	U.mult(&R23);
	U.mult(&R14);
	U.mult(&R13);
	U.mult(&R12);


	//std::cout<<"PRINT U"<<std::endl;

	for(int i=0; i<4; i++){
	//std::cout<<"("<<U.real(i,0)<<" "<<U.imag(i,0)<<") ("<<U.real(i,1)<<" "<<U.imag(i,1)<<") ("<<U.real(i,2)<<" "<<U.imag(i,2)<<") ("<<U.real(i,3)<<" "<<U.imag(i,3)<<")"<<std::endl; 

	}

	Uconj=U;
	Uconj.hermitianConjugate();

	std::vector<double> masses = {0,0.5*Dm21,0.5*Dm31,0.5*Dm41};
	//std::cout<<"Msq "<<0<<" "<<Dm21<<" "<<Dm31<<" "<<Dm41<<std::endl;
	hamil_kin.setDiagonal(masses);

	//So we have U, Udagger and H

	U_H0_Ut = U;
	U_H0_Ut.mult(&hamil_kin);
	U_H0_Ut.mult(&Uconj);

return 0;
}

double SBNprob::probabilityMatterExact(int a, int b, double E, double L){
	double ans;
	if(useAntiNeutrino){
		ans = probabilityMatterExact(a,b,-1,E,L);
	}else{
		ans = probabilityMatterExact(a,b,1,E,L);
	}
	return ans; 
}

double SBNprob::probabilityMatterExact(int a, int b, int nuornubar, double Energy, double Length ){
	
	for(int i=0;i<4; i++){
	for(int j=0;j<4; j++){
	//	std::cout<<"UHoTUtr "<<i<<" "<<j<<" real "<<U_H0_Ut.real(i,j)<<" imag "<<U_H0_Ut.imag(i,j)<<std::endl;
	}
	}

	double E = Energy*1e9; //in eV
	double L = Length*conversion_parameter;// in eV^-1 also 
	
	//std::cout<<"LCONVR "<<L<<std::endl;
	
	complex_matrix S(Nneutrino);

	hamiltonian = U_H0_Ut;

	//Are we working with antineutrinos here?
	if(nuornubar<0){
		hamiltonian.conj();
	
		potential.real(0,0)= -1.0*Vcc;
		potential.real(dimension-1,dimension-1)=-Vnc;
	}else{
		potential.real(0,0)= Vcc;
		potential.real(dimension-1,dimension-1)=-Vnc;
	}

	hamiltonian.mult(1.0/E,1.0/E);


	
	for(int i=0;i<4; i++){
	for(int j=0;j<4; j++){
	//	std::cout<<"Hamil "<<i<<" "<<j<<" real "<<hamiltonian.real(i,j)<<" imag "<<hamiltonian.imag(i,j)<<std::endl;
	}
	}

	if(useMatterEffect){
		hamiltonian.add(&potential);
	}
	//std::cout<<"Hamil after adding V, "<<hamiltonian.real(0,0)<<" "<<hamiltonian.imag(0,0)<<std::endl;
	
	//Blarg, hermitian->antihermitian...  Using the fact Exp[-I M] = Cos[M]-I Sin[M], even for matricies
	std::vector<double> eigenval;
	complex_matrix eigenvec(Nneutrino);
	complex_matrix eigenvecTr(Nneutrino);




	//hamiltonian has units of inverse km at this point I believe. 
	// Going to find the eigenvalues and eigenvectors of this hamiltonian. This destorys temp4eigen at the moment.
	complex_matrix temp4eigen(Nneutrino);
	temp4eigen = hamiltonian;
	temp4eigen.getEigenStuff(&eigenval, &eigenvec);

	//std::cout<<"Eigen: "<<eigenval.at(0)<<" "<<eigenval.at(1)<<" "<<eigenval.at(2)<<" "<<eigenval.at(3)<<std::endl;
	
	eigenvecTr = eigenvec;
	eigenvecTr.hermitianConjugate();


	//Now calculate the S matrix in the mass basis in batter. Its diagonal here by definitoon;
	//its zero frm its constructer already
	for(int i=0; i<S.dimension; i++){
		double cphase = -L*eigenval.at(i);
		S.real(i,i) = cos(cphase); 
		S.imag(i,i) = sin(cphase); 
	}

	//Now lets transform back to the flavour basis using Q:=eigenvec
	//T0 is just a multplicative aid
	complex_matrix T0(Nneutrino);
	T0 = eigenvec;
	T0.mult(&S);
	T0.mult(&eigenvecTr);

	S=T0;

	double re = S.real(abs(a),abs(b));
	double im = S.imag(abs(a),abs(b));

	return re*re+im*im;


};




double SBNprob::probabilityVacuumExact(int a, int b, double E, double L ){
	return probabilityVacuumExact(a,b,1,E,L);
}

double SBNprob::probabilityVacuumExact(int a, int b, int nuornubar, double E, double L ){
	useMatterEffect =false;
	double ans = probabilityMatterExact(a,b,nuornubar,E,L);
	useMatterEffect = true;
	return ans;
}



double SBNprob::probabilityMatterExactSmear(int a, int b, double E, double L ,double percen, double n){



	double sigma = percen*E/sqrt(E)+0.05*E;

	double low = E-4*sigma;
	if (low<0) low=0.01;

	double high = E+4*sigma;

	double step = fabs(low-high)/n;

	//std::cout<<E<<" low: "<<low<<" high: "<<high<<" step: "<<step<<" sigma: "<<sigma<<std::endl;
	double avg_prob = 0.0;
	for(double et = low; et<=high; et+=step){
		double tmp1 =	probabilityMatterExact(a,b,et,L);
		double tmp2 =	gaussian(et, E, sigma);
		double tmp3 =	tmp1*tmp2*step;
		//std::cout<<E<<" "<<et<<" "<<tmp1<<" "<<tmp2<<" "<<tmp3<<" "<<avg_prob<<std::endl;
		avg_prob += tmp3;	
	}

	return avg_prob;

};


double SBNprob::gaussian(double x, double mean, double sigma){
	double pi=3.14159;
	return 1.0/(sqrt(2*pi)*sigma)*exp( -pow(x-mean,2.0)/(2.0*sigma*sigma));


}

int SBNprob::plotProbabilityMatter(int a, int b, double EminT, double EmaxT, double L, double percen, double n, std::ofstream *filestream){

	std::cout<<"Starting "<<a<<" "<<b<<" plot. L: "<<L<<" VNC: "<<Vnc<<" VCC: "<<Vcc<<" useMatter? "<<useMatterEffect<<" Neutrino or Anti?: "<<useAntiNeutrino<<std::endl;

	double Emin = 1.5*EminT;//(1-0.5*percen)*EminT;
	double Emax = 1.5*EmaxT;//(1.5*percen)*EmaxT;

	int ua=abs(a);
	int ub=abs(b);

	std::vector<double> prob_vec_exact;
	std::vector<double> E_vec;
	std::vector<double> step_vec;
	double step = fabs(Emin-Emax)/n;

	for(double ee=Emin; ee<=Emax; ee+=step){
		double Ee = pow(10,ee);
		prob_vec_exact.push_back( probabilityMatterExact(ua-1,ub-1,Ee,L));
		E_vec.push_back(Ee);
	}

	step_vec.push_back(  fabs( E_vec.at(0)-E_vec.at(1)));
	for(int i=1; i< E_vec.size(); i++){
		step_vec.push_back( fabs( E_vec.at(i)-E_vec.at(i-1)));
	}

	for(int i=0; i< prob_vec_exact.size(); i++){
		if(E_vec.at(i) < pow(10,EminT) || E_vec.at(i)> pow(10,EmaxT)) continue;

		double avg_prob = 0.0;
		double sigma = 	percen*E_vec.at(i)+0.05;///sqrt(E);

		for(int j=0; j< prob_vec_exact.size(); j++){
			
			double tmp1 =   prob_vec_exact.at(j);
			double tmp2 =	gaussian(E_vec.at(j), E_vec.at(i), sigma);
			avg_prob += tmp1*tmp2*step_vec.at(j);	
		//	std::cout<<E_vec.at(i)<<" "<<prob_vec_exact.at(i)<<" "<<E_vec.at(j)<<" "<<prob_vec_exact.at(j)<<" "<<sigma<<" "<<tmp1<<" "<<tmp2<<" "<<avg_prob<<std::endl;
		}
		*filestream<<a<<b<<" "<<L<<" "<<E_vec.at(i)<<" "<<prob_vec_exact.at(i)<<" "<<avg_prob<<std::endl;
	}


	return 0;
}


/*
double SBNprob::probabilityGlobes(int a, int b, int panti, double E, double L ){
	double ans = 0;
	{
	char* dunno;
	glbInit(dunno);


	glb_params globes_params = glbAllocParams();
	glbDefineParams(globes_params,t12,t13, t23, d13,Dm21,Dm31);
	glbSetDensityParams(globes_params,1.0,GLB_ALL);
	glbSetOscillationParameters(globes_params);

  	glbSetRates();
  
	ans =  glbConstantDensityProbability(a+1,b+1,panti, E,L, rho);

	glbFreeParams(globes_params);

	}
	return ans;

}
*/

int SBNprob::setMatterEffect(bool in){
	useMatterEffect = in;
	this->init();
	return 0;
}

int SBNprob::setNCMatterEffect(bool in){
	useNCMatterEffect = in;
	this->init();
	return 0;
}

int SBNprob::setAntiNeutrinoMode(bool in){
       useAntiNeutrino = in;
       //d13 = -d13;
       //d24 = -d24;
       //d34 = -d34;
       this->init();
       return 0;

}




SBNprob::SBNprob(int dim) : hamiltonian(dim), hamil_kin(dim), potential(dim), UtVU(dim), U(dim), Uconj(dim), U_H0_Ut(dim) {
	dimension=dim;
	Nneutrino = dim;
	degree = 3.14159/180.0;
	conversion_parameter = 5.06842;

	t12 = 30*degree;
	t23 = 44*degree;
	t13 = 8*degree;
	t14 = 15*degree;
	t24 = 10*degree;
	t34 = 20*degree;

	d34=0;
	d24=0;
	d13=0;

	Dm21=7.5*pow(10,-5);//7.5e-5;
	Dm31=2.552*pow(10,-3);//e-3;
	Dm41=1;

	rho = 2.8;

	useMatterEffect = true;
	useNCMatterEffect = true;
	this->init();
}


