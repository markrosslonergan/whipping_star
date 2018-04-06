#include "SBNprob.h"
using namespace sbn;

SBNprob::SBNprob(int dim,std::vector<double> angles, std::vector<double> phases, std::vector<double> mass ) :hamiltonian(dim), hamil_kin(dim), potential(dim), UtVU(dim), U(dim), Uconj(dim){
	dimension=dim;
	Nneutrino = dim;
	degree = 3.14159/180.0;
	conversion_parameter = 5.06842; //this is Gev->inv KM of course

	this->setParameters(angles,phases,mass);
	rho = 2.8;

	useMatterEffect = true;
	useNCMatterEffect = true;
	useAntiNeutrino = false;
	this->init();
}

SBNprob::SBNprob(int dim) : hamiltonian(dim), hamil_kin(dim), potential(dim), UtVU(dim), U(dim), Uconj(dim) {
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
	this->init();

	return 0;

}

int SBNprob::init(){

	//eV2GeV = 1e-9;	
	
	double Ye=0.4957;
	//formula from KOPP theisis
	Vcc= conversion_parameter*7.56e-14*1e9*rho*Ye;// want potential to be in inv km also.


	if(useAntiNeutrino){
		Vcc=-Vcc;
	}

	//does the VNC flip slight also with antineutrino, check. I believe so.
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



	Uconj=U;
	Uconj.hermitianConjugate();
	std::vector<double> mass_splittings = {0.0,Dm21,Dm31,Dm41};
	hamil_kin.setDiagonal(mass_splittings);

	potential.real(0,0)= Vcc;
	potential.real(dimension-1,dimension-1)=-Vnc;

	UtVU=U;
	UtVU.mult(&hamil_kin);
	UtVU.mult(&Uconj);


return 0;
}

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
       d13 = -d13;
       d24 = -d24;
       d34 = -d34;
       this->init();
       return 0;

}

double SBNprob::probabilityVacuumExact(int a, int b, double E, double L ){

	complex_matrix S0(Nneutrino);

	hamiltonian = UtVU;
	
	hamiltonian.mult(conversion_parameter/(2.0*E));
	
	//Blarg, hermitian->antihermitian...  Using the fact Exp[-I M] = Cos[M]-I Sin[M], even for matricies
	//And calculate the matrix exponant 
	std::vector<double> eigenval;
	complex_matrix eigenvec(Nneutrino);
	complex_matrix eigenvecTr(Nneutrino);

	//hamiltonian has units of inverse km at this point I believe. 
	S0=hamiltonian;
	S0.matrixExpTest(L, &eigenval, &eigenvec);

	//exponant is now S0 who cares how it got there

	complex_matrix ans(Nneutrino);
	eigenvecTr = eigenvec;
	eigenvecTr.hermitianConjugate();

	ans=S0;
	//ans = Uconj; //should be eigenvecTr
	//ans.mult(&S0);
//	ans.mult(&U); //should be eigenvec ?? 

	double re = ans.real(b,a);
	double im = ans.imag(b,a);

	return re*re+im*im;


};



double SBNprob::probabilityMatterExact(int a, int b, double E, double L ){

	complex_matrix S0(Nneutrino);

	hamiltonian = UtVU;
	
	hamiltonian.mult(conversion_parameter/(2.0*E));
	
	if(useMatterEffect){
		hamiltonian.add(&potential);
	}

	//Blarg, hermitian->antihermitian...  Using the fact Exp[-I M] = Cos[M]-I Sin[M], even for matricies
	//And calculate the matrix exponant 
	std::vector<double> eigenval;
	complex_matrix eigenvec(Nneutrino);
	complex_matrix eigenvecTr(Nneutrino);

	//hamiltonian has units of inverse km at this point I believe. 
	S0=hamiltonian;
	S0.matrixExpTest(L, &eigenval, &eigenvec);

	//exponant is now S0 who cares how it got there

	complex_matrix ans(Nneutrino);
	eigenvecTr = eigenvec;
	eigenvecTr.hermitianConjugate();

	ans=S0;
	//ans = Uconj; //should be eigenvecTr
	//ans.mult(&S0);
//	ans.mult(&U); //should be eigenvec ?? 

	double re = ans.real(b,a);
	double im = ans.imag(b,a);

	return re*re+im*im;


};

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






