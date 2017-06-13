#include "SBNprob.h"
using namespace sbn;


SBNprob::SBNprob(int dim) : hamiltonian(dim), hamil_kin(dim), potential(dim), UtVU(dim), U(dim), Uconj(dim) {
	
	//eV2GeV = 1e-9;	
	//km2invGeV = 5.06842e18;
	rho = 2.8;
	degree = 3.14159/180.0;

	conversion_parameter = 5.06842;
	dimension=dim;

	Vcc= 2.0*0.27*rho/(1000.0);
	Vnc= -0.5*Vcc;
	Nneutrino = dim;

	t12 = 30*degree;
	t23 = 44*degree;
	t13 = 8*degree;
	t14 = 20*degree;
	t24 = 10*degree;
	t34 = 15*degree;

	d34=34;
	d24=24;
	d13=13;

	Dm21=7.5*pow(10,-5);//7.5e-5;
	Dm31=2.552*pow(10,-3);//e-3;
	Dm41=1;

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
	
	
	potential.real(0,0) = Vcc;
	potential.real(dim-1,dim-1)=-Vnc;


	UtVU=U;
	UtVU.mult(&potential);
	UtVU.mult(&Uconj);

	//Checked this far


}

double SBNprob::probabilityMatterExact(int a, int b, double E, double L ){

	complex_matrix S0(Nneutrino);

	hamiltonian = hamil_kin;
	hamiltonian.mult(conversion_parameter/(2.0*E));
//	UtVU.multI();//artifact
	hamiltonian.add(&UtVU);



	//Blarg, hermitian->antihermitian...  Using the fact Exp[-I M] = Cos[M]-I Sin[M], even for matricies
	//And calculate the matrix exponant 
	std::vector<double> eigenval;
	complex_matrix eigenvec(Nneutrino);
	complex_matrix eigenvecTr(Nneutrino);

	S0=hamiltonian;
	S0.matrixExpTest(L, &eigenval, &eigenvec);

	//exponant is now S0 who cares how it got there

	complex_matrix ans(Nneutrino);
	eigenvecTr = eigenvec;
	eigenvecTr.hermitianConjugate();

	ans = Uconj; //should be eigenvecTr
	ans.mult(&S0);
	ans.mult(&U); //should be eigenvec

	double re = ans.real(b,a);
	double im = ans.imag(b,a);

	return re*re+im*im;
	

};

