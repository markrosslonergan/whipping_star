#ifndef PROB_H_
#define PROB_H_

#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <TRandom3.h>
#include <TMatrixT.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


namespace sbn{

/*************************************************************
 *************************************************************
 *	ToDO:
 *	Actually:  add tau to neutrinoModel. Simple and quick, but have to adjust constructors
 ************************************************************
 ************************************************************/

double smearEnergy(double E, double percen, TRandom3 * rangen);


struct complex_matrix{
	complex_matrix(int dim);
	
	int dimension;
	TMatrixT<double> real;
	TMatrixT<double> imag;
	int mult(complex_matrix* in);
	int multI();
	int mult(double );
	int mult(double ,double);
	int add(complex_matrix* in );
	int conj();



	std::vector<double> matrixExp();
	std::vector<double> matrixExpTest(double L, std::vector<double> *, complex_matrix *);


	std::vector<double> getEigenStuff(std::vector<double> *, complex_matrix *);

	int setRotation(int,int, double);
	int setComplexRotation(int,int,double,double);
	int setDiagonal(std::vector<double> ms);
	int setIdentity();

	int hermitianConjugate();
	int transpose();

	int print();
};
struct neutrinoModel{
	double mNu[3], Ue[3], Um[3], phi[3];
	double dm41Sq, dm51Sq, dm61Sq, dm54Sq, dm64Sq, dm65Sq;
	std::vector< std::vector < std::complex<double> > >  U; 
	//std::vector<std::vector<double>> dm;

	int numsterile;

	//constructors!! Should overload these immensely for  3+1, 3+2, 3+3 and NULL
	neutrinoModel();
	neutrinoModel(double * mn, double * ue, double *um, double *ph);
	neutrinoModel(double m4, double ue4, double um4);

	double UUem;
	double UUme;
	double UUmm;
	double UUee;

 	void printall();

	void zero();
	void difference();

	//Might as well have the oscProb in the neutrioModel class as all the things are aready here.
	double oscProb(int init, int fin, double Ev, double L);
	double oscProb_dis(int a, double  Ev, double L);
	double oscProb_app(int a,int b, double  Ev, double L);

	double oscProbSin(double Ev, double L);
	double oscProbSinSq(double Ev, double L);


	double oscAmp(int a, int b, int which_dm, int sqornot);
	double oscAmp_dis(int a, int which_dm);
	double oscAmp_app(int a, int b, int which_dm, int sqornot);
};


// Simple 3+1 testing probabilities, completely redundant now.
double Pmue(double L, double E, double dm, double sin2);
double Pmm(double L, double Ev, double Dm, double sinSq2thmm);


//arbitraty U matrix version overload?
double oscProb(int a, int b, double Ev, double L,  std::vector< std::vector < std::complex<double> > >    U, std::vector<std::vector<double> > dm);



};
#endif
