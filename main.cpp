#include <math.h>
#include <iostream>
#include <Windows.h>
#include <vector>

template<typename T>
T SIGN(T& a, T& b) { return b >= 0.0 ? abs(a) : -abs(b); }

double SIGN(double& a, double& b) { return (b) >= 0.0 ? abs(a) : -abs(a); }
double SQR(double& a) { return a*a; }
double MAX(double& a, double& b) { return a > b ? a : b; }
double MIN(double& a, double& b) { return a < b ? a : b; }

#include <fstream>
#include "gaussians.h"
#include "mrg32.h"

#include "tape.h"
Tape *GlobalTape = new Tape;
#include "vars.h"

Variable MAX(Variable& a, Variable& b) { return a.value > b.value ? a : b; }
Variable MIN(Variable& a, Variable& b) { return a.value < b.value ? a : b; }
Variable MAX(Variable& a, double& b) { return a.value > b ? a : b; }
Variable MAX(double& a, Variable& b) { return MAX(b, a); }
Variable MIN(Variable& a, double& b) { return a.value < b ? a : b; }
Variable MIN(double& a, Variable& b) { return MIN(b, a); }
Variable SQR(Variable& a) { return a*a; }

Variable SIGN(Variable& a, Variable& b) { return b.value >= 0.0 ? abs(a) : -abs(a); }
Variable SIGN(Variable& a, double& b) { return (b) >= 0.0 ? abs(a) : -abs(a); }
Variable SIGN(double& a, Variable& b) {
	Variable aa = a;
	return b.value >= 0.0 ? abs(aa) : -abs(aa);
}

#include "matrix.h"				// General matrix class with inversion
#include "svd.h"				// Algorithm for singular value decomposition

#include "spline.h"				// Defines the spline object for hermitic spline interpolation
#include "timeline.h"			

#include "baseratemodel.h"		// Interface class for short rate models
#include "vasicek.h"			// Vasicek object class
#include "hullwhite.h"			// Hull White object class

#include "derivatives.h"		// Defines a fixed rate bond, floating rate bond and IRS
#include "basisclasses.h"		// Defines the basis functions needed for LSM
#include "cds.h"				// Credit object with piecewise constant hazard rate

#include "value_simulations.h"  // Contains routines for simulation of value matrices of interest rate 
								// swaps with the different methods described in the thesis
#include "full_routines.h"

#include "examples.h"			//This header contains functions which calculates
								//the results needed in the examples of the thesis




int main() {
	// Uncomment to run code examples.
	std::cout << "Running main(). Press enter." << std::endl;
	std::cin.ignore();
	
	//////// Example chapter 2 Hull-White ////////
	//hullWhiteFittingExample();

	//////// Example chapter 3 ////////
	// The amounts of paths has to be changed depending on which example it is
	//example_after_proxy();
	//example_after_poi();
	
	//////// Example Chapter 4 ////////
	// Example IRS delta vector //
	//example_aad_swap_delta_vector();
	//final_example_naive_full_sens_IRS10Y();
	//final_example_naive_cashflow_sens_IRS10Y();
	
	
	//To measure FD performance
	//long time_mean = 0.0;
	//for (size_t i = 0; i < 100; i++)	time_mean += delta_vector_only_finite_difference(false);
	//for (size_t i = 0; i < 100; i++)	time_mean += final_example_naive_full_sens_IRS10Y(false);
	//for (size_t i = 0; i < 100; i++)	time_mean += final_example_naive_cashflow_sens_IRS10Y(false);
	//std::cout << time_mean / 100;

	//Testing non-naive cashflow 8192 performance
	/*
	std::vector<size_t> idx = { 1,8,16,32,64,128,256,512,1024 };
	for (size_t kk = 0; kk < idx.size(); kk++)
	{
		long time_mean = 0.0;
		for (size_t i = 0; i < 100; i++)
		{
			time_mean += final_example_pathwise_cashflow_sens_IRS10Y(idx[kk], false);
		}
		std::cout << time_mean / 100 << std::endl;
	}
	*/
	

	// Portfolio 2 with 5 derivatives
	//portfolio_two_npv_cva_with_doubles();
	//portfolio_two_cashflow(256);
	

	GlobalTape->nodes.clear();
	delete GlobalTape;
	std::cin.ignore();
	return 0;
}