#pragma once
#include <memory>
#include <vector>
template <class T, class U>
T modulo(T x, U y) {
	T r = x%y;
	return r < 0 ? r + y : r;
}


typedef long long ullong;

class mrg32
{
	// Random generator object with the mrg32k3a methodology. Main method is next() which
	//  returns the next unif[0,1] in the stream of random numbers
public:

	mrg32(long long s11, long long s12, long long s13, long long s21, long long s22, long long s23) {
		setseeds(s11, s12, s13, s21, s22, s23);
	};

	virtual std::unique_ptr<mrg32> clone() const {
		std::unique_ptr<mrg32> rng = std::unique_ptr<mrg32>(new mrg32(12345, 12345, 12345, 12346, 12346, 12346));
		return rng;
	}; //Makes a unique_ptr clone of the object, such that the threads don't "have their fingers on the same one"
	   //Thanks to raii semantics the use of unique_ptr's is easy and "cleans-up after itself"

	void setseeds(long long s11, long long s12, long long s13, long long s21, long long s22, long long s23) {
		x = { s11, s12, s13 };
		y = { s21, s22, s23 };
	};

	double next() 
	{
		//Outputs a single unif[0,1] from the mrg32-stream of pseudo-random
		x.push_back(modulo(aOneTwo * x[1] + aOneThree * x[0], m1));
		y.push_back(modulo(aTwoOne * y[2] + aTwoThree * y[0], m2));
		x.erase(x.begin());
		y.erase(y.begin());

		return (double)(modulo(x[2] - y[2], m1)) / (double)m1;
	};
	
	
	// Below methods are not necessary for this thesis as no parallization is done. 
	// code can be required by request on simonboepedersen@gmail.com

	//void nextSkip(int nSkips);          //Very bad manual skipping by calling next() nSkips times. Used for testing fastSkip!
	//std::vector<double> nextVec(int nSteps); //outputs a vector full of next
	//void powerskip(int nSkips);			//This method skip-ahead, but only in lengths of 2^i for i = 0,1,...
	//void skipAhead(int nSkips);			//This method does fast skip-ahead using the powerskip method and a clever binary representation


private:

	std::vector<long long> x, y; //state variable
	long long aOneTwo = 1403580;
	long long aOneThree = -810728;
	long long aTwoOne = 527612;
	long long aTwoThree = -1370589;
	long long m1 = 4294967087;
	long long m2 = 4294944443;

};
