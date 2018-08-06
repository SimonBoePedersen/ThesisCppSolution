#pragma once

template<typename T>
class creditObject
{
	// Cox model setup with piecewise constant hazard rates
	// Main method is default_prob which give default probability for a given time interval
public:
	std::vector<T> bps;
	std::vector<T> spreads;
	std::vector<T> tenor;
	std::vector<T> surv;
	std::vector<T> hazard_rates;

	creditObject(std::vector<T>& spreads, std::vector<T>& tenor, BaseRateModel<T>* mdl) 
		: spreads(spreads), tenor(tenor) {
		surv.resize(tenor.size());
		hazard_rates.resize(tenor.size());
		T deduce = 0.0;
		for (size_t i = 0; i < tenor.size(); i++)
		{
			bps.push_back(mdl->bond_price(0.0, tenor[i], mdl->r0));
		}
		surv[0] = 1.0;
		for (size_t ii = 1; ii < tenor.size(); ii++)
		{
			T upper_sum = 0.0;
			T stepconst = spreads[ii] / (1.0 - 0.4);
			for (size_t kk = 1; kk < ii; kk++)
			{
				upper_sum = upper_sum +  stepconst * bps[kk] * (tenor[kk] - tenor[kk - 1])*surv[kk - 1] - bps[kk] * (surv[kk - 1] - surv[kk]);
			}
			T denom = bps[ii] * (1 + (tenor[ii] - tenor[ii - 1]) * stepconst);
			surv[ii] = (surv[ii - 1] * bps[ii] - upper_sum) / denom;

		}
		surv_to_hazard(deduce);
	};
	~creditObject() {};
	
	template<typename T>
	T default_prob(T& start, T& mat) {
		// Returns the porbability of default in the interval (start, mat)
		T intmat = 0;	T intsta = 0;

		int i = 1;
		while (tenor[i] < mat) {
			intmat = intmat + (tenor[i] - tenor[i - 1]) * hazard_rates[i];
			i++;
		}
		intmat = intmat + (mat - tenor[i - 1]) * hazard_rates[i];

		int j = 1;
		while (tenor[j] < start) {
			intsta = intsta + (tenor[j] - tenor[j - 1]) * hazard_rates[j];
			j++;
		}
		intsta = intsta + (start - tenor[j - 1]) * hazard_rates[j];

		return exp(-intsta) - exp(-intmat);
	}

	template<typename T>
	void surv_to_hazard(T deduce)
	{
		// deduce is a poor solution to let the compiler know which type T is.
		hazard_rates[0] = 0.0;
		for (size_t i = 1; i < hazard_rates.size(); i++)
		{
			hazard_rates[i] = -log(surv[i] / surv[i - 1]) / (tenor[i] - tenor[i - 1]);
		}
	}

};

