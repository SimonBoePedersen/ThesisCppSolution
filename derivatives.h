#pragma once

// NOTE: These classes are not poorly build for VERY specific use
// Currently they can only be used with 1Y fixed tenor and 6M Float tenor

double TOLERANCE = 0.000001;

template<typename T>
class FixedRateBond
{
	// Fixed rate bond used in swap
	// Fixed rate bond is specified by rate, tenor and maturity
	// Note: poorly specific code. Needs to be generalized.
public:
	FixedRateBond() {};
	FixedRateBond(T& fixed_rate, T& maturity, T& tenor, BaseRateModel<T>* mdl)
		: fixed_rate(fixed_rate), mdl(mdl)
	{
		Timeline<T> time(tenor, maturity + tenor, tenor);

		dates = time;
	};
	~FixedRateBond() { };

	T npv(T& current_time, T& short_rate) {
		// Calculate net present value at time "current_time" and at short rate "short_rate"
		if (current_time > dates.line.back() || abs(current_time - dates.line.back()) < TOLERANCE)
		{
			return 0.0;
		}
		std::vector<T> cvg = remaining_cvg(current_time, dates);
		std::vector<T> disc_times = cumsum(cvg);

		T out = 0;
		for (size_t i = 0; i < cvg.size(); i++)
		{
			out = out + mdl->bond_price(current_time, current_time + disc_times[i], short_rate) * fixed_rate;
		}
		return out;
	};

	BaseRateModel<T>* mdl;
	Timeline<T> dates;
	T fixed_rate;
};



template<typename T>
class FloatingRateBond
{
	// Floating rate bond used in swap
	// Floating rate bond is specified by rate, tenor and maturity
	// Note: poorly specific code. Needs to be generalized.
public:
	FloatingRateBond() {};
	FloatingRateBond(T& maturity, T& tenor, BaseRateModel<T>* mdl, T& short_rate)
		: mdl(mdl)
	{
		Timeline<T> time(tenor, maturity + tenor, tenor);
		dates = time;
		T first_payment = (dates[0]) * mdl->spot_libor(0,dates[0],short_rate); 
		std::vector<T> prep(time.line.size());
		prep[0] = first_payment;
		payments = prep;
	};
	~FloatingRateBond() {};


	void update(T& current_time, T& short_rate) {
		//Updates the payment if current_time %in% dates
		for (size_t i = 0; i < dates.size() - 1; i++)
		{
			if (abs(current_time - dates[i]) < TOLERANCE) {
				payments[i + 1] = (dates[i+1]-dates[i]) *
					mdl->spot_libor(dates[i], dates[i + 1], short_rate);
				break;
			}
		}
	}

	T npv(T& current_time, T& short_rate) {
		// Calculate net present value at time "current_time" and at short rate "short_rate"
		// If we are at a fixing point in time we update the next payment
		if (current_time > dates.line.back() ||
			abs(current_time - dates.line.back()) < TOLERANCE)
		{
			return 0.0;
		}
		std::vector<T> cvg = remaining_cvg(current_time, dates);
		std::vector<T> disc_times = cumsum(cvg);
		for (int kk = 0; kk < dates.line.size(); kk++)
		{
			if (abs(dates[kk] - current_time) < TOLERANCE) {
				update(current_time, short_rate);
				break;
			}
		}
		T out = mdl->bond_price(current_time, current_time + disc_times[0], short_rate) * payments[dates.size()-cvg.size()];
		
		for (size_t i = 1; i < cvg.size(); i++)
		{
			out = out + mdl->bond_price(current_time, current_time + disc_times[i - 1], short_rate) -
				mdl->bond_price(current_time, current_time + disc_times[i], short_rate);
		}
		return out;
	};


	std::vector<T> payments;
	BaseRateModel<T>* mdl;
	Timeline<T> dates;
};


template<typename T>
class InterestRateSwap
{
	// IRS consists of a fixed rate bond and floating leg bond and payer type
	// It also contains a pointer to a model which it is valued on.
	// Note: poorly specific code. Needs to be generalized.
public:
	FixedRateBond<T> fixed_leg;
	FloatingRateBond<T> floating_leg;
	BaseRateModel<T>* mdl;
	T rate;
	bool payer;

	InterestRateSwap() {};
	InterestRateSwap(T& fixed_tenor, T& fixed_maturity, T& float_tenor,
		T& float_maturity, T& short_rate, BaseRateModel<T>* mdl,
		bool payer, T _rate_ = 1.0)
		: mdl(mdl), payer(payer), rate(_rate_)
	{
		rate = _rate_;
		FixedRateBond<T> fixd(rate, fixed_maturity, fixed_tenor, mdl);
		FloatingRateBond<T> floating(float_maturity, float_tenor, mdl, short_rate);
		if (abs(rate - 1.0) < 0.000001) {
			T z = 0.0;
			rate = floating.npv(z, short_rate) / fixd.npv(z, short_rate);
			fixd.fixed_rate = rate;
		}
		//rate = r;
		fixed_leg = fixd;
		floating_leg = floating;
	};
	~InterestRateSwap() {};

	T npv(T& current_time, T& short_rate) {
		// Calculate net present value at time "current_time" and at short rate "short_rate"
		T out = floating_leg.npv(current_time, short_rate) - fixed_leg.npv(current_time, short_rate);
		if (!payer) { out = -out; }
		
		return out;
	}
};
