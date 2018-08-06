#pragma once

template<typename T>
class BaseRateModel
{
	// Interface class for the short rate models
	// All short rate models are set up to simulate future scenarios
public:
	BaseRateModel() {};
	virtual T forward_rate(T start, T maturity, T& short_rate) = 0;
	virtual T disc_factor(T start, T end, T& short_rate) = 0;
	virtual T bond_price(T start, T maturity, T& short_rate) = 0;
	virtual Matrix<T> simulate(Timeline<T>& timeline, const size_t paths, const T& r0, mrg32& rng, bool byrow = false) = 0;
	virtual T spot_libor(T current_time, T maturity, T& short_rate) = 0;
	T r0;
};
