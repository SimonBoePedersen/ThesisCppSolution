#pragma once
// This header defines Basis function classes needed for LSM and proxy generation.
// Currently this code is too specific and has to be generalized more (Currently can only take in 
//  a swap, needs to be generalized to loads of other instruments).
// Four classes are defined. We only use the PowerShortRatePhi and PowerParRatePhi


template<typename T>
class ParentPhi
{
	// Interface class for basis functions needed in LSM
public:
	ParentPhi() {};
	~ParentPhi() {};
	virtual std::vector<T> func(InterestRateSwap<T>& swap, T short_rate, T current_time) = 0;
	virtual Matrix<T> mfunc(InterestRateSwap<T>& swap, std::vector<T> short_rate, T current_time) = 0;
	virtual size_t get_basis_funcs() = 0;
};

template<typename T>
class PowerShortRatePhi : public ParentPhi<T>
{
	// Basis function class which uses returns r^0, r^1, r^2, etc. with r = short rate
public:
	PowerShortRatePhi(size_t basis_funcs) : basis_funcs(basis_funcs) {};
	~PowerShortRatePhi() {};

	virtual std::vector<T> func(InterestRateSwap<T>& swap, T short_rate,
		T current_time) {
		std::vector<T> out(basis_funcs + 1);
		out[0] = 1;
		for (size_t i = 1; i < basis_funcs + 1; i++) { out[i] = out[i - 1] * short_rate; }
		return out;
	}

	virtual Matrix<T> mfunc(InterestRateSwap<T>& swap, std::vector<T> short_rate,
		T current_time) {

		Matrix<T> out(short_rate.size(), basis_funcs + 1);
		for (size_t path = 0; path < short_rate.size(); path++)
		{
			out[path][0] = 1;
			for (size_t i = 1; i < basis_funcs + 1; i++) {
				out[path][i] = out[path][i - 1] * short_rate[path];
			}
		}
		return out;
	}

	virtual size_t get_basis_funcs() {
		return basis_funcs;
	}

	size_t basis_funcs;
};

template<typename T>
class PowerParRatePhi : public ParentPhi<T>
{
	// Basis function class which uses returns r^0, r^1, r^2, etc. with r = par rate
public:
	PowerParRatePhi(size_t basis_funcs) : basis_funcs(basis_funcs) {};
	~PowerParRatePhi() {};

	virtual std::vector<T> func(InterestRateSwap<T>& swap, T short_rate,
		T current_time) {
		std::vector<T> out(basis_funcs + 1);
		std::vector<T> safekeeping = swap.floating_leg.payments;
		InterestRateSwap<T> copy = swap;

		for (size_t i = 0; i < copy.fixed_leg.dates.size(); i++) {
			copy.fixed_leg.dates[i] = swap.fixed_leg.dates[i] + current_time;
		}
		for (size_t i = 0; i < copy.floating_leg.dates.size(); i++) {
			copy.floating_leg.dates[i] = swap.floating_leg.dates[i] + current_time;
		}

		copy.floating_leg.payments[0] = (copy.floating_leg.dates[0] - current_time) *
			swap.mdl->spot_libor(current_time, copy.floating_leg.dates[0], short_rate);
		T floating_npv = copy.floating_leg.npv(current_time, short_rate);
		T fixed_npv = copy.fixed_leg.npv(current_time, short_rate) / swap.fixed_leg.fixed_rate;
		T par_rate = floating_npv / fixed_npv;

		out[0] = 1.0;
		for (size_t i = 1; i < basis_funcs + 1; i++) {

			out[i] = out[i - 1] * par_rate;
		}

		swap.floating_leg.payments = safekeeping;
		return out;
	}

	virtual Matrix<T> mfunc(InterestRateSwap<T>& swap, std::vector<T> short_rate,
		T current_time) {

		Matrix<T> out(short_rate.size(), basis_funcs + 1);

		InterestRateSwap<T> copy = swap;
		std::vector<T> safekeeping = swap.floating_leg.payments;
		std::vector<T> timelineFixed = swap.fixed_leg.dates.line;
		std::vector<T> timelineFloat = swap.floating_leg.dates.line;

		for (size_t i = 0; i < copy.fixed_leg.dates.size(); i++) {
			copy.fixed_leg.dates[i] = swap.fixed_leg.dates[i] + current_time;}
		for (size_t i = 0; i < copy.floating_leg.dates.size(); i++) {
			copy.floating_leg.dates[i] = swap.floating_leg.dates[i] + current_time;}


		for (size_t path = 0; path < short_rate.size(); path++)
		{
			out[path][0] = 1;
			
			copy.floating_leg.payments[0] = (copy.floating_leg.dates[0]-current_time) * 
				swap.mdl->spot_libor(current_time, copy.floating_leg.dates[0], short_rate[path]);
			T floating_npv = copy.floating_leg.npv(current_time, short_rate[path]);
			T fixed_npv = copy.fixed_leg.npv(current_time, short_rate[path]) / swap.fixed_leg.fixed_rate;
			T par_rate = floating_npv / fixed_npv;
			for (size_t i = 1; i < basis_funcs + 1; i++) {
				out[path][i] = out[path][i - 1] * par_rate;
			}
		}

		swap.floating_leg.payments = safekeeping;
		return out;
	}

	virtual size_t get_basis_funcs() {
		return basis_funcs;
	}

	size_t basis_funcs;
};


template<typename T>
class Power6MLiborPhi : public ParentPhi<T>
{
public:
	Power6MLiborPhi(size_t basis_funcs, T tenor) : tenor(tenor), basis_funcs(basis_funcs) {};
	~Power6MLiborPhi() {};

	virtual std::vector<T> func(InterestRateSwap<T>& swap, T short_rate,
		T current_time) {
		std::vector<T> out(basis_funcs + 1);
		out[0] = 1;
		for (size_t i = 1; i < basis_funcs + 1; i++) { out[i] = out[i - 1] * 
			swap.mdl->spot_libor(current_time,current_time +  tenor, short_rate); }
		return out;
	}

	virtual Matrix<T> mfunc(InterestRateSwap<T>& swap, std::vector<T> short_rate,
		T current_time) {

		Matrix<T> out(short_rate.size(), basis_funcs + 1);
		for (size_t path = 0; path < short_rate.size(); path++)
		{
			out[path][0] = 1;
			T sixmlibor = swap.mdl->spot_libor(current_time, current_time + tenor, short_rate[path]);
			for (size_t i = 1; i < basis_funcs + 1; i++) {
				out[path][i] = out[path][i - 1] * sixmlibor;
			}
		}
		return out;
	}

	virtual size_t get_basis_funcs() {
		return basis_funcs;
	}
	T tenor;
	size_t basis_funcs;
};

template<typename T>
class PowerDoubleLiborPhi : public ParentPhi<T>
{
public:
	PowerDoubleLiborPhi(size_t basis_funcs, T tenor) : tenor(tenor), basis_funcs(basis_funcs) {};
	~PowerDoubleLiborPhi() {};

	virtual std::vector<T> func(InterestRateSwap<T>& swap, T short_rate,
		T current_time) {
		std::vector<T> out(basis_funcs + 1);
		out[0] = 1;
		out[1] = swap.mdl->spot_libor(current_time, current_time + tenor, short_rate);
		out[2] = swap.mdl->spot_libor(current_time, current_time + 2 * tenor, short_rate);
		return out;
	}

	virtual Matrix<T> mfunc(InterestRateSwap<T>& swap, std::vector<T> short_rate,
		T current_time) {

		Matrix<T> out(short_rate.size(), basis_funcs + 1);
		for (size_t path = 0; path < short_rate.size(); path++)
		{
			out[path][0] = 1;
			out[path][1] = swap.mdl->spot_libor(current_time, current_time + tenor, short_rate[path]);
			out[path][2] = swap.mdl->spot_libor(current_time, current_time + 2 * tenor, short_rate[path]);
		}
		return out;
	}

	virtual size_t get_basis_funcs() {
		return basis_funcs;
	}
	T tenor;
	size_t basis_funcs;
};



template<typename T>
PowerParRatePhi<double> as_double(PowerParRatePhi<T>& phi) {
	PowerParRatePhi<double> out(phi.basis_funcs);
	return out;
}

template<typename T>
PowerShortRatePhi<double> as_double(PowerShortRatePhi<T>& phi) {
	PowerShortRatePhi<double> out(phi.basis_funcs);
	return out;
}

template<typename T>
Power6MLiborPhi<double> as_double(Power6MLiborPhi<T>& phi) {
	Power6MLiborPhi<double> out(phi.basis_funcs);
	return out;
}

template<typename T>
PowerDoubleLiborPhi<double> as_double(PowerDoubleLiborPhi<T>& phi) {
	PowerDoubleLiborPhi<double> out(phi.basis_funcs);
	return out;
}

