#pragma once
template<typename T>
Matrix<T> IRS_bruteforce(InterestRateSwap<T>& swap, Matrix<T>& r,
	Timeline<T>& timeline) 
{
	// As the IRS has implemented a npv method which returns npv of the remaning cashflows
	// the brute force method is simple, yet still quite slow.
	size_t paths = r.nrows;
	size_t steps = timeline.size();
	Matrix<T> v(paths, steps);

	InterestRateSwap<T> swap_copy;
	Matrix<T> bank(paths, steps);
	for (size_t path = 0; path < paths; path++)
	{
		bank[path][0] = 1.0;
		for (size_t date = 1; date < steps; date++)
		{
			bank[path][date] = bank[path][date - 1] * exp(r[path][date] * (timeline[date]-timeline[date-1]));
		}
	}

	for (size_t path = 0; path < paths; path++)
	{
		swap_copy = swap;
		for (size_t step = 0; step < steps; step++)
		{
			v[path][step] = swap_copy.npv(timeline[step], r[path][step]) / bank[path][step];
		}
	}
	return v;
}


template<typename T>
Matrix<T> IRS_cash_flows(InterestRateSwap<T>& swap, Matrix<T>& r,
	Timeline<T>& timeline) 
{
	// Returns a matrix of realized cashflows for swap given interest rate paths r.
	// VERY Specific code. NEEDS to be generalized a lot to encompass swaps with different 
	// floating and fixed tenors
	size_t paths = r.nrows;
	size_t steps = timeline.size();
	Matrix<T> cf(paths, steps);
	for (size_t path = 0; path < paths; path++)
	{
		for (size_t step = 0; step < steps; step++)
		{
			cf[path][step] = 0.0;
		}
	}

	InterestRateSwap<T> swap_copy = swap;

	double type = -1.0;
	if (swap.payer) {
		type = 1.0;
	}

	for (size_t path = 0; path < paths; path++)
	{
		cf[path][5] = type * swap_copy.floating_leg.payments[0];
		for (size_t step = 1; step < swap_copy.floating_leg.dates.size(); step++)
		{
			cf[path][(step + 1) * 5] = type * swap.mdl->spot_libor(swap_copy.floating_leg.dates[step],
				swap_copy.floating_leg.dates[step] + 0.5, r[path][step * 5]) * 0.5;
		}
	}

	for (size_t path = 0; path < paths; path++)
	{
		for (size_t step = 0; step < swap_copy.fixed_leg.dates.size(); step++)
		{
			cf[path][(step + 1) * 10] = cf[path][(step + 1) * 10] -
				type * swap_copy.fixed_leg.fixed_rate;
		}
	}

	return cf;
}

template<typename T>
void IRS_cash_flows(InterestRateSwap<T>& swap, Matrix<T>& r,
	Timeline<T>& timeline, Matrix<T>& cf) {
	// Alters  the cf matrix to be the realized cashflows for swap given interest rate paths r.
	// VERY Specific code. NEEDS to be generalized a lot to encompass swaps with different 
	// floating and fixed tenors

	// Clear cf and adds swaps payments
	size_t paths = r.nrows;
	size_t steps = timeline.size();
	cf.resize(paths, steps);
	for (size_t path = 0; path < paths; path++)
	{
		for (size_t step = 0; step < steps; step++)
		{
			cf[path][step] = 0.0;
		}
	}

	InterestRateSwap<T> swap_copy = swap;

	double type = -1.0;
	if (swap.payer) {
		type = 1.0;
	}

	for (size_t path = 0; path < paths; path++)
	{
		cf[path][5] = type * swap_copy.floating_leg.payments[0];
		for (size_t step = 1; step < swap_copy.floating_leg.dates.size(); step++)
		{
			cf[path][(step + 1) * 5] = type * swap.mdl->spot_libor(swap_copy.floating_leg.dates[step],
				swap_copy.floating_leg.dates[step] + 0.5, r[path][step * 5]) * 0.5;
		}
	}

	for (size_t path = 0; path < paths; path++)
	{
		for (size_t step = 0; step < swap_copy.fixed_leg.dates.size(); step++)
		{
			cf[path][(step + 1) * 10] = (cf[path][(step + 1) * 10] -
				type * swap_copy.fixed_leg.fixed_rate);
		}
	}
}


template<typename T>
Matrix<T> IRS_proxy_beta_calculator(InterestRateSwap<T>& swap, ParentPhi<T>& phi,
	Matrix<T>& r, Matrix<T>& cash_flows,
	Timeline<T>& timeline) 
{
	// Returns a matrix containing all betas for the cashflow matrix.

	// Due to the specific nature of the phi class we also need to input the swap.
	// Also needs to be generalized for timelines with different increments

	size_t steps = timeline.size();
	size_t paths = cash_flows.nrows;

	Matrix<T> v(paths, steps);
	for (size_t path = 0; path < paths; path++)
	{
		v[path][steps - 1] = 0.0;
	}

	Matrix<T> beta(phi.get_basis_funcs() + 1, steps - 1);

	Matrix<T> phis;
	std::vector<T> X(paths);
	Matrix<T> Y(paths, 1);

	Matrix<T> V;
	std::vector<T> sig;
	Matrix<T> Sig;

	Matrix<T> tmp_beta;
	Matrix<T> svd_phis;
	
	
	Matrix<T> bank(paths, steps);
	for (size_t path = 0; path < paths; path++) 
	{		
		bank[path][0] = 1.0;	
		for (size_t date = 1; date < timeline.size(); date++)		
		{
			bank[path][date] = bank[path][date - 1] * exp(r[path][date] * 0.1);		
		}
	}
	
	for (int date = steps - 2; date > 0; date--)
	{
		T t = timeline[date];

		//Gather state variables 
		for (size_t path = 0; path < paths; path++)
		{
			X[path] = r[path][date];
			Y[path][0] = 0.0;
			for (size_t ii = date + 1; ii < timeline.size(); ii++)
			{
				Y[path][0] = Y[path][0] + (cash_flows[path][ii] / bank[path][ii]);
			}
			Y[path][0] = Y[path][0] * bank[path][date];
		}
		 
		//Perform svd regression
		phis = phi.mfunc(swap, X, t);
		svd_phis = phis;
		svd(svd_phis, sig, V);
		Sig = diag(sig, sig.size(), sig.size());
		diaginverse(Sig);
		tmp_beta = V * Sig * transpose(svd_phis) * Y; 

		for (size_t i = 0; i < tmp_beta.nrows; i++)
		{
			beta[i][date] = tmp_beta[i][0];
		}
	}
	return beta;
}

template<typename T>
Matrix<T> proxy_beta_calculator(InterestRateSwap<T>& swap, ParentPhi<T>& phi,
	Matrix<T>& r, Matrix<T>& cash_flows,
	Timeline<T>& timeline) {
	// Same function as IRS_proxy_beta_calculator. Not fully developed to not depend on swap


	size_t steps = timeline.size();
	size_t paths = cash_flows.nrows;

	Matrix<T> v(paths, steps);
	for (size_t path = 0; path < paths; path++)
	{
		v[path][steps - 1] = 0.0;
	}

	Matrix<T> beta(phi.get_basis_funcs() + 1, steps - 1);

	Matrix<T> phis;
	std::vector<T> X(paths);
	Matrix<T> Y(paths, 1);

	Matrix<T> V;
	std::vector<T> sig;
	Matrix<T> Sig;

	Matrix<T> tmp_beta;
	Matrix<T> svd_phis;


	Matrix<T> bank(paths, steps);
	for (size_t path = 0; path < paths; path++)
	{
		bank[path][0] = 1.0;
		for (size_t date = 1; date < timeline.size(); date++)
		{
			bank[path][date] = bank[path][date - 1] * exp(r[path][date] * 0.1);
		}
	}

	for (int date = steps - 2; date > 0; date--)
	{
		T t = timeline[date];

		//Gather state variables 
		for (size_t path = 0; path < paths; path++)
		{
			X[path] = r[path][date];
			Y[path][0] = 0.0;
			for (size_t ii = date + 1; ii < timeline.size(); ii++)
			{
				Y[path][0] = Y[path][0] + (cash_flows[path][ii] / bank[path][ii]);
			}
			Y[path][0] = Y[path][0] * bank[path][date];
		}

		//Perform svd regression
		phis = phi.mfunc(swap, X, t);
		svd_phis = phis;
		svd(svd_phis, sig, V);
		Sig = diag(sig, sig.size(), sig.size());
		diaginverse(Sig);
		tmp_beta = V * Sig * transpose(svd_phis) * Y; 

		for (size_t i = 0; i < tmp_beta.nrows; i++)
		{
			beta[i][date] = tmp_beta[i][0];
		}
	}
	return beta;
}


template<typename T, typename Y>
Matrix<T> IRS_proxy_value(InterestRateSwap<T>& swap, ParentPhi<T>& phi,
	Matrix<Y>& beta, Matrix<T>& r, Matrix<T>& cash_flows,
	Timeline<T>& timeline)
{
	// Returns the proxy value \tilde{V} given a Matrix of realized cash flows and beta matrix
	size_t steps = timeline.size();
	size_t paths = cash_flows.nrows;
	Matrix<T> v(paths, steps);
	for (size_t path = 0; path < paths; path++){ v[path][steps - 1] = 0.0; }

	Matrix<T> bank(paths, steps);
	for (size_t path = 0; path < paths; path++) {
		bank[path][0] = 1.0;
		for (size_t date = 1; date < timeline.size(); date++) { bank[path][date] = bank[path][date - 1] * exp(r[path][date] * 0.1); }
	}

	for (int date = steps - 2; date > 0; date--)
	{
		T t = timeline[date];
		//Gather state variables and update V
		for (size_t path = 0; path < paths; path++)
		{
			std::vector<T> phis = phi.func(swap, r[path][date], t);
			T tmp = 0.0;
			for (size_t j = 0; j < beta.nrows; j++)
			{
				tmp = tmp + phis[j] * beta[j][date];
			}
			v[path][date] = tmp / bank[path][date];
		}
	}

	T mean = 0.0;
	for (size_t path = 0; path < paths; path++)	{	mean = mean + v[path][1];	}
	mean = swap.mdl->bond_price(timeline[0], timeline[1], swap.mdl->r0) *  mean / paths;
	for (size_t path = 0; path < paths; path++)	{	v[path][0] = mean; }

	return v;
}


template<typename T>
std::vector<T> EPE_from_value_matrix(Matrix<T>& v, Timeline<T>& timeline) 
{
	// C++ alternative to the R-version of colMeans(pmax(v,0))
	std::vector<T> EPE(timeline.size());
	T timestep_value = 0.0;
	for (size_t time = 0; time < timeline.size(); time++)
	{
		timestep_value = 0.0;
		for (size_t path = 0; path < v.nrows; path++)
		{
			timestep_value = timestep_value + smoothed_indicator(v[path][time]) * v[path][time];
		}
		EPE[time] = timestep_value / static_cast<double>(v.nrows);
	}
	return EPE;
}

template<typename T>
T SD_CVA_from_value_matrix(Matrix<T>& v, Timeline<T>& timeline) 
{
	// Not fully developed. Needs to have Credit model added

	std::vector<T> CVA(v.nrows);

	for (size_t pp = 0; pp < v.nrows; pp++)
	{
		CVA[pp] = 0.0;
		for (size_t ii = 0; ii < v.ncols; ii++)
		{
			CVA[pp] += (v[pp][ii] > 0 ? v[pp][ii] : 0);
		}
		CVA[pp] *= 1e8 * (1 - 0.4) * 0.1 * 0.02; // 2% flat yearly PD and 0.1 timesteps. Poor code. 
	}

	T sum1 = 0.0;
	T sum2 = 0.0;
	for (size_t uu = 0; uu < v.nrows; uu++)
	{
		sum1 += CVA[uu] * CVA[uu];
		sum2 += CVA[uu];
	}
	sum1 = sum1 / double(v.nrows);
	sum2 = sum2 * sum2 / double(v.nrows*v.nrows);
	T SE = sqrt(sum1 - sum2);
	return SE;
}

template<typename T>
T smoothed_indicator(T& inp) {
	double eps = 1e-9;
	bool ind1 = inp > -(eps*0.5) ? true : false;
	bool ind2 = inp > (eps*0.5) ? true : false;

	if (ind1 && ind2) {
		return inp / eps - inp / eps + 1.0;
	}

	if (ind1) {
		return inp / eps + 0.5;
	}

	return 0.0;
}

template<typename T, typename Y>
T IRS_POI(InterestRateSwap<T>& swap, ParentPhi<T>& phi,
	Matrix<Y>& beta, Matrix<T>& r, Matrix<T>& cash_flows,
	Matrix<T>& proxy_v,
	Timeline<T>& timeline)
{
	// Applies POI method under the assumptions of flat 2% yearly default probability
	size_t steps = timeline.size();
	size_t paths = cash_flows.nrows;

	Matrix<T> POI(paths, steps);
	for (size_t path = 0; path < paths; path++)
	{
		POI[path][0] = 0.0;
	}

	T indicator;
	std::vector<T> phis;

	for (size_t date = 1; date < timeline.size(); date++)
	{
		T t = timeline[date];
		for (size_t path = 1; path < paths; path++)
		{
			POI[path][date] = POI[path][date - 1] +
				0.02 * 0.1 * smoothed_indicator(proxy_v[path][date]);
		}
	}

	T cva = 0.0;
	T path_bank;
	for (size_t path = 0; path < paths; path++)
	{
		path_bank = 1.0;
		for (size_t date = 1; date < timeline.size(); date++)
		{
			path_bank = path_bank * exp(r[path][date] * (timeline[date] - timeline[date - 1]));
			cva = cva + POI[path][date] * cash_flows[path][date] / path_bank;
		}
	}

	cva = cva / static_cast<double>(paths);

	return cva;
}

template<typename T, typename Y>
T SD_POI_Example(InterestRateSwap<T>& swap, ParentPhi<T>& phi,
	Matrix<Y>& beta, Matrix<T>& r, Matrix<T>& cash_flows,
	Matrix<T>& proxy_v,
	Timeline<T>& timeline)
{
	// Applies POI method under the assumptions of flat 2% yearly default probability on each
	// path and gives standard error
	size_t steps = timeline.size();
	size_t paths = cash_flows.nrows;

	Matrix<T> POI(paths, steps);
	for (size_t path = 0; path < paths; path++)
	{
		POI[path][0] = 0.0;
	}

	T indicator;
	std::vector<T> phis;

	for (size_t date = 1; date < timeline.size(); date++)
	{
		T t = timeline[date];
		for (size_t path = 1; path < paths; path++)
		{
			POI[path][date] = POI[path][date - 1] +
				 0.02 * 0.1 * smoothed_indicator(proxy_v[path][date]);
		}
	}

	T cva = 0.0;
	std::vector<T> path_cva(paths);
	T path_bank;
	for (size_t path = 0; path < paths; path++)
	{
		path_bank = 1.0;
		path_cva[path] = 0.0;
		for (size_t date = 1; date < timeline.size(); date++)
		{
			path_bank = path_bank * exp(r[path][date] * (timeline[date] - timeline[date - 1]));
			path_cva[path] = path_cva[path] + POI[path][date] * cash_flows[path][date] / path_bank;
		}
		path_cva[path] *= 1e8 * (1 - 0.4);
	}

	T sum1 = 0.0;
	T sum2 = 0.0;
	for (size_t uu = 0; uu < paths; uu++)
	{
		sum1 += path_cva[uu] * path_cva[uu];
		sum2 += path_cva[uu];
	}
	sum1 = sum1 / double(paths);
	sum2 = sum2 * sum2 / double(paths*paths);
	return sqrt(sum1-sum2);
}



template<typename T, typename Y>
T POI_from_value_matrix(InterestRateSwap<T>& swap, ParentPhi<T>& phi,
	Matrix<Y>& beta, Matrix<T>& r, Matrix<T>& cash_flows,
	Matrix<T>& proxy_v,	Timeline<T>& timeline, creditObject<T>& creditObj)
{
	// OBS: Does not ues phi or beta as it already has proxy V. Poor code.

	// Actual POI method which uses the cox credit model. Returns the CVA value given
	// realized cashflows, realized rates and the proxy values
	size_t steps = timeline.size();
	size_t paths = cash_flows.nrows;

	Matrix<T> POI(paths, steps);
	for (size_t path = 0; path < paths; path++)
	{
		POI[path][0] = 0.0;
	}

	T indicator;
	std::vector<T> phis;

	for (size_t date = 1; date < timeline.size(); date++)
	{
		T t = timeline[date];
		T default_prob = creditObj.default_prob(timeline[date - 1], timeline[date]);

		for (size_t path = 0; path < paths; path++)
		{
			POI[path][date] = POI[path][date - 1] + default_prob * smoothed_indicator(proxy_v[path][date]);
		}
	}

	T cva = 0.0;
	T path_bank;

	for (size_t path = 0; path < paths; path++)
	{
		path_bank = 1.0;
		for (size_t date = 1; date < timeline.size(); date++)
		{
			path_bank = path_bank * exp(r[path][date] * (timeline[date] - timeline[date - 1]));
			cva = cva + POI[path][date] * cash_flows[path][date] / path_bank;
		}
	}
	cva = cva / static_cast<double>(paths);

	return cva;
}
