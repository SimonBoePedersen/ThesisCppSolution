#pragma once

template<typename T>
void matrixwriter2(std::string location, Matrix<T>& matrix) {
	std::ofstream myFile;
	myFile.open(location);

	for (size_t path = 0; path < matrix.nrows; path++)
	{
		for (size_t date = 0; date < matrix.ncols; date++)
		{
			myFile << value(matrix[path][date]) << ",";
		}
		myFile << "\n";
	}

	myFile.flush();
	myFile.close();

}


void hullWhiteFittingExample() {
	// Example used after Hull-White section
	// Takes the curve, fits it and prints all model prices and new model parameters
	std::vector<double> mats = { 0.1, 3.0 / 12.0, 0.5,1,2,3,	5,8,10,12,15,20,30 };
	std::vector<double> rats = { 0.035877, 0.036576, 0.037741,0.038560, 0.038484,
		0.038218,0.038186, 0.038787,0.039264, 0.039695,0.040222,0.040828,0.041475 };
	HullWhiteModel<double>* mdl = new HullWhiteModel<double>(mats, rats);

	Matrix<double> mdl_prices(mats.size(), 1);
	for (size_t i = 0; i < mats.size(); i++)
	{
		mdl_prices[i][0] = mdl->bond_price(0, mats[i], mdl->r0);
		std::cout << mats[i] << " Model prices " << mdl_prices[i][0] << "  Market prices  " << exp(-rats[i] * mats[i]) << std::endl;
	}
	std::cout << "Model parameters: a = " << mdl->a << ", sigma = " << mdl->sigma << std::endl;
}


void example_after_proxy() {
	// This is used to generate the values needed in the example after the Pure Proxy section
	// Varying paths both pre and main are used to generate the numbers in the thesis
	// We also write the results to disk to use third-party plotting software

	VasicekModel<double> *mdl = new VasicekModel<double>(0.04, 0.1, 0.015, 0.04);

	size_t pre_paths = 2 * 256;
	size_t main_paths = 4 * 256;
	double maturity = 10.0;

	double fixed_rate = 1.0;
	double fixed_tenor = 1;
	double float_tenor = 0.5;

	Timeline<double> timeline(0.0, maturity, 0.1);
	InterestRateSwap<double> swap(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, fixed_rate);
	PowerParRatePhi<double> phi(1);
	

	mrg32 rng(123, 1234, 1234, 12345, 12345, 12346);
	Matrix<double> r_pre = mdl->simulate(timeline, pre_paths, mdl->r0, rng);
	Matrix<double> cf_pre = IRS_cash_flows(swap, r_pre, timeline);
	Matrix<double> beta = IRS_proxy_beta_calculator(swap, phi, r_pre, cf_pre, timeline);

	// Mainsimulation
	rng.setseeds(1, 2, 3, 4, 5, 6);
	Matrix<double> r = mdl->simulate(timeline, main_paths, mdl->r0, rng);
	Matrix<double> cash_flows = IRS_cash_flows(swap, r, timeline);

	Matrix<double> pure_proxy_value = IRS_proxy_value(swap, phi,
		beta, r, cash_flows, timeline);									  //Pure proxy
	Matrix<double> brute_force_value = IRS_bruteforce(swap, r, timeline); //Brute

	// OVERSKRIV ALDRIG DISSE 
	// Laptop:
	//matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_pure_proxy/pure_proxy256.txt", pure_proxy_value);
	//matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_pure_proxy/brute_force256.txt", brute_force_value);
	// Desktop:
	//matrixwriter("C:/Users/Simon/Dropbox/Uni/Speciale/TeX/document/examples/example_after_pure_proxy/pure_proxy512.txt", pure_proxy_value);
	//matrixwriter("C:/Users/Simon/Dropbox/Uni/Speciale/TeX/document/examples/example_after_pure_proxy/brute_force512.txt", brute_force_value);


	std::vector<double> EPE_proxy = EPE_from_value_matrix(pure_proxy_value, timeline);
	std::vector<double> EPE_brute = EPE_from_value_matrix(brute_force_value, timeline);
	
	double cva_proxy = 0.0;
	double cva_brute = 0.0;
	for (size_t i = 0; i < timeline.size(); i++)
	{
		cva_proxy += EPE_proxy[i] * 0.02 * 0.1;
		cva_brute += EPE_brute[i] * 0.02 * 0.1;
	}
	std::cout << "Fixed Rate " << swap.fixed_leg.fixed_rate << std::endl;
	std::cout << "CVA of Proxy: " << 1e8 * cva_proxy * (1 - 0.4) << " SE: " 
		<< SD_CVA_from_value_matrix(pure_proxy_value, timeline) / sqrt(double(main_paths)) << std::endl;
	std::cout << "CVA of Brute: " << 1e8 * cva_brute * (1 - 0.4) << " SE: " 
		<< SD_CVA_from_value_matrix(brute_force_value, timeline) / sqrt(double(main_paths)) << std::endl;
	delete mdl;
}


void example_after_poi() {
	// This is used to generate the values needed in the example after the POI section
	// Varying paths both pre and main are used to generate the numbers in the thesis
	// Note that this ALSO includes the calculation of expected exposures with POI - this is ONLY for plotting purposes
	VasicekModel<double> *mdl = new VasicekModel<double>(0.04, 0.1, 0.015, 0.04);

	size_t pre_paths = 2 * 256;
	size_t main_paths = 4 * 256;
	double maturity = 10.0;

	double fixed_rate = 1.0;
	double fixed_tenor = 1;
	double float_tenor = 0.5;

	Timeline<double> timeline(0.0, maturity, 0.1);
	InterestRateSwap<double> swap(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, fixed_rate);
	PowerParRatePhi<double> phi(1);

	mrg32 rng(123, 1234, 1234, 12345, 12345, 12346);
	Matrix<double> r_pre = mdl->simulate(timeline, pre_paths, mdl->r0, rng);
	Matrix<double> cf_pre = IRS_cash_flows(swap, r_pre, timeline);
	Matrix<double> beta = IRS_proxy_beta_calculator(swap, phi, r_pre, cf_pre, timeline);

	// Mainsimulation
	rng.setseeds(1, 2, 3, 4, 5, 6);
	Matrix<double> r = mdl->simulate(timeline, main_paths, mdl->r0, rng);
	Matrix<double> cash_flows = IRS_cash_flows(swap, r, timeline);

	size_t steps = timeline.size();
	size_t paths = cash_flows.nrows;

	Matrix<double> INDICATOR(paths, steps);
	Matrix<double> POI(paths, steps);
	Matrix<double> value(paths, steps);
	for (size_t path = 0; path < paths; path++)
	{
		value[path][steps - 1] = 0.0;
		POI[path][0] = 0.0;
	}

	// Reference value matrices
	Matrix<double> pure_proxy_value = IRS_proxy_value(swap, phi,
		beta, r, cash_flows, timeline);
	Matrix<double> brute_force_value = IRS_bruteforce(swap, r, timeline);
	
	Matrix<double> bank(paths, steps);
	for (size_t path = 0; path < paths; path++) {
		bank[path][0] = 1.0;
		for (size_t date = 1; date < timeline.size(); date++) { bank[path][date] = bank[path][date - 1] * exp(r[path][date] * 0.1); }
	}

	double indicator;
	std::vector<double> phis;
	double t;
	for (size_t date = 1; date < timeline.size() - 1; date++)
	{
		t = timeline[date];
		for (size_t path = 0; path < paths; path++)
		{
			// Calculate indicator
			phis = phi.func(swap, r[path][date], t);
			indicator = 0.0;
			for (size_t j = 0; j < beta.nrows; j++)
			{
				indicator = indicator + beta[j][date] * phis[j];
			}
			indicator = indicator * bank[path][date];
			INDICATOR[path][date] = smoothed_indicator(indicator);
		}
	}

	for (size_t path = 0; path < paths; path++)
	{
		for (size_t date = 1; date < timeline.size() - 1; date++)
		{
			if (pure_proxy_value[path][date] > 0.0) {
				t = timeline[date];
				POI[path][date] = POI[path][date - 1] + 0.02 * 0.1;
				double v = 0.0;
				for (size_t innersum = date + 1; innersum < cash_flows.ncols; innersum++)
				{	
					v = v + cash_flows[path][innersum] / bank[path][innersum];
				}
				value[path][date] = v;
			}
			else
			{
				value[path][date] = 0;
				POI[path][date] = POI[path][date - 1];
 			}
		}
	}

	double v0 = 0.0;
	double disc = swap.mdl->bond_price(0, timeline[1], mdl->r0);
	for (size_t path = 0; path < paths; path++)
	{
		v0 += value[path][1];
	}
	for (size_t path = 0; path < paths; path++)
	{
		value[path][0] = v0 * disc / paths;
	}

	std::cout << "CVA of POI:   " << IRS_POI(swap, phi, beta, r, cash_flows, pure_proxy_value, timeline) * 1e8 * (1 - 0.4) << std::endl;
	std::cout << "Standard Error " << SD_POI_Example(swap, phi, beta, r, cash_flows, pure_proxy_value, timeline) / double(sqrt(main_paths)) << std::endl;


	// Desktop:
	//matrixwriter("C:/Users/Simon/Dropbox/Uni/Speciale/TeX/document/examples/example_after_poi/pure_proxy16_4096.txt", pure_proxy_value);
	//matrixwriter("C:/Users/Simon/Dropbox/Uni/Speciale/TeX/document/examples/example_after_poi/brute_force16_4096.txt", brute_force_value);
	//matrixwriter("C:/Users/Simon/Dropbox/Uni/Speciale/TeX/document/examples/example_after_poi/poi16_4096.txt", value);

	/*
	matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_poi/poi512.txt", value);
	matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_poi/pure_proxy512.txt", pure_proxy_value);
	matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_poi/brute_force512.txt", brute_force_value);
	*/
	/*
	matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_poi/poi16.txt", value);
	matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_poi/pure_proxy16.txt", pure_proxy_value);
	matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_poi/brute_force16.txt", brute_force_value);
	*/
/*
	matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_poi/poi16_4096.txt", value);
	matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_poi/pure_proxy16_4096.txt", pure_proxy_value);
	matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_poi/brute_force16_4096.txt", brute_force_value);
*/
	//matrixwriter("C:/Users/Simon/Documents/Dropbox/Uni/Speciale/TeX/structured_models/examples/example_after_poi/poi_function.txt", POI);
	delete mdl;
}


template<typename T>
T swap_formula(InterestRateSwap<T>& swp, std::vector<T>& mats, std::vector<T>& rates) {
	HullWhiteModel<T> *mdl = new HullWhiteModel<T>(mats, rates);
	// Used for example in AAD section
	T sumsum = 0.0;
	for (size_t i = 0; i < swp.fixed_leg.dates.size(); i++)
		sumsum = sumsum + 1 * swp.rate * exp(-swp.fixed_leg.dates[i] * mdl->zeroCurve_obj.interpolate(swp.fixed_leg.dates[i]));
	T npvformula = (1 - exp(- 10 * rates[9]) - sumsum); // Poor hardcoded 10Y swap formula
	
	delete mdl;
	return npvformula;
}

void example_aad_swap_delta_vector() {
	// Produces the delta vector of a 10Y par payer swap in Hull-White model with
	// finite differencing and AAD.
	// AAD instrumented part.

	std::vector<Variable> mats = { 0.0,	0.1, 3.0 / 12.0, 0.5,1,2,3,5,8,10,12,15,20,30 };
	std::vector<Variable> rats = { 0.035877,0.035877,  0.036576, 0.037741,
		0.038560, 0.038484,	0.038218, 0.038186,	0.038787,
		0.039264, 0.039695,	0.040222, 0.040828, 0.041475 };

	HullWhiteModel<Variable> *mdl = new HullWhiteModel<Variable>(mats, rats);

	Variable fixed_rate = 0.0398985;
	Variable maturity = 10;
	Variable fixed_tenor = 1;
	Variable float_tenor = 0.5;
	InterestRateSwap<Variable> swp(fixed_tenor, maturity, float_tenor, maturity,
		mdl->r0, mdl, true, fixed_rate);
	Variable res = swap_formula(swp, mats, rats);
	std::vector<double> derivs = res.calc_derivatives(); // Backwards sweep, stores all adjoints in res
	GlobalTape->nodes.clear();
	delete mdl;

	// Finite difference part
	std::vector<double> mats_d = { 0.0,	0.1, 3.0 / 12.0, 0.5,1,2,3,
		5,8,10,12,15,20,30 };
	std::vector<double> rats_d = { 0.035877,0.035877,  0.036576, 0.037741,
		0.038560, 0.038484,	0.038218, 0.038186,	0.038787,
		0.039264, 0.039695,	0.040222, 0.040828, 0.041475 };

	HullWhiteModel<double> *mdl_d = new HullWhiteModel<double>(mats_d, rats_d);

	double fixed_rate_d = 0.0398985;
	double maturity_d = 10;
	double fixed_tenor_d = 1;
	double float_tenor_d = 0.5;
	InterestRateSwap<double> swp_d(fixed_tenor_d, maturity_d, float_tenor_d, maturity_d,
		mdl_d->r0, mdl_d, true, fixed_rate_d);

	double eps = 1e-8;
	double multiplier = 1e8 * 1e-4; // Notional 100M and we want bp values.
	std::cout << "Maturity / AAD / Forward finite difference" << std::endl;
	for (size_t i = 1; i < mats_d.size(); i++)
	{
		double unbumped = swap_formula(swp_d, mats_d, rats_d);
		// Bump
		rats_d[i] += eps;
		double bumped = swap_formula(swp_d, mats_d, rats_d);
		std::cout << mats_d[i] << " / " <<
			multiplier * wrt(rats[i], derivs) << " / " << 
			multiplier * (bumped - unbumped) / eps << std::endl;

		rats_d[i] -= eps;

	}
	delete mdl_d;
}

void CVA_delta_bumps(double idx, double eps) {

	std::vector<double> mats = {0.1, 3.0 / 12.0, 0.5,1,2,3,5,8,10,12,15,20,30 };
	std::vector<double> rats = {0.035877, 0.036576, 0.037741,0.038560, 0.038484,0.038218,0.038186, 0.038787,
		0.039264, 0.039695,0.040222,0.040828,0.041475 };
	rats[idx] -= eps;

	HullWhiteModel<double> *mdl = new HullWhiteModel<double>(mats, rats);

	double fixed_rate = 0.04;
	double maturity = 10;
	double fixed_tenor = 1;
	double float_tenor = 0.5;
	InterestRateSwap<double> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, fixed_rate);

	size_t pre_paths = 256;
	size_t main_paths = 1024;

	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);
	Timeline<double> timeline(0.0, maturity, 0.1);
	Matrix<double> r = mdl->simulate(timeline, pre_paths, mdl->r0, rng);
	Matrix<double> cf = IRS_cash_flows(swp, r, timeline);

	PowerShortRatePhi<double> basis(1);
	Matrix<double> beta_proxy = IRS_proxy_beta_calculator(swp, basis, r, cf, timeline);
	rng.setseeds(1, 2, 3, 4, 5, 6);
	Matrix<double> r2 = mdl->simulate(timeline, main_paths, mdl->r0, rng, true);
	Matrix<double> cf2 = IRS_cash_flows(swp, r2, timeline);
	Matrix<double> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);

	double cva_poi1 =  IRS_POI(swp, basis, beta_proxy, r2, cf2, v2, timeline);

	rats[idx] += eps;
	std::vector<double> rats2 = rats;
	rats2[idx] += eps;
	
	HullWhiteModel<double> *mdl2 = new HullWhiteModel<double>(mats, rats2);
	InterestRateSwap<double> swp2(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl2, true, fixed_rate);

	rng.setseeds(1234, 1234, 1234, 12345, 12346, 12346);
	r = mdl2->simulate(timeline, pre_paths, mdl2->r0, rng);
	cf = IRS_cash_flows(swp2, r, timeline);

	beta_proxy = IRS_proxy_beta_calculator(swp2, basis, r, cf, timeline);
	rng.setseeds(1, 2, 3, 4, 5, 6);
	r2 = mdl2->simulate(timeline, main_paths, mdl2->r0, rng, true);
	cf2 = IRS_cash_flows(swp2, r2, timeline);
	v2 = IRS_proxy_value(swp2, basis, beta_proxy, r2, cf2, timeline);

	double cva_poi2 =  IRS_POI(swp2, basis, beta_proxy, r2, cf2, v2, timeline);
	std::cout <<  "   "  << (100000000.0 / 10000.0) * (cva_poi2 - cva_poi1)/(2*eps) << std::endl;
	delete mdl, mdl2;
}


void CVA_delta_vector() {

	long before = GetTickCount();
	std::vector<Variable> mats = { 
		//0.0,
		0.1, 3.0 / 12.0, 0.5,
		1,2,3,
		5,8,
		10,12,
		15,
		20,
		30 };
	std::vector<Variable> rats = { 
		//0.035877,
		0.035877,  0.036576,	 0.037741,
		0.038560,	 0.038484,		0.038218,
		0.038186,	 0.038787,
		0.039264,	  0.039695,
		0.040222,
		0.040828,
		0.041475 };
	HullWhiteModel<Variable> *mdl = new HullWhiteModel<Variable>(mats, rats);
	

	Variable fixed_rate = 0.03;
	Variable maturity = 10;
	Variable fixed_tenor = 1;
	Variable float_tenor = 0.5;
	InterestRateSwap<Variable> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, fixed_rate);

	size_t pre_paths = 256;
	size_t main_paths = 1024;
	
	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);
	Timeline<Variable> timeline(0.0, maturity, 0.1);
	Matrix<Variable> r = mdl->simulate(timeline, pre_paths, mdl->r0, rng);
	Matrix<Variable> cf = IRS_cash_flows(swp, r, timeline);

	PowerShortRatePhi<Variable> basis(1);
	Matrix<Variable> beta_proxy = IRS_proxy_beta_calculator(swp, basis, r, cf, timeline);
	rng.setseeds(1, 2, 3, 4, 5, 6);
	Matrix<Variable> r2 = mdl->simulate(timeline, main_paths, mdl->r0, rng);
	Matrix<Variable> cf2 = IRS_cash_flows(swp, r2, timeline);
	Matrix<Variable> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);

	Variable cva_poi = IRS_POI(swp, basis, beta_proxy, r2, cf2, v2, timeline);
	std::cout << "CVA POI: " << 100000000 * cva_poi.value << std::endl;
	std::vector<double> derivs = cva_poi.calc_derivatives();
	
	long after = GetTickCount();
	std::cout << after - before << "ms in aad" << std::endl;
	
	before = GetTickCount();
	for (size_t i = 0; i < mats.size(); i++)
	{
		std::cout << as_double(mats[i]) << " " << 100000000 * wrt(rats[i], derivs) / 10000;
		CVA_delta_bumps(i, 0.0000001);
	}
	after = GetTickCount();
	std::cout << after - before << "ms in fd" << std::endl;
	delete mdl;
}


void CVA_delta_vector_no_diff_presim() {
	std::vector<Variable> mats = {
		//0.0,
		0.1, 3.0 / 12.0, 0.5,
		1,2,3,
		5,8,
		10,12,
		15,
		20,
		30 };
	std::vector<Variable> rats = {
		//0.035877,
		0.035877,  0.036576,	 0.037741,
		0.038560,	 0.038484,		0.038218,
		0.038186,	 0.038787,
		0.039264,	  0.039695,
		0.040222,
		0.040828,
		0.041475 };

	//Variable a = 0.164883; Variable sigma = 0.0280143;
	HullWhiteModel<Variable> *mdl = new HullWhiteModel<Variable>(mats, rats);


	Variable fixed_rate = 0.03;
	Variable maturity = 10;
	Variable fixed_tenor = 1;
	Variable float_tenor = 0.5;
	InterestRateSwap<Variable> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, fixed_rate);

	size_t pre_paths = 256;
	size_t main_paths = 1024;

	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);
	Timeline<Variable> timeline(0.0, maturity, 0.1);
	PowerShortRatePhi<Variable> basis(1);
	
	// Proxy generation
	double maturity_double = as_double(maturity);
	std::vector<double> mats_double = as_double(mats);
	std::vector<double> rats_double = as_double(rats);
	HullWhiteModel<double> *mdl_double = new HullWhiteModel<double>(mats_double, rats_double);
	Timeline<double> timeline_double = as_double(timeline);
	InterestRateSwap<double> swp_double(fixed_tenor.value, maturity_double, float_tenor.value, maturity_double, mdl_double->r0, mdl_double, true, 1.0);

	Matrix<double> r = mdl_double->simulate(timeline_double, pre_paths, mdl_double->r0, rng);
	Matrix<double> cf = IRS_cash_flows(swp_double, r, timeline_double);
	PowerShortRatePhi<double> basis_double(basis.basis_funcs);
	Matrix<double> beta_proxy = IRS_proxy_beta_calculator(swp_double, basis_double, r, cf, timeline_double);

	// Main simulation
	rng.setseeds(1, 2, 3, 4, 5, 6);
	Matrix<Variable> r2 = mdl->simulate(timeline, main_paths, mdl->r0, rng);
	Matrix<Variable> cf2 = IRS_cash_flows(swp, r2, timeline);
	Matrix<Variable> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);

	Variable cva_poi = IRS_POI(swp, basis, beta_proxy, r2, cf2, v2, timeline);
	std::cout << "CVA POI: " << 100000000.0 * cva_poi.value << std::endl;
	std::vector<double> derivs = cva_poi.calc_derivatives();
	for (size_t i = 0; i < mats.size(); i++)
	{
		std::cout << as_double(mats[i]) << " " << (100000000.0 / 10000.0) * wrt(rats[i], derivs);
		CVA_delta_bumps(i, 0.0000001);
	}

	delete mdl;
}

void CVA_no_diff_pre_tape_per_path() {

	std::vector<Variable> mats = {
		//0.0,
		0.1, 3.0 / 12.0, 0.5,
		1,2,3,
		5,8,
		10,12,
		15,
		20,
		30 };
	std::vector<Variable> rats = {
		//0.035877,
		0.035877,  0.036576,	 0.037741,
		0.038560,	 0.038484,		0.038218,
		0.038186,	 0.038787,
		0.039264,	  0.039695,
		0.040222,
		0.040828,
		0.041475 };

	//Variable a = 0.164883; Variable sigma = 0.0280143;
	HullWhiteModel<Variable> *mdl = new HullWhiteModel<Variable>(mats, rats);


	Variable fixed_rate = 0.04;
	Variable maturity = 10;
	Variable fixed_tenor = 1;
	Variable float_tenor = 0.5;
	InterestRateSwap<Variable> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, fixed_rate);

	size_t pre_paths = 256;
	size_t main_paths = 1024;

	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);
	Timeline<Variable> timeline(0.0, maturity, 0.1);
	PowerShortRatePhi<Variable> basis(1);

	// Proxy generation
	double maturity_double = as_double(maturity);
	std::vector<double> mats_double = as_double(mats);
	std::vector<double> rats_double = as_double(rats);
	HullWhiteModel<double> *mdl_double = new HullWhiteModel<double>(mats_double, rats_double);
	Timeline<double> timeline_double = as_double(timeline);
	InterestRateSwap<double> swp_double(fixed_tenor.value, maturity_double, float_tenor.value, maturity_double, mdl_double->r0, mdl_double, true, 1.0);

	Matrix<double> r = mdl_double->simulate(timeline_double, pre_paths, mdl_double->r0, rng);
	Matrix<double> cf = IRS_cash_flows(swp_double, r, timeline_double);
	PowerShortRatePhi<double> basis_double(basis.basis_funcs);
	Matrix<double> beta_proxy = IRS_proxy_beta_calculator(swp_double, basis_double, r, cf, timeline_double);


	std::vector<Node> copy = GlobalTape->nodes;
	double cva_poi = 0.0;
	std::vector<double> path_derivs(rats.size());
	for (size_t i = 0; i < rats.size(); i++)
	{
		path_derivs[i] = 0.0;
	}

	
	// Main simulation
	rng.setseeds(1, 2, 3, 4, 5, 6);
	for (size_t i = 0; i < main_paths; i++)
	{
		Matrix<Variable> r2 = mdl->simulate(timeline, 1, mdl->r0, rng, true);
		Matrix<Variable> cf2 = IRS_cash_flows(swp, r2, timeline);
		Matrix<Variable> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);

		Variable cva = IRS_POI(swp, basis, beta_proxy, r2, cf2, v2, timeline);
		
		std::vector<double> derivs = cva.calc_derivatives();
		for (size_t i = 0; i < rats.size(); i++)
		{
			path_derivs[i] += wrt(rats[i], derivs);
		}

		cva_poi += cva.value;
		GlobalTape->nodes = copy;
	}
	std::cout <<  1e8 * cva_poi / double(main_paths) << std::endl;
	for (size_t i = 0; i < rats.size(); i++)
	{
		std::cout << mats[i].value << "   " <<  (1e8 / 10000.0) *  path_derivs[i] / double(main_paths);
		//CVA_delta_bumps_brute(i, 1e-8);
		CVA_delta_bumps(i, 1e-8);
	}

	delete mdl;
}

void pathwise_diff_in_main() {
	std::vector<Variable> mats = {
		//0.0,
		0.1, 3.0 / 12.0, 0.5,
		1,2,3,
		5,8,
		10,12,
		15,
		20,
		30 };
	std::vector<Variable> rats = {
		//0.035877,
		0.035877,  0.036576,	 0.037741,
		0.038560,	 0.038484,		0.038218,
		0.038186,	 0.038787,
		0.039264,	  0.039695,
		0.040222,
		0.040828,
		0.041475 };

	//Variable a = 0.164883; Variable sigma = 0.0280143;
	HullWhiteModel<Variable> *mdl = new HullWhiteModel<Variable>(mats, rats);


	Variable fixed_rate = 0.03;
	Variable maturity = 10;
	Variable fixed_tenor = 1;
	Variable float_tenor = 0.5;
	InterestRateSwap<Variable> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, fixed_rate);

	size_t pre_paths = 256;
	size_t main_paths = 1024;

	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);
	Timeline<Variable> timeline(0.0, maturity, 0.1);
	Matrix<Variable> r = mdl->simulate(timeline, pre_paths, mdl->r0, rng);
	Matrix<Variable> cf = IRS_cash_flows(swp, r, timeline);

	PowerShortRatePhi<Variable> basis(1);
	Matrix<Variable> beta_proxy = IRS_proxy_beta_calculator(swp, basis, r, cf, timeline);
	rng.setseeds(1, 2, 3, 4, 5, 6);

	std::vector<Node> copy = GlobalTape->nodes;

	double cva_poi = 0.0;
	std::vector<double> derivs;
	double of_interest = 0.0;
	rng.setseeds(1, 2, 3, 4, 5, 6);

	for (size_t i = 0; i < main_paths; i++)
	{
		Matrix<Variable> r2 = mdl->simulate(timeline, 1, mdl->r0, rng);
		Matrix<Variable> cf2 = IRS_cash_flows(swp, r2, timeline);
		Matrix<Variable> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);

		Variable cva = IRS_POI(swp, basis, beta_proxy, r2, cf2, v2, timeline);
		derivs = cva.calc_derivatives();
		of_interest += wrt(rats[4], derivs);

		cva_poi += cva.value;
		GlobalTape->nodes = copy;
	}
	std::cout << (100000000.0 / 10000.0) * of_interest / main_paths << "   " << cva_poi * 100000000.0 / main_paths; 


	delete mdl;
}



double final_example_full_sens_IRS10Y_d(size_t idx, bool cds_bump, bool fwd, double eps)
{
	std::vector<double> mats = { 0.0,0.1, 3.0 / 12.0, 0.5,1,2,3,	5,8,10,12,15,20, 30 };
	std::vector<double> rats = { 0.035877,0.035877,  0.036576, 0.037741,0.038560,0.038484,0.038218,0.038186, 0.038787,0.039264, 0.039695,0.040222,0.040828, 0.041475 };
	std::vector<double> tenor = { 0, 0.5, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30 };
	std::vector<double> spreads = { 0, 0.00154, 0.00186, 0.00277, 0.00477, 0.00488, 0.00543, 0.00751, 0.00831, 0.00955, 0.01058, 0.01188 };
	
	if (cds_bump) { 
		if (fwd) { spreads[idx] += eps; } else { spreads[idx] -= eps; } 
	}
	else {
		if (fwd) { rats[idx] += eps; }	else { rats[idx] -= eps; }
	}
	HullWhiteModel<double> *mdl = new HullWhiteModel<double>(mats, rats);
	creditObject<double> CreditModel(spreads, tenor, mdl);

	double fixed_rate = 0.0398985;
	double maturity = 10;
	double fixed_tenor = 1;
	double float_tenor = 0.5;
	InterestRateSwap<double> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, fixed_rate);
	Timeline<double> timeline(0.0, maturity, 0.1); // Time grid
	PowerShortRatePhi<double> basis(1); // Basis function

	size_t pre_paths = 256;
	size_t main_paths = 1024;

	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);

	// Pre simulation
	Matrix<double> r = mdl->simulate(timeline, pre_paths, mdl->r0, rng);
	Matrix<double> cf = IRS_cash_flows(swp, r, timeline);
	Matrix<double> beta_proxy = IRS_proxy_beta_calculator(swp, basis, r, cf, timeline);

	// Main simulation
	rng.setseeds(1, 2, 3, 4, 5, 6); // Resetting seeds
	Matrix<double> r2 = mdl->simulate(timeline, main_paths, mdl->r0, rng, true);
	Matrix<double> cf2 = IRS_cash_flows(swp, r2, timeline);
	Matrix<double> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);
	double cva_poi = 1e8 * (1 - 0.4) * POI_from_value_matrix(swp, basis, beta_proxy,
		r2, cf2, v2, timeline, CreditModel);
	delete mdl;
	return cva_poi;
}

long final_example_naive_full_sens_IRS10Y(bool print = true) 
{
	long before = GetTickCount64();
	std::vector<Variable> mats = {	0.0,0.1, 3.0 / 12.0, 0.5,1,2,3,	5,8,10,12,15,20, 30 };
	std::vector<Variable> rats = {	0.035877,0.035877,  0.036576, 0.037741,0.038560,0.038484,0.038218,0.038186, 0.038787,0.039264, 0.039695,0.040222,0.040828, 0.041475 };
	std::vector<Variable> tenor = { 0, 0.5, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30 }; 
	std::vector<Variable> spreads = { 0, 0.00154, 0.00186, 0.00277, 0.00477, 0.00488, 0.00543, 0.00751, 0.00831, 0.00955, 0.01058, 0.01188 };
	HullWhiteModel<Variable> *mdl = new HullWhiteModel<Variable>(mats, rats);
	creditObject<Variable> CreditModel(spreads, tenor, mdl);

	Variable fixed_rate = 0.0398985;
	Variable maturity = 10.0;
	Variable fixed_tenor = 1;
	Variable float_tenor = 0.5;
	InterestRateSwap<Variable> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, fixed_rate);

	Timeline<Variable> timeline(0.0, maturity, 0.1); // Time grid
	PowerShortRatePhi<Variable> basis(1); // Basis function

	size_t pre_paths = 256;
	size_t main_paths = 1024;

	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);

	// Pre simulation
	Matrix<Variable> r = mdl->simulate(timeline, pre_paths, mdl->r0, rng);
	Matrix<Variable> cf = IRS_cash_flows(swp, r, timeline);
	Matrix<Variable> beta_proxy = IRS_proxy_beta_calculator(swp, basis, r, cf, timeline);

	// Main simulation
	rng.setseeds(1, 2, 3, 4, 5, 6); // Resetting seeds
	Matrix<Variable> r2 = mdl->simulate(timeline, main_paths, mdl->r0, rng, true);
	Matrix<Variable> cf2 = IRS_cash_flows(swp, r2, timeline);
	Matrix<Variable> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);
	
	Variable cva_poi = 1e8 * (1-0.4) * POI_from_value_matrix(swp, basis, beta_proxy,
		r2, cf2, v2, timeline, CreditModel);
	
	std::vector<double> derivs = cva_poi.calc_derivatives();
	long after = GetTickCount64();
	
	long AAD_time_in_ms = after - before;
	
	// Cleanup
	delete mdl;
	GlobalTape->nodes.clear();
	
	if (print) {
		// Finite difference
		before = GetTickCount();
		std::vector<double> derivs_fd;
		double eps = 1e-10;
		for (size_t kk = 1; kk < mats.size(); kk++)
		{
			derivs_fd.push_back((final_example_full_sens_IRS10Y_d(kk, false, true, eps) - final_example_full_sens_IRS10Y_d(kk, false, false, eps)) / (2 * eps) / 10000);
		}

		for (size_t jj = 1; jj < tenor.size(); jj++)
		{
			derivs_fd.push_back((final_example_full_sens_IRS10Y_d(jj, true, true, eps) - final_example_full_sens_IRS10Y_d(jj, true, false, eps)) / (2 * eps) / 10000);
		}
		after = GetTickCount();

		std::cout << "CVA: " << cva_poi.value << std::endl;
		std::cout << "Mat   /  AAD Full    /  FD " << std::endl;
		for (size_t ii = 1; ii < mats.size(); ii++)
		{
			std::cout << mats[ii].value << "  /  " << wrt(rats[ii], derivs) / 10000 << "  / " << derivs_fd[ii - 1] << std::endl;
		}
		std::cout << "CDS sensitivities" << std::endl;
		for (size_t jj = 1; jj < tenor.size(); jj++)
		{
			std::cout << tenor[jj].value << " /  " << wrt(spreads[jj], derivs) / 10000 << " / " << derivs_fd[jj + mats.size() - 2] << std::endl;
		}
		std::cout << "AAD time in ms " << AAD_time_in_ms << " FD time in ms " << after - before << std::endl;
	}

	return AAD_time_in_ms;
}


long final_example_naive_cashflow_sens_IRS10Y(bool print = true)
{
	long before = GetTickCount64();
	std::vector<Variable> mats = { 0.0,0.1, 3.0 / 12.0, 0.5,1,2,3,	5,8,10,12,15,20, 30 };
	std::vector<Variable> rats = { 0.035877,0.035877,  0.036576, 0.037741,0.038560,0.038484,0.038218,0.038186, 0.038787,0.039264, 0.039695,0.040222,0.040828, 0.041475 };
	std::vector<Variable> tenor = { 0, 0.5, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30 };
	std::vector<Variable> spreads = { 0, 0.00154, 0.00186, 0.00277, 0.00477, 0.00488, 0.00543, 0.00751, 0.00831, 0.00955, 0.01058, 0.01188 };
	HullWhiteModel<Variable> *mdl = new HullWhiteModel<Variable>(mats, rats);
	creditObject<Variable> CreditModel(spreads, tenor, mdl);

	Variable fixed_rate = 0.0398985;
	Variable maturity = 10;
	Variable fixed_tenor = 1;
	Variable float_tenor = 0.5;
	InterestRateSwap<Variable> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, fixed_rate);
	Timeline<Variable> timeline(0.0, maturity, 0.1); // Time grid
	PowerShortRatePhi<Variable> basis(1); // Basis function

	size_t pre_paths = 256;
	size_t main_paths = 1024;

	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);

	// Pre simulation using doubles
	double maturity_double = as_double(maturity);
	std::vector<double> mats_double = as_double(mats);
	std::vector<double> rats_double = as_double(rats);
	HullWhiteModel<double> *mdl_double = new HullWhiteModel<double>(mats_double, rats_double);
	Timeline<double> timeline_double = as_double(timeline);
	InterestRateSwap<double> swp_double(fixed_tenor.value, maturity_double, float_tenor.value, maturity_double, mdl_double->r0, mdl_double, swp.payer, fixed_rate.value);

	Matrix<double> r = mdl_double->simulate(timeline_double, pre_paths, mdl_double->r0, rng);
	Matrix<double> cf = IRS_cash_flows(swp_double, r, timeline_double);
	PowerShortRatePhi<double> basis_double(basis.basis_funcs);
	Matrix<double> beta_proxy = IRS_proxy_beta_calculator(swp_double, basis_double, r, cf, timeline_double);

	// Main simulation
	rng.setseeds(1, 2, 3, 4, 5, 6); // Resetting seeds
	Matrix<Variable> r2 = mdl->simulate(timeline, main_paths, mdl->r0, rng, true);
	Matrix<Variable> cf2 = IRS_cash_flows(swp, r2, timeline);
	Matrix<Variable> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);

	Variable cva_poi = 1e8 * (1 - 0.4) * POI_from_value_matrix(swp, basis, beta_proxy,
		r2, cf2, v2, timeline, CreditModel);

	std::vector<double> derivs = cva_poi.calc_derivatives();
	long after = GetTickCount64();

	if (print) {
		std::cout << after - before << "ms" << std::endl;

		std::cout << "CVA: " << cva_poi.value << std::endl;
		for (size_t ii = 1; ii < mats.size(); ii++)
		{
			std::cout << mats[ii].value << "   " << wrt(rats[ii], derivs) / 10000 << std::endl;
		}
		std::cout << "Credit sensitivities " << std::endl;
		for (size_t ii = 1; ii < tenor.size(); ii++)
		{
			std::cout << tenor[ii].value << "   " << wrt(spreads[ii], derivs) / 10000 << std::endl;
		}
		std::cout << "Time in ms " << after - before << std::endl;
	}
	delete mdl;

	// Clean up
	GlobalTape->nodes.clear();

	return after - before;
}

long final_example_pathwise_cashflow_sens_IRS10Y(size_t paths_per_sim, bool print = true)
{
	long before = GetTickCount();
	std::vector<Variable> mats = { 0.0,0.1, 3.0 / 12.0, 0.5,1,2,3,5,8,10,12 ,15,20, 30 };
	std::vector<Variable> rats = { 0.035877,0.035877,  0.036576, 0.037741,0.038560,0.038484,0.038218,0.038186, 0.038787,0.039264, 0.039695,0.040222,0.040828, 0.041475 };
	std::vector<Variable> tenor = { 0, 0.5, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30 };
	std::vector<Variable> spreads = { 0, 0.00154, 0.00186, 0.00277, 0.00477, 0.00488, 0.00543, 0.00751, 0.00831, 0.00955, 0.01058, 0.01188 };
	HullWhiteModel<Variable> *mdl = new HullWhiteModel<Variable>(mats, rats);
	creditObject<Variable> CreditModel(spreads, tenor, mdl);

	Variable fixed_rate = 0.0398985;
	Variable maturity = 10;
	Variable fixed_tenor = 1;
	Variable float_tenor = 0.5;
	InterestRateSwap<Variable> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, fixed_rate);
	Timeline<Variable> timeline(0.0, maturity, 0.1); // Time grid
	PowerShortRatePhi<Variable> basis(1); // Basis function

	size_t pre_paths = 256;
	size_t main_paths = 1024 * 8;

	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);

	// Pre simulation using doubles
	double maturity_double = as_double(maturity);
	std::vector<double> mats_double = as_double(mats);
	std::vector<double> rats_double = as_double(rats);
	HullWhiteModel<double> *mdl_double = new HullWhiteModel<double>(mats_double, rats_double);
	Timeline<double> timeline_double = as_double(timeline);
	InterestRateSwap<double> swp_double(fixed_tenor.value, maturity_double, float_tenor.value, maturity_double, mdl_double->r0, mdl_double, swp.payer, fixed_rate.value);

	Matrix<double> r = mdl_double->simulate(timeline_double, pre_paths, mdl_double->r0, rng);
	Matrix<double> cf = IRS_cash_flows(swp_double, r, timeline_double);
	PowerShortRatePhi<double> basis_double(basis.basis_funcs);
	Matrix<double> beta_proxy = IRS_proxy_beta_calculator(swp_double, basis_double, r, cf, timeline_double);

	// Main simulation
	rng.setseeds(1, 2, 3, 4, 5, 6); // Resetting seeds

	std::vector<Node> tape_copy = GlobalTape->nodes;
	size_t simuls_needed = size_t(double(main_paths) / double(paths_per_sim));

	std::vector<double> zero_rate_derivatives(mats.size());
	std::vector<double> credit_derivatives(tenor.size());
	double cva = 0.0;
	for (size_t i = 0; i < simuls_needed; i++)
	{
		Matrix<Variable> r2 = mdl->simulate(timeline, paths_per_sim, mdl->r0, rng, true);
		Matrix<Variable> cf2 = IRS_cash_flows(swp, r2, timeline);
		Matrix<Variable> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);
		Variable cva_poi = (1 - 0.4) * POI_from_value_matrix(swp, basis, beta_proxy,
			r2, cf2, v2, timeline, CreditModel);
		cva += cva_poi.value;
		std::vector<double> derivs = cva_poi.calc_derivatives();
		
		for (size_t ii = 0; ii < mats.size(); ii++)	zero_rate_derivatives[ii] += wrt(rats[ii], derivs);
		for (size_t ii = 0; ii < tenor.size(); ii++) credit_derivatives[ii] += wrt(spreads[ii], derivs);
		
		GlobalTape->nodes = tape_copy;
	}


	for (size_t ii = 0; ii < mats.size(); ii++)	zero_rate_derivatives[ii] *= 1e4 / double(simuls_needed);
	for (size_t ii = 0; ii < tenor.size(); ii++) credit_derivatives[ii] *=  1e4 / double(simuls_needed);

	long after = GetTickCount();
	long time_in_ms = after - before;
	if (print) {
		std::cout << "CVA " << 1e8 / double(simuls_needed) * cva << "  " << paths_per_sim << "  Time in ms " << time_in_ms << std::endl;

		// Print
		for (size_t i = 0; i < rats.size(); i++)std::cout << mats[i].value << "  " << zero_rate_derivatives[i] << std::endl;
		for (size_t i = 0; i < tenor.size(); i++) std::cout << tenor[i].value << "  " << credit_derivatives[i] << std::endl;
	}

	GlobalTape->nodes.clear();
	delete mdl;
	return time_in_ms;
}

long delta_vector_only_finite_difference() {
	double eps = 1e-10;
	std::vector<double> derivs_fd;
	long before = GetTickCount();
	for (size_t kk = 1; kk < 14; kk++)
		derivs_fd.push_back((final_example_full_sens_IRS10Y_d(kk, false, true, eps) - final_example_full_sens_IRS10Y_d(kk, false, false, eps)) / (2 * eps) / 10000);

	for (size_t jj = 1; jj < 12; jj++)
		derivs_fd.push_back((final_example_full_sens_IRS10Y_d(jj, true, true, eps) - final_example_full_sens_IRS10Y_d(jj, true, false, eps)) / (2 * eps) / 10000);
	long after = GetTickCount();

	long time_in_ms = after - before;
	//std::cout << "Time in ms " << time_in_ms << std::endl;
	return time_in_ms;
}


void portfolio_two_npv_cva_with_doubles() {
	std::vector<double> mats = { 0.0,0.1, 3.0 / 12.0, 0.5,1,2,3,	5,8,10,12,15,20, 30 };
	std::vector<double> rats = { 0.035877,0.035877,  0.036576, 0.037741,0.038560,0.038484,0.038218,0.038186, 0.038787,0.039264, 0.039695,0.040222,0.040828, 0.041475 };
	std::vector<double> tenor = { 0, 0.5, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30 };
	std::vector<double> spreads = { 0, 0.00154, 0.00186, 0.00277, 0.00477, 0.00488, 0.00543, 0.00751, 0.00831, 0.00955, 0.01058, 0.01188 };
	HullWhiteModel<double> *mdl = new HullWhiteModel<double>(mats, rats);
	creditObject<double> CreditModel(spreads, tenor, mdl);


	std::vector<double> fixed_rates = { 0.03, 0.04, 0.0398985, 0.02, 0.05 };
	std::vector<double> maturity = { 30, 21, 10, 5, 2 }; 
	std::vector<bool>   payer = { true, false, true, false, false };
	std::vector<double> notionals = { 0.10, 0.20, 0.20, 0.25, 0.25 }; // This portfolio has close to NPV its actually around 1000
	for (size_t i = 0; i < notionals.size(); i++) notionals[i] *= 5e8; // Now the 10Y par swap has 100M notional
	double fixed_tenor = 1;
	double float_tenor = 0.5;
	std::vector<Timeline<double>> timelines;
	for (size_t i = 0; i < maturity.size(); i++) timelines.push_back(Timeline<double>(0.0, maturity[i], 0.1));

	std::vector<InterestRateSwap<double>> swaps; // Should build Portfolio class to hold many different derivatives
	double npv = 0.0;
	for (size_t jj = 0; jj < payer.size(); jj++)
	{
		swaps.push_back(InterestRateSwap<double>(
			fixed_tenor, maturity[jj], float_tenor,
			maturity[jj], mdl->r0, mdl, payer[jj], fixed_rates[jj]));
		npv += notionals[jj]*swaps[jj].npv(mats[0], mdl->r0);
	}
	PowerShortRatePhi<double> basis(2);
	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);

	// Pre simulation
	size_t pre_paths = 256;
	Matrix<double> r = mdl->simulate(timelines[0], pre_paths, mdl->r0, rng, true);
	Matrix<double> cf(pre_paths, timelines[0].size());
	for (size_t kk = 0; kk < cf.internal_vector.size(); kk++) cf.internal_vector[kk] = 0.0;

	Matrix<double> tmp_cf;
	for (size_t i = 0; i < payer.size(); i++)
	{
		IRS_cash_flows(swaps[i], r, timelines[i], tmp_cf); // Fills tmp_cf with swaps[i] cashflows
		for (size_t path = 0; path < tmp_cf.nrows; path++)
		{
			for (size_t date = 0; date < tmp_cf.ncols; date++)
			{
				cf[path][date] += notionals[i] * tmp_cf[path][date];
			}
		}
	}
	Matrix<double> beta = proxy_beta_calculator(swaps[0], basis, r, cf, timelines[0]); // Currently this needs a swap for the phi
																					   // This should be altered.

	// Main simulation
	size_t main_paths = 8196;
	rng.setseeds(1, 2, 3, 4, 5, 6);
	r = mdl->simulate(timelines[0], main_paths, mdl->r0, rng, true);

	Matrix<double> tmp_cf_main;
	Matrix<double> cf_main(main_paths, timelines[0].size()); 
	for (size_t kk = 0; kk < cf_main.internal_vector.size(); kk++) { cf_main.internal_vector[kk] = 0.0; }

	for (size_t i = 0; i < payer.size(); i++)
	{
		IRS_cash_flows(swaps[i], r, timelines[i], tmp_cf_main); // Fills tmp_cf with swaps[i] cashflows
		for (size_t path = 0; path < tmp_cf_main.nrows; path++)
		{
			for (size_t date = 0; date < tmp_cf_main.ncols; date++)
			{
				cf_main[path][date] += notionals[i] * tmp_cf_main[path][date];
			}
		}
	}

	Matrix<double> value = IRS_proxy_value(swaps[0], basis, beta, r, cf_main, timelines[0]);
	
	double cva_poi = (1.0 - 0.4) * POI_from_value_matrix(swaps[0], basis, beta, r, cf_main, value, timelines[0], CreditModel);
	std::cout << "NPV " << npv << std::endl;
	std::cout << "CVA " << cva_poi << std::endl;
	delete mdl;
}


void portfolio_two_cashflow(size_t paths_per_tape, bool print = true) {
	std::vector<Variable> mats = { 0.0,0.1, 3.0 / 12.0, 0.5,1,2,3,	5,8,10,12,15,20, 30 };
	std::vector<Variable> rats = { 0.035877,0.035877,  0.036576, 0.037741,0.038560,0.038484,0.038218,0.038186, 0.038787,0.039264, 0.039695,0.040222,0.040828, 0.041475 };
	std::vector<Variable> tenor = { 0, 0.5, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30 };
	std::vector<Variable> spreads = { 0, 0.00154, 0.00186, 0.00277, 0.00477, 0.00488, 0.00543, 0.00751, 0.00831, 0.00955, 0.01058, 0.01188 };
	HullWhiteModel<Variable> *mdl = new HullWhiteModel<Variable>(mats, rats);
	creditObject<Variable> CreditModel(spreads, tenor, mdl);
 
	std::vector<Variable> fixed_rates = { 0.03, 0.04, 0.0398985, 0.02, 0.05 };
	std::vector<Variable> maturity = { 30, 21, 10, 5, 2 };
	std::vector<bool>   payer = { true, false, true, false, false };
	std::vector<Variable> notionals = { 0.10, 0.20, 0.20, 0.25, 0.25 }; // This portfolio has close to zero NPV
	//for (size_t i = 0; i < notionals.size(); i++) notionals[i] *= 5e8; // Now the 10Y par swap has 100M notional
	Variable fixed_tenor = 1;
	Variable float_tenor = 0.5;
	std::vector<Timeline<Variable>> timelines;
	for (size_t i = 0; i < maturity.size(); i++) timelines.push_back(Timeline<Variable>(0.0, maturity[i], 0.1));

	std::vector<InterestRateSwap<Variable>> swaps; // Should build Portfolio class to hold many different derivatives
	for (size_t jj = 0; jj < payer.size(); jj++)
	{
		swaps.push_back(InterestRateSwap<Variable>(
			fixed_tenor, maturity[jj], float_tenor,
			maturity[jj], mdl->r0, mdl, payer[jj], fixed_rates[jj]));
	}
	PowerShortRatePhi<Variable> basis(1);
	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);

	///////////////// Make everything into doubles... /////////////////
	///////////////// Poor practice but does the trick/////////////////
	std::vector<double> mats_double = as_double(mats);
	std::vector<double> rats_double = as_double(rats);
	HullWhiteModel<double> *mdl_doubles = new HullWhiteModel<double>(mats_double, rats_double);
	std::vector<double> fixed_rates_doubles = as_double(fixed_rates);
	std::vector<double> maturity_doubles = as_double(maturity);
	std::vector<double> notionals_doubles = as_double(notionals);
	double fixed_tenor_doubles = as_double(fixed_tenor); double float_tenor_doubles = as_double(float_tenor);
	std::vector<Timeline<double>> timelines_doubles;
	for (size_t i = 0; i < maturity.size(); i++) timelines_doubles.push_back(Timeline<double>(0.0, maturity_doubles[i], 0.1));
	std::vector<InterestRateSwap<double>> swaps_doubles;
	double npvs = 0.0;
	for (size_t jj = 0; jj < payer.size(); jj++)
	{
		swaps_doubles.push_back(InterestRateSwap<double>(
			fixed_tenor_doubles, maturity_doubles[jj], float_tenor_doubles,
			maturity_doubles[jj], mdl_doubles->r0, mdl_doubles, payer[jj], fixed_rates_doubles[jj]));
		npvs += notionals_doubles[jj] * swaps_doubles[jj].npv(mats_double[0], mdl_doubles->r0);
	}
	PowerShortRatePhi<double> basis_doubles(basis.basis_funcs);

	// Begin pre-simulation
	size_t pre_paths = 256;

	Matrix<double> r_pre = mdl_doubles->simulate(timelines_doubles[0], pre_paths, mdl_doubles->r0, rng, true);
	Matrix<double> cf(pre_paths, timelines_doubles[0].size());
	for (size_t kk = 0; kk < cf.internal_vector.size(); kk++) cf.internal_vector[kk] = 0.0;

	Matrix<double> tmp_cf;
	for (size_t i = 0; i < payer.size(); i++)
	{
		IRS_cash_flows(swaps_doubles[i], r_pre, timelines_doubles[i], tmp_cf); // Fills tmp_cf with swaps[i] cashflows
		for (size_t path = 0; path < tmp_cf.nrows; path++)
		{
			for (size_t date = 0; date < tmp_cf.ncols; date++)
			{
				cf[path][date] += notionals_doubles[i] * tmp_cf[path][date];
			}
		}
	}
	Matrix<double> beta = proxy_beta_calculator(swaps_doubles[0], basis_doubles, r_pre, cf, timelines_doubles[0]);
	////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////

	// Main simulation
	rng.setseeds(1, 2, 3, 4, 5, 6); // Resetting seeds

	size_t main_paths = 8196;
	std::vector<Node> tape_copy = GlobalTape->nodes;
	size_t simuls_needed = size_t(double(main_paths) / double(paths_per_tape));

	std::vector<double> zero_rate_derivatives(mats.size());
	std::vector<double> credit_derivatives(tenor.size());
	double cva = 0.0;
	for (size_t i = 0; i < simuls_needed; i++)
	{
		Matrix<Variable> r2 = mdl->simulate(timelines[0], paths_per_tape, mdl->r0, rng, true);

		Matrix<Variable> cf(paths_per_tape, timelines_doubles[0].size());
		for (size_t kk = 0; kk < cf.internal_vector.size(); kk++) cf.internal_vector[kk] = 0.0;
		Matrix<Variable> tmp_cf;
		for (size_t i = 0; i < payer.size(); i++)
		{
			IRS_cash_flows(swaps[i], r2, timelines[i], tmp_cf); // Fills tmp_cf with swaps[i] cashflows
			for (size_t path = 0; path < tmp_cf.nrows; path++)
			{
				for (size_t date = 0; date < tmp_cf.ncols; date++)
				{
					cf[path][date] += notionals[i] * tmp_cf[path][date];
				}
			}
		}
		Matrix<Variable> v2 = IRS_proxy_value(swaps[0], basis, beta, r2, cf, timelines[0]);
		Variable cva_poi = (1 - 0.4) * POI_from_value_matrix(swaps[0], basis, beta,
			r2, cf, v2, timelines[0], CreditModel);

		cva += cva_poi.value;
		std::vector<double> derivs = cva_poi.calc_derivatives();

		for (size_t ii = 0; ii < mats.size(); ii++)	zero_rate_derivatives[ii] += wrt(rats[ii], derivs);
		for (size_t ii = 0; ii < tenor.size(); ii++) credit_derivatives[ii] += wrt(spreads[ii], derivs);

		GlobalTape->nodes = tape_copy;
	}

	for (size_t ii = 0; ii < mats.size(); ii++)	zero_rate_derivatives[ii]	/= double(simuls_needed);
	for (size_t ii = 0; ii < tenor.size(); ii++) credit_derivatives[ii]		/= double(simuls_needed);
	
	if (print) {
		std::cout << "NPV " << 5e8 * npvs << std::endl;
		std::cout << "CVA " << 5e8 / double(simuls_needed) * cva << std::endl;

		// Print
		for (size_t i = 1; i < rats.size(); i++)std::cout << mats[i].value << "  " << (5e8 * 1e-4)  * zero_rate_derivatives[i] << std::endl;
		for (size_t i = 1; i < tenor.size(); i++) std::cout << tenor[i].value << "  " << (5e8 * 1e-4) * credit_derivatives[i] << std::endl;
	}
	delete mdl_doubles;
	delete mdl;
}
