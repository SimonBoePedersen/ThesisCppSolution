#pragma once


void CVA(size_t inp, size_t paths) {
	std::vector<Variable> mats = { 0.0, 0.1, 3.0 / 12.0, 0.5, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 };
	std::vector<Variable> rats = { 0.035877, 0.035877,  0.036576,	 0.037741,	 0.038560,	 0.038484,		0.038218,	 0.038126,	 0.038186,	 0.038342,	 0.038552,		0.038787,	 0.039028,	 0.039264,	 0.039487,	 0.039695,		0.039887,	 0.040062,	 0.040222,	 0.040368,	 0.040500,		0.040619,	 0.040729,	 0.040828,	 0.040919,	 0.041002,		0.041079,	 0.041149,	 0.041214,	 0.041274,	 0.041330,		0.041382,	 0.041430,	 0.041475 };

	HullWhiteModel<Variable> *mdl = new HullWhiteModel<Variable>(mats, rats);

	size_t pre_paths = 2 * 256;
	size_t main_paths = paths;
	Variable maturity = 3.0;

	Variable fixed_tenor = 1;
	Variable float_tenor = 0.5;

	InterestRateSwap<Variable> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, 1.0);

	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);
	Timeline<Variable> timeline(0.0, maturity, 0.1);
	Matrix<Variable> r = mdl->simulate(timeline, pre_paths, mdl->r0, rng);
	Matrix<Variable> cf = IRS_cash_flows(swp, r, timeline);

	PowerShortRatePhi<Variable> basis(2);
	Matrix<Variable> beta_proxy = IRS_proxy_beta_calculator(swp, basis, r, cf, timeline);
	rng.setseeds(1, 2, 3, 4, 5, 6);
	Matrix<Variable> r2 = mdl->simulate(timeline, main_paths, mdl->r0, rng);
	Matrix<Variable> cf2 = IRS_cash_flows(swp, r2, timeline);
	Matrix<Variable> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);

	Variable cva_poi = IRS_POI(swp, basis, beta_proxy, r2, cf2, v2, timeline);
	std::cout << "CVA POI: " << cva_poi.value << std::endl;
	std::vector<double> derivs = cva_poi.calc_derivatives();
	double wrttwo = wrt(rats[inp], derivs);
	std::cout << "dr2 " << wrttwo << std::endl;
	delete mdl;
}


double CVA_POI_all_doubles(size_t inp, size_t paths, double eps) {
	//std::vector<double> mats = { 0.0,0.1, 3.0 / 12.0, 0.5,1,2,3,5,8,10,12,15,20, 30 };
	//std::vector<double> rats = { 0.035877,0.035877,  0.036576, 0.037741,0.038560,0.038484,0.038218,0.038186, 0.038787,0.039264, 0.039695,0.040222,0.040828, 0.041475 };
	//rats[inp] += eps;
	//HullWhiteModel<double> *mdl = new HullWhiteModel<double>(mats, rats);
	VasicekModel<double> *mdl = new VasicekModel<double>(0.04, 0.1, 0.015, 0.04);

	size_t pre_paths = 256*2;
	size_t main_paths = paths;
	double maturity = 10.0;
	double fixed_tenor = 1;
	double float_tenor = 0.5;

	InterestRateSwap<double> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, 1.0);

	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);
	Timeline<double> timeline(0.0, maturity, 0.1);
	Matrix<double> r = mdl->simulate(timeline, pre_paths, mdl->r0, rng);
	Matrix<double> cf = IRS_cash_flows(swp, r, timeline);

	PowerShortRatePhi<double> basis(1);
	Matrix<double> beta_proxy = IRS_proxy_beta_calculator(swp, basis, r, cf, timeline);
	rng.setseeds(1, 2, 3, 4, 5, 6);
	Matrix<double> r2 = mdl->simulate(timeline, main_paths, mdl->r0, rng);
	Matrix<double> cf2 = IRS_cash_flows(swp, r2, timeline);
	Matrix<double> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);

	std::vector<double> epe = EPE_from_value_matrix(v2,timeline);
	double pcva = 0.0;
	for (size_t i = 1; i < epe.size(); i++)
	{
		pcva += epe[i] * 0.02*0.1;
	}
	std::cout << "pur prox " << pcva * 1e8 * 0.6 << std::endl;

	double cva_poi = IRS_POI(swp, basis, beta_proxy, r2, cf2, v2, timeline);
	delete mdl;
	return cva_poi;
}


template<typename T>
T IRS_CVA_BRUTE(std::vector<T>& mats,
	std::vector<T>& rats, size_t paths) {
	//HullWhiteModel<T> *mdl = new HullWhiteModel<T>(mats, rats);
	VasicekModel<T> *mdl = new VasicekModel<T>(0.04, 0.1, 0.015, 0.04);

	size_t main_paths = paths;
	T maturity = 10.0;
	T fixed_tenor = 1.0;
	T float_tenor = 0.5;
	InterestRateSwap<T> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, 1.0);
	mrg32 rng(1, 2, 3, 4, 5, 6);
	Timeline<T> timeline(0.0, maturity, 0.1);
	Matrix<T> r = mdl->simulate(timeline, main_paths, mdl->r0, rng);
	Matrix<T> v = IRS_bruteforce(swp, r, timeline);
	std::vector<T> EE_brute = EPE_from_value_matrix(v, timeline);
	T cva_brute = 0.0;
	for (size_t i = 1; i < EE_brute.size(); i++)
	{
		cva_brute = cva_brute + EE_brute[i] * 0.02 * 0.1;
	}
	delete mdl;
	return cva_brute;
}


void CVA_doubles_in_proxy_generation(size_t inp, size_t paths) {
	std::vector<Variable> mats = { 0.0, 0.1, 3.0 / 12.0, 0.5, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 };
	std::vector<Variable> rats = { 0.035877, 0.035877,  0.036576,	 0.037741,	 0.038560,	 0.038484,		0.038218,	 0.038126,	 0.038186,	 0.038342,	 0.038552,		0.038787,	 0.039028,	 0.039264,	 0.039487,	 0.039695,		0.039887,	 0.040062,	 0.040222,	 0.040368,	 0.040500,		0.040619,	 0.040729,	 0.040828,	 0.040919,	 0.041002,		0.041079,	 0.041149,	 0.041214,	 0.041274,	 0.041330,		0.041382,	 0.041430,	 0.041475 };

	HullWhiteModel<Variable> *mdl = new HullWhiteModel<Variable>(mats, rats);
	
	size_t pre_paths = 2 * 256;
	size_t main_paths = paths;
	Variable maturity = 3.0;

	Variable fixed_tenor = 1;
	Variable float_tenor = 0.5;

	InterestRateSwap<Variable> swp(fixed_tenor, maturity, float_tenor, maturity, mdl->r0, mdl, true, 1.0);
	Timeline<Variable> timeline(0.0, maturity, 0.1);
	PowerShortRatePhi<Variable> basis(2);
	mrg32 rng(1234, 1234, 1234, 12345, 12346, 12346);

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
	delete mdl_double;

	// Main simulation
	rng.setseeds(1, 2, 3, 4, 5, 6);
	Matrix<Variable> r2 = mdl->simulate(timeline, main_paths, mdl->r0, rng);
	Matrix<Variable> cf2 = IRS_cash_flows(swp, r2, timeline);
	Matrix<Variable> v2 = IRS_proxy_value(swp, basis, beta_proxy, r2, cf2, timeline);

	Variable cva_poi = IRS_POI(swp, basis, beta_proxy, r2, cf2, v2, timeline);
	std::cout << "CVA POI: " << cva_poi.value << std::endl;
	std::vector<double> derivs = cva_poi.calc_derivatives();
	double wrttwo = wrt(rats[inp], derivs);
	std::cout << "dr2 " << wrttwo << std::endl;
	
}