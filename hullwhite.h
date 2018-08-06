#pragma once

template <typename T>
class HullWhiteModel : public BaseRateModel<T>
{
	// Hull White short rate model object.
	// Fits a Hull-White model to the input curve with the nelder-mead method
	// Most imporant are the methods which returns forward rates, bond prices and simulate short rates
public:
	HullWhiteModel() {};
	HullWhiteModel(std::vector<T>& maturities, std::vector<T>& zero_rates);
	HullWhiteModel(std::vector<T>& maturities, std::vector<T>& zero_rates, T& _a_, T& _sigma_);
	~HullWhiteModel() {};

	T forward_rate(T start, T maturity, T& short_rate) {
		T df_start;
		double TOL = 0.000001;
		if (abs(start - 0.0) < TOL) {
			df_start = 1.0;
		}
		else {
			df_start = disc_factor(0, start, short_rate);
		}
		T df_mat = disc_factor(0, maturity, short_rate);
		T out = (df_start / df_mat - 1.0) / (maturity - start);
		return out;
	}

	T theta(T t) {
		T out = fM_spline_obj.interpolate_derivative(t) +
			a * fM_spline_obj.interpolate(t) + sigma*sigma / (2.0 * a) * (1.0 - exp(-2.0 * a*t));
		return out;
	}

	T disc_factor(T start, T end, T& short_rate) {
		return bond_price(start, end, short_rate);
	}

	T bond_price(T start, T maturity, T& short_rate) {
		T BBB = B(start, maturity);
		T AAA = A(start, maturity, BBB);
		return exp(AAA + BBB * short_rate);
	}

	T spot_libor(T current_time, T maturity, T& short_rate) {
		T ptT = bond_price(current_time, maturity, short_rate);
		return (1 - ptT) / ((maturity - current_time)*ptT);
	}

	T B(T start, T mat) {
		T out = (exp(-a*(mat - start)) - 1.0) / a;
		return out;
	}

	T A(T start, T mat, T b) {
		T pMstart = exp(-start * zeroCurve_obj.interpolate(start));
		T pMmat = exp(-mat * zeroCurve_obj.interpolate(mat));

		T f = fM_spline_obj.interpolate(mat) * b;
		T lg = log(pMmat / pMstart);
		T last = sigma*sigma / (4.0*a) * b * b * (exp(-2.0*a*mat) - 1.0);
		T test = -f + lg + last;
		T out = -fM_spline_obj.interpolate(mat) * b + log(pMmat / pMstart) +
			sigma*sigma / (4.0*a) * b * b * (exp(-2.0*a*mat) - 1.0);
		return out;
	}

	T sumsq_for_calibration(T _a_, T _sigma_) {
		T out = 0.0;
		T holder_a = a;
		T holder_sig = sigma;
		a = _a_; sigma = _sigma_;
		for (size_t i = 0; i < mats.size(); i++)
		{
			T tmp_p = bond_price(0, mats[i], r0);
			out = out + (tmp_p - pM[i]) * (tmp_p - pM[i]);
		}
		a = holder_a; sigma = holder_sig;
		return out;
	}

	Matrix<T> simulate(Timeline<T>& timeline,
		const size_t paths, const T& r0, mrg32& rng, bool byrow = false);

	void twodimsimplex(); // This method is the Nelder-Mead optimization on the variables

	T r0;
	T a;
	T sigma;
	std::vector<T> fM; //Foward curve
	CubicSpline<T> fM_spline_obj; //Splined forward curve obj
	CubicSpline<T> zeroCurve_obj; //Zero curve, for discounting from t_0
	std::vector<T> pM; //Observed prices in market.
	std::vector<T> mats;
};

template<typename T>
inline HullWhiteModel<T>::HullWhiteModel(std::vector<T>& maturities, std::vector<T>& zero_rates)
{
	//Just a guess on parameter. Arbitrary. Has to be sensible however
	a = 0.15;
	sigma = 0.03;

	mats = maturities;
	r0 = zero_rates[0];
	zeroCurve_obj = CubicSpline<T>(mats, zero_rates);
	pM.resize(mats.size());
	fM.resize(mats.size());
	for (size_t i = 0; i < mats.size(); i++)
	{
		pM[i] = exp(-mats[i] * zero_rates[i]);
		T bumped = mats[i] + 0.001;
		T df = exp(-mats[i] * zeroCurve_obj.interpolate(mats[i]));
		T df2 = exp(-(bumped)* zeroCurve_obj.interpolate(bumped));
		fM[i] = (df / df2 - 1) / 0.001;
	}
	fM_spline_obj = CubicSpline<T>(mats, fM);
	twodimsimplex(); // Fit the model
}

template<typename T>
inline HullWhiteModel<T>::HullWhiteModel(std::vector<T>& maturities, std::vector<T>& zero_rates, T& _a_, T& _sigma_)
{
	//Just a guess on parameter
	a = _a_;
	sigma = _sigma_;

	mats = maturities;
	r0 = zero_rates[0];
	zeroCurve_obj = CubicSpline<T>(mats, zero_rates);
	pM.resize(mats.size());
	fM.resize(mats.size());
	for (size_t i = 0; i < mats.size(); i++)
	{
		pM[i] = exp(-mats[i] * zero_rates[i]);
		T bumped = mats[i] + 0.001;
		T df = exp(-mats[i] * zeroCurve_obj.interpolate(mats[i]));
		T df2 = exp(-(bumped)* zeroCurve_obj.interpolate(bumped));
		fM[i] = (df / df2 - 1) / 0.001;
	}
	fM_spline_obj = CubicSpline<T>(mats, fM);

}

template<typename T>
class Vertex
{
public:

	Vertex() {};
	Vertex(T _a_, T _sigma_) : a(_a_), sigma(_sigma_) {};

	void update_value(HullWhiteModel<T>& mdl) { value = mdl.sumsq_for_calibration(a, sigma); }
	~Vertex() {};
	T value;
	T a;
	T sigma;

	bool operator < (const Vertex<T>& str) const
	{
		return (value < str.value);
	}
};

template<typename T>
bool sortVertex(Vertex<T>& i, Vertex<T>& j) {
	return (i.value < j.value);
}

template<typename T>
void HullWhiteModel<T>::twodimsimplex() 
{
	// Create simplex
	std::vector<Vertex<T>> pts = { Vertex<T>(a, sigma),
		Vertex<T>(a + 0.01, sigma),
		Vertex<T>(a, sigma + 0.001) };
	for (size_t i = 0; i < 3; i++)
	{
		pts[i].update_value(*this);
	}
	Vertex<T> Reflection;
	Vertex<T> Extension;
	Vertex<T> Centroid;
	Vertex<T> Contraction;

	bool mybool = false;
	int n = 1;
	while (!mybool && n < 100) {
		//Start by sorting
		std::sort(pts.begin(), pts.end());

		if (pts[2].value - pts[0].value < 0.000000001) { mybool = true; break; }

		Centroid.a = (pts[0].a + pts[1].a)*0.5; Centroid.sigma = (pts[0].sigma + pts[1].sigma)*0.5;

		Reflection.a = 2.0 * Centroid.a - pts[2].a;
		Reflection.sigma = 2.0 * Centroid.sigma - pts[2].sigma;
		Reflection.update_value(*this);

		if (Reflection < pts[1] && pts[0] < Reflection) {
			pts[2] = Reflection;
		}
		else if (Reflection < pts[0]) {
			//Extending
			Extension.a = 2 * Reflection.a - Centroid.a;
			Extension.sigma = 2 * Reflection.sigma - Centroid.sigma;
			Extension.update_value(*this);
			if (Extension < Reflection) {
				pts[2] = Extension;
			}
			else {
				pts[2] = Reflection;
			}
		}
		else if (pts[1] < Reflection) {
			if (pts[2] < Reflection) {
				Contraction.a = -0.5 * pts[2].a + 0.5 * Centroid.a;
				Contraction.sigma = -0.5 * pts[2].sigma + 0.5 * Centroid.sigma;
				Contraction.update_value(*this);
				if (Contraction < Reflection) {
					pts[3] = Contraction;
				}
				else {
					//Shrinking
					pts[0].a = (pts[0].a + pts[2].a) * 0.5;
					pts[1].a = (pts[1].a + pts[2].a) * 0.5;
					pts[0].sigma = (pts[0].sigma + pts[2].sigma) * 0.5;
					pts[1].sigma = (pts[1].sigma + pts[2].sigma) * 0.5;
					pts[0].update_value(*this);
					pts[1].update_value(*this);
				}
			}
			else {
				Contraction.a = 0.5 * Centroid.a - 0.5 * Reflection.a;
				Contraction.sigma = 0.5 * Centroid.sigma - 0.5 * Reflection.sigma;
				Contraction.update_value(*this);
				if (Contraction < Reflection) {
					pts[2] = Contraction;
				}
				else {
					//Shrinking
					pts[0].a = (pts[0].a + pts[2].a) * 0.5;
					pts[1].a = (pts[1].a + pts[2].a) * 0.5;
					pts[0].sigma = (pts[0].sigma + pts[2].sigma) * 0.5;
					pts[1].sigma = (pts[1].sigma + pts[2].sigma) * 0.5;
					pts[0].update_value(*this);
					pts[1].update_value(*this);
				}
			}
		}
		n++;
	}
	if (mybool) {
		a = (pts[0].a + pts[1].a + pts[2].a) / 3.0;
		sigma = (pts[0].sigma + pts[1].sigma + pts[2].sigma) / 3.0;
		//std::cout << "Calibration complete after " << n << " iterations." << std::endl;
	}
	else {
		std::cout << "Calibration failed." << std::endl;
	}
}

template<typename T>
inline Matrix<T> HullWhiteModel<T>::simulate(Timeline<T>& timeline, const size_t paths, const T& r0, mrg32& rng, bool byrow = false)
{
	// Returns a matrix of simulated short rate paths from the Hull-White model.
	// In the returned matrix each row represent a path
	// Can be filled out either by row or by column. This has impact on the order which 
	// the random noise is used (necessary to control the random noise for parallel computations)
	size_t steps = timeline.line.size();
	Matrix<T> out(paths, steps);

	for (size_t i = 0; i < paths; i++) { out[i][0] = r0; }

	if (byrow) {
		for (size_t path = 0; path < paths; path++)
		{
			for (size_t date = 1; date < steps; date++)
			{
				double U = rng.next();
				T dt = timeline[date] - timeline[date - 1];
				out[path][date] = out[path][date - 1] +
					(theta(timeline[date]) - a * out[path][date - 1]) * dt +
					sigma * sqrt(dt) * invNormalCdf(U);
			}
		}
	}
	else 
	{
		for (size_t date = 1; date < steps; date++)
		{
			T dt = timeline[date] - timeline[date - 1];
			for (size_t path = 0; path < paths; path++)
			{
				double U = rng.next();
				out[path][date] = out[path][date - 1] +
					(theta(timeline[date]) - a * out[path][date - 1]) * dt +
					sigma * sqrt(dt) * invNormalCdf(U);
			}
		}
	}
	return out;
}



template<typename T>
HullWhiteModel<double> as_double(HullWhiteModel<T>& self) {
	HullWhiteModel<double> out;
	out.r0 = as_double(self.r0);
	out.a = as_double(self.a);
	out.sigma = as_double(self.sigma);
	out.fM = as_double(self.fM);
	out.fM_spline_obj = as_double(self.fM_spline_obj);
	out.zeroCurve_obj = as_double(self.zeroCurve_obj);
	out.pM = as_double(self.pM);
	out.mats = as_double(self.mats);
	return out;
}

template<typename T>
HullWhiteModel<double> as_double(HullWhiteModel<T>* self, HullWhiteModel<double>* newself) {
	newself->r0 = as_double(self->r0);
	newself->a = as_double(self->a);
	newself->sigma = as_double(self->sigma);
	newself->fM = as_double(self->fM);
	newself->fM_spline_obj = as_double(self->fM_spline_obj);
	newself->zeroCurve_obj = as_double(self->zeroCurve_obj);
	newself->pM = as_double(self->pM);
	newself->mats = as_double(self->mats);
	return *newself;
}

