#pragma once

template <typename T>
class VasicekModel : public BaseRateModel<T>
{
	// Vasicek short rate model object.
	// Most imporant are the methods which returns forward rates, bond prices and simulate short rates
public:
	VasicekModel() {};
	VasicekModel(T myMu, T myBeta, T myVol, T r0) :
		mu(myMu), beta(myBeta), vol(myVol), r0(r0) {};
	~VasicekModel() {};

	T bond_price(T time, T mat, T& shortrate) {
		T timetoexp = mat - time;
		T B = 1 / beta * (1 - exp(-beta*timetoexp));
		T A = (B - timetoexp) * (mu - vol*vol / (2.0*beta*beta)) - B * B * vol * vol / (4 * beta);
		return(exp(A - B*shortrate));
	}

	T yield(T time, T mat, T& shortrate) {
		return(-log(bond_price(time, mat, shortrate)) / (mat - time));
	}


	T spot_libor(T current_time, T maturity, T& short_rate) {
		T ptT = bond_price(current_time, maturity, short_rate);
		return (1 - ptT) / ((maturity - current_time)*ptT);
	}

	T disc_factor(T start, T end, T& short_rate) {
		return  bond_price(start, end, short_rate);
	}

	T forward_rate(T start, T maturity, T& short_rate) {
		T df_start;
		double TOL = 0.000001;
		if (abs(start - 0.0) < TOL) {
			df_start = 1;
		}
		else {
			df_start = disc_factor(0, start, short_rate);
		}

		T df_mat = disc_factor(0, maturity, short_rate);

		T out = (df_start / df_mat - 1) / (maturity - start);
		return out;
	}

	Matrix<T> simulate(Timeline<T>& timeline, const size_t paths, const T& r0, mrg32& rng, bool byrow = false) {
		// Returns a matrix of simulated short rate paths from the Hull-White model.
		// In the returned matrix each row represent a path

		size_t steps = timeline.line.size();
		Matrix<T> out(paths, steps);

		for (size_t i = 0; i < paths; i++)
		{
			out[i][0] = r0;
		}

		for (size_t date = 1; date < steps; date++)
		{
			T dt = timeline[date] - timeline[date - 1];
			for (size_t path = 0; path < paths; path++)
			{
				double U = rng.next();
				out[path][date] = out[path][date - 1] +
					beta * (mu - out[path][date - 1]) * dt +
					vol * sqrt(dt) * invNormalCdf(U);
			}
		}
		return out;
	}

	T r0;
	T mu;
	T beta;
	T vol;
};


template<typename T>
VasicekModel<double> as_double(VasicekModel<T>& mdl) {
	VasicekModel<double> out(value(mdl.mu), value(mdl.beta),
		value(mdl.vol), value(mdl.r0));
	return out;
}
