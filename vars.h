#pragma once


#include <vector>
#include <math.h>

extern Tape *GlobalTape;

class Variable
{
	// Main AAD object class
	// Needs to have all possible operations overloaded to be used as a value type
	
public:
	Variable() {}
	Variable(double val)
		: tape(GlobalTape), value(val) {
		idx = tape->push0();
	}

	Variable(double val, size_t index)
		: tape(GlobalTape), value(val), idx(index) {}

	std::vector<double> calc_derivatives() {
		// Returns a vector of the adjoints wrt. to the variable it is called on
		size_t len = tape->nodes.size();
		std::vector<double> derivs(len);
		derivs[len - 1] = 1.0; //With respect to itself.

		for (int i = len - 1; i > 0; i--)
		{
			Node node = tape->nodes[i];
			double deriv = derivs[i];

			derivs[node.deps[0]] += node.weights[0] * deriv;
			derivs[node.deps[1]] += node.weights[1] * deriv;
		}
		return derivs;
	}

	Tape* tape;
	size_t idx;
	double value;
	
	Variable& operator+=(const Variable& rhs) {
		idx = tape->push2(idx, 1.0, rhs.idx, 1.0);
		value += rhs.value;
		return *this;
	}

	Variable operator+(const Variable& other) {
		return Variable(*this) += other;
	}
	
	Variable& operator-=(const Variable& rhs) {
		idx = tape->push2(idx, 1.0, rhs.idx, -1.0);
		value -= rhs.value;
		return *this;
	}

	Variable operator-(const Variable& other) {
		return Variable(*this) -= other;
	}

	Variable& operator*=(const Variable& rhs) {
		idx = tape->push2(idx, rhs.value, rhs.idx, value);
		value *= rhs.value;
		return *this;
	}

	Variable operator*(const Variable& other) {
		return Variable(*this) *= other;
	}

	Variable& operator/=(const Variable& rhs) {
		idx = tape->push2(idx, 1.0/rhs.value, rhs.idx, -value/(rhs.value*rhs.value));
		value /= rhs.value;
		return *this;
	}

	 Variable operator/(const Variable& other) {
		return Variable(*this) /= other;
	}
	
}; ///////// End of Variable class


double wrt(Variable& var, std::vector<double>& Derivatives) {
	// Picks out the derivative from the derivatives vector
	return Derivatives[var.idx];
}

   // Operators en masse

inline Variable operator+(Variable& lhs, double rhs) {
	return Variable(lhs.value + rhs,
		lhs.tape->push1(lhs.idx, 1.0));
}

inline Variable operator+(double lhs, Variable& rhs) {
	return Variable(rhs.value + lhs,
		rhs.tape->push1(rhs.idx, 1.0));
}


inline Variable operator-(Variable& var) {
	return Variable(-var.value,
		var.tape->push1(var.idx, -1.0));
}


inline Variable operator-(Variable& lhs, double rhs) {
	return Variable(lhs.value - rhs,
		lhs.tape->push1(lhs.idx, 1.0));
}

inline Variable operator-(double lhs, Variable& rhs) {
	return Variable(lhs - rhs.value,
		rhs.tape->push1(rhs.idx, -1.0));
}

inline Variable operator*(Variable& lhs, double rhs) {
	return Variable(lhs.value * rhs,
		lhs.tape->push1(lhs.idx, rhs));
}

inline Variable operator*(double lhs, Variable& rhs) {
	return Variable(rhs.value * lhs,
		rhs.tape->push1(rhs.idx, lhs));
}

inline Variable operator/(Variable& lhs, double rhs) {
	return lhs * (1.0 / rhs);
}

inline Variable operator/(double lhs, Variable& rhs) {
	return Variable(lhs / rhs.value,
		rhs.tape->push1(rhs.idx, -lhs / (rhs.value*rhs.value)));
}

inline bool operator>(const Variable& lhs, const Variable& rhs) {
	return lhs.value > rhs.value ? true : false;
}

inline bool operator>=(const Variable& lhs, const Variable& rhs) {
	return lhs.value >= rhs.value ? true : false;
}

inline bool operator<(const Variable& lhs, const  Variable& rhs) {
	return lhs.value < rhs.value ? true : false;
}

inline bool operator<=(const Variable& lhs, const  Variable& rhs) {
	return lhs.value <= rhs.value ? true : false;
}

inline bool operator>(const Variable& lhs, const double& rhs) {
	return lhs.value > rhs ? true : false;
}

inline bool operator>=(const Variable& lhs, const double& rhs) {
	return lhs.value >= rhs ? true : false;
}

inline bool operator<(const Variable& lhs, const double& rhs) {
	return lhs.value < rhs ? true : false;
}

inline bool operator<=(const Variable& lhs, const double& rhs) {
	return lhs.value <= rhs ? true : false;
}

inline bool operator>(const double& lhs, const  Variable& rhs) {
	return lhs > rhs.value ? true : false;
}

inline bool operator>=(const double& lhs, const Variable& rhs) {
	return lhs >= rhs.value ? true : false;
}

inline bool operator<(const double& lhs, const Variable& rhs) {
	return lhs < rhs.value ? true : false;
}

inline bool operator==(const Variable& lhs, const Variable& rhs) {
	return abs(lhs.value - rhs.value) < 0.000001 ? true : false;
}


inline bool operator==(const Variable& lhs, const double& rhs) {
	return abs(lhs.value - rhs) < 0.000001 ? true : false;
}

inline bool operator==(const double& lhs, const Variable& rhs) {
	return (rhs == lhs);
}

inline bool operator!=(const Variable& lhs, const Variable& rhs) {
	return abs(lhs.value - rhs.value) < 0.000001 ? false : true;
}

inline bool operator!=(const Variable& lhs, const double& rhs) {
	return abs(lhs.value - rhs) < 0.000001 ? false : true;
}

inline bool operator!=(const double& lhs, const Variable& rhs) {
	return (rhs != lhs);
}

// Some basic funcs
inline Variable sin(Variable self) {
	return Variable(sin(self.value),
		(self.tape)->push1(self.idx, cos(self.value)));
}

inline Variable exp(Variable self) {
	return Variable(exp(self.value),
		(self.tape)->push1(self.idx, exp(self.value)));
}

inline Variable sqrt(Variable self) {
	return Variable(sqrt(self.value),
		(self.tape)->push1(self.idx, 1 / (2 * sqrt(self.value))));
}

inline Variable abs(Variable self) {
	return Variable(abs(self.value),
		(self.tape)->push1(self.idx, self.value > 0 ? 1 : -1));
}

inline Variable log(Variable self) {
	return Variable(log(self.value),
		(self.tape)->push1(self.idx, 1 / (self.value)));
}

double value(double& k) {
	return k;
}

double value(Variable& k) {
	return k.value;
}

template<typename T>
double value(T& const k) {
	return value(k);
}

template<typename T>
double as_double(T& k) {
	return value(k);
}

double as_double(double k) {
	return k;
}

double as_double(Variable k) {
	return k.value;
}