#pragma once

template<typename T>
class Timeline
{
	// Simple extension of vector only used to generate a time-grid
	// Clearly not necessary but can be coded much more general to encompass the merging
	// of two time-grids
public:
	std::vector<T> line;

	Timeline() {};
	Timeline(T from, T to, T by) {
		T d = from;
		while (d < to) {
			line.push_back(d);
			d = d + by;
		}
	};
	~Timeline() {};

	size_t size() { return line.size(); }


	T& operator[] (const size_t i) { return line[i]; };
	const T& operator[] (const size_t i) const { return line[i]; }

};

template<typename T>
std::vector<T> remaining_cvg(T current_time, Timeline<T> timeline) {

	std::vector<T> out;

	int it = 0;
	double TOL = 0.00001;
	while (current_time > timeline[it] || abs(current_time - timeline[it]) < TOL) {
		it++;
	}
	out.push_back(timeline[it] - current_time);
	it++;

	while (it < timeline.line.size()) {
		out.push_back(timeline[it] - timeline[it - 1]);
		it++;
	}

	return out;
}

template<typename T>
std::vector<T> cumsum(std::vector<T> vec) {
	std::vector<T>  out(vec.size());
	out[0] = vec[0];
	for (size_t i = 1; i < vec.size(); i++)
	{
		out[i] = out[i - 1] + vec[i];
	}
	return out;
}

template<typename T>
Timeline<double> as_double(Timeline<T>& timeline) {
	Timeline<double> out;
	out.line = as_double(timeline.line);
	return out;
}
