#pragma once
#include <vector>
#include <math.h>

class Node 
{
	// Represents a single point in the a computational tree. As a default we allow for 
	// two dependencies and two weights.
public:
	Node(double w0, double w1, size_t dep0, size_t dep1)
		: weights{ w0, w1 }, deps{ dep0, dep1 } {}
	double weights[2];
	size_t deps[2];
}; ///////// End of Node class

class Tape 
{
	// Vector of Nodes used to store the previous computations.
	// Set up such that a pointer to a tape has to be present in the global environment
public:
	std::vector<Node> nodes;

	size_t push0() {
		size_t len = nodes.size();
		nodes.push_back(Node(0.0, 0.0, len, len));
		return len;
	}

	size_t push1(size_t dep, double weight)	{
		size_t len = nodes.size();
		nodes.push_back(Node(weight, 0.0, dep, len));
		return len;
	}
	
	size_t push2(size_t dep0, double w0, size_t dep1, double w1) {
		size_t len = nodes.size();
		nodes.push_back(Node(w0, w1, dep0, dep1));
		return len;
	}

}; ///////// End of Tape class

