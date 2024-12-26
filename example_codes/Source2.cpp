#include"Header2.h"

int main() {
	finite_difference(1.0, 3.0, 0.1,[](const double& x)->const double
		{ return std::exp(2 * x); }
	);
}
