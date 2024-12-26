#include"Header5.h"

int main() {

    std::cout << simpson(-1.0, 1.0, 0.01,[](const double &x){  return std::exp(2 * x); });


}
