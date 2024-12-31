#include"Header5.h"

int main() {

     if(auto u=simpson(-1.0, 1.0, 0.01,[](const double &x){  return std::exp(2 * x); });u.has_value()){
      std::cout<<u.value();
    }


}
