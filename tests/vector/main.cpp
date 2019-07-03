#include <iostream>
#include "../../arithmeticvector.h"

int main(){
    celerium::ArithmeticVector A(1.0,2.0,3.0);
    celerium::ArithmeticVector B(-2.0,-3.0,-1.0);
    std::cout<<"  A = "<< A<<std::endl;
    std::cout<<" -A = "<<-A<<std::endl;
    std::cout<<"  B = "<< B<<std::endl;
    std::cout<<" -B = "<<-B<<std::endl;
    std::cout<<"A+B = "<<A+B<<std::endl;
    std::cout<<"A-B = "<<A-B<<std::endl;
    std::cout<<"A*B = "<<A*B<<std::endl;

    return 0;
}
