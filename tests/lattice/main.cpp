#include <iostream>
#include <array>
#include "../../spinlattice.h"
#include "../../arithmeticvector.h"

int main(){
    typedef celerium::ArithmeticVector Vector;
    constexpr size_t N = 49;
    std::array<Vector,N>   lattice;
    metropolis::Solver<N> spinLatticeSolver;
    for (size_t i = 0; i<N; ++i){
       lattice[i]                       = Vector(i%int(sqrt(N)),i%int(sqrt(N)),0.0)*3.0;
       spinLatticeSolver.get_nodes().nodes[i] = Vector(1.0,0.0,0.0);
    }
    for (size_t i = 0; i<N; ++i)
        for (size_t j = 0; j<i; ++j)
            if(Vector::norm(lattice[i]- lattice[j]) <= 3.0)
                spinLatticeSolver.add_interaction(i,j,-1.0);
    std::cout<<spinLatticeSolver;
    auto E0 = spinLatticeSolver.energy();
    std::cout<<"E = "<<E0<<std::endl;
//    spinLatticeSolver.random_state();

    for(int i =0; i<100000; ++i){
        auto rmb = spinLatticeSolver;
        spinLatticeSolver.randomize_state();
//        std::cout<<std::endl<<spinLatticeSolver;
        auto Et = spinLatticeSolver.energy();
        if (Et >= E0)
           spinLatticeSolver = rmb;
        else
            E0=Et;
        std::cout<<"E = "<<Et<<std::endl;
    }
//    std::cout<<spinLatticeSolver;
    E0 = spinLatticeSolver.energy();
    std::cout<<"E = "<<E0<<std::endl;

    return 0;
}

