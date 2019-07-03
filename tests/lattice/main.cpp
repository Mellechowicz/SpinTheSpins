#include <iostream>
#include <array>
#include "../../spinlattice.h"
#include "../../arithmeticvector.h"

int main(int argc, char** argv){
    constexpr size_t N = 144;

    double J     = -1.0; //ferro
    double delta = 0.98; // sigma mult
    if (argc > 1) J     = std::atof(argv[1]);
    if (argc > 2) delta = std::atof(argv[2]);
    typedef celerium::ArithmeticVector Vector;
    std::normal_distribution<>::param_type p(0.0,1.0);
    std::array<Vector,N>   lattice;
    metropolis::Solver<N> spinLatticeSolver;
    for (size_t i = 0; i<N; ++i){
       lattice[i]                       = Vector(i%int(sqrt(N)),i%int(sqrt(N)),0.0)*3.0;
       spinLatticeSolver.get_nodes().nodes[i] = Vector(1.0,0.0,0.0);
    }
    for (size_t i = 0; i<N; ++i)
        for (size_t j = 0; j<i; ++j)
            if(Vector::norm(lattice[i]- lattice[j]) <= 3.1)
                spinLatticeSolver.add_interaction(i,j,J);
    spinLatticeSolver.random_state();
    auto E0 = spinLatticeSolver.energy();
    for(int i =0; i<10000; ++i){
        auto rmb = spinLatticeSolver;
        spinLatticeSolver.randomize_state2(p);
        auto Et = spinLatticeSolver.energy();
        if (Et >= E0)
           spinLatticeSolver = rmb;
        else
            E0=Et;
    }
    std::cout<<"E = "<<E0<<std::endl;
    for(int i =0; i<10000000; ++i){
        auto rmb = spinLatticeSolver;
        spinLatticeSolver.randomize_state2(p);
        auto Et = spinLatticeSolver.energy();
        if (Et >= E0)
           spinLatticeSolver = rmb;
        else{
            p = std::normal_distribution<>::param_type(0.0,delta*p.stddev());
            auto dE = Et-E0;
            std::cout<<"E_"<<i<<"= "<<Et<<"\tdE= "<<dE<<"\t"<<p.stddev()<<std::endl;
            E0=Et;
            if (p.stddev() < 1e-5 && abs(dE) < 1e-13) break;
        }
    }
    std::cout<<"E = "<<E0<<std::endl;
    auto Moment = spinLatticeSolver.mean();
    std::cout<<"M = "<<Moment<<" <M> = "<<Vector::norm(Moment)<<std::endl;

    return 0;
}

