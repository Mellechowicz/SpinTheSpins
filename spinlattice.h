#ifndef _SPINLATTICE_H
#define _SPINLATTICE_H

#include <iostream>

#include <tuple>
#include <array>
#include <vector>
#include <bitset>

#include <algorithm>
#include <random>
#include <string>

#include <memory>
#include <limits>
#include <type_traits>
#include <algorithm>

#include "asa.h"
#include "arithmeticvector.h"
#include "quaternion.h"

namespace metropolis{

typedef celerium::ArithmeticVector   Vector;
typedef rotation::Quaternion<Vector> Quaternion;

template<size_t N> double spinEnergy (void*  state);
template<size_t N> double spinMeasure(void*  stateI,
                                      void*  stateJ);
template<size_t N> void   spinStep   (const  gsl_rng* random,
                                      void*  state,
                                      double step);
template<size_t N> void   spinPrint  (void*  state);

template<size_t N>
struct SpinLattice{
    SpinLattice(){};
    SpinLattice(const SpinLattice& rhs):nodes(rhs.nodes){};
    SpinLattice(const std::array<Vector,N>& _nodes){
         std::copy(_nodes.begin(),_nodes.end(),nodes.begin());
    }
    SpinLattice& operator=(const SpinLattice& rhs){
        nodes = rhs.nodes;
        return *this;
    }
    void set_nodes(const std::array<Vector,N>& _nodes){
         std::copy(_nodes.begin(),_nodes.end(),nodes.begin());
    }
    std::array<Vector,N> nodes;
    friend std::ostream& operator<<(std::ostream& stream, const SpinLattice& lattice){
        for (const auto& node : lattice.nodes)
            stream<<node<<std::endl;
        return stream;
    }
};

template<size_t N>
class Solver{
public:
    Solver(const Solver& rhs):Solver(){
        this->lattice = rhs.get_nodes();
    }

    Solver(double _T = 273.0, double average = 0.0, double sigma = 0.1):
        randomDevice(),
        randomEngine(randomDevice()),
        gauss(average,sigma),uniform(-1.0,1.0),all_angles(-M_PI,M_PI),T(_T){
        hamiltonian.reserve(N*N/2); // maximum number of interactions ( each site with each )

//            solver.set_energy (  isingEnergy<N> );
//            solver.set_measure( isingMeasure<N> );
//            solver.set_step   (    isingStep<N> );
#ifdef _VERBOSE
//            solver.set_print  (   isingPrint<N> );
#endif
        }

    Solver& operator=(const Solver& rhs){
        lattice = rhs.get_nodes();
        return *this;
    }
    ~Solver(){}

    typedef std::tuple<unsigned,unsigned,double>  TwoSiteInteraction;
    typedef std::vector<TwoSiteInteraction>       HamiltonianType;

protected:
    std::random_device                              randomDevice;
    std::mt19937                   randomEngine;
    std::normal_distribution<>                      gauss;
    std::uniform_real_distribution<>                uniform;
    std::uniform_real_distribution<>::param_type    all_angles;

    double T;

    SpinLattice<N>                                  lattice;

//    gsl::SimulatedAnnealing                         solver;
    HamiltonianType                                 hamiltonian;

public:
//    typename gsl::SimulatedAnnealing::Parameters&
//    set_parameters(const typename gsl::SimulatedAnnealing::Parameters& _params){
//            return solver.set_parameters(_params);
//    }

protected:
    // For SFINAE compiler-time evaulation
    template<class T>
    T tester(T t)const{
        if(std::is_integral<T>::value) return static_cast<unsigned>(t);
        return t;
    }

public:
    void add_interaction(...){
        std::cerr<<"Wrong input for auxiliary::IsingModel::add_interaction:"<<std::endl;
        std::cerr<<"       Either non-iterable or iterable of non <int,int,float> tuples."<<std::endl;
        std::cerr<<"       Nothing was added!"<<std::endl;
    }

    template<class intlike, class floatlike>
    auto add_interaction(intlike i, intlike j, floatlike J) -> decltype((unsigned)(tester<intlike>)(i),void()){
        hamiltonian.push_back(std::make_tuple(i,j,J));
    }

    template<class T>
    auto add_interaction(T interaction) -> decltype((TwoSiteInteraction&&)(tester<T>)(interaction),void()){
        hamiltonian.push_back(interaction);
    }

    template<class T>
    auto add_interaction(T interaction) -> decltype((TwoSiteInteraction&)(tester<T>)(interaction),void()){
        hamiltonian.push_back(interaction);
    }

    template<class Iterable>
    auto add_interaction(const Iterable& interactions) -> decltype((decltype(interactions.begin()))(std::begin)(interactions),void()){
        for(auto& interaction : interactions)
            hamiltonian.push_back(interaction);
    }

    const SpinLattice<N>& get_nodes()const{
        return lattice;
    }

    SpinLattice<N>& get_nodes(){
        return lattice;
    }

    SpinLattice<N>* get_nodes_ptr(){
        return &lattice;
    }
    
    void set_nodes(const std::array<Vector,N>& _nodes){
         lattice.set_nodes(_nodes); 
    }

    void clear_hamiltonian(){
        hamiltonian.clear();
    }

    void reset(){
        this->clear_hamiltonian();
        this->randomize_state();
    }

    void randomize_state(){
        for (auto& node : lattice.nodes){
            Quaternion::rotate(node,uniform(randomEngine),uniform(randomEngine),uniform(randomEngine),gauss(randomEngine));
        }
    }

    void random_state(){
        for (auto& node : lattice.nodes){
            Quaternion::rotate(node,uniform(randomEngine),uniform(randomEngine),uniform(randomEngine),uniform(randomEngine,all_angles));
        }
    }

    const HamiltonianType& get_hamiltonian() const{
        return hamiltonian;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Solver& model){
        return stream<<model.get_nodes();
    }
    
    double energy(){
        double E = 0.0;
        for(const auto& interaction : hamiltonian)
            E += std::get<2>(interaction)*(lattice.nodes[std::get<0>(interaction)]*lattice.nodes[std::get<1>(interaction)]);
        return E;
    }
/*
    static double measure(const std::bitset<N>& stateI, const std::bitset<N>& stateJ){
        std::bitset<N> output = ~stateI & stateJ;
        return output.count();
    }    

    std::bitset<N> run(){
        this->mask = nullptr;
#ifdef _VERBOSE
        std::cout<<"Starting from: "<<lattice.nodes<<std::endl;
#endif
        solver.run<LatticeType<N,IsingModel<N>>>(lattice,sizeof(lattice));
#ifdef _VERBOSE
        std::cout<<"Solution:      ";
        std::cout<<lattice.nodes<<std::endl;
#endif
        return lattice.nodes;
    }

#ifdef _VERBOSE
    std::bitset<N> run(std::bitset<N>* mask){
        this->mask = mask;
        lattice.nodes &= *mask;
        std::cout<<"Starting from: "<<lattice.nodes<<std::endl;
        solver.run<LatticeType<N,IsingModel<N>>>(lattice,sizeof(lattice));
        std::cout<<"Solution:      ";
        std::cout<<lattice.nodes<<std::endl;
        return lattice.nodes;
    }
#else
    std::bitset<N> run(std::bitset<N>* mask){
        this->mask = mask;
        lattice.nodes &= *mask;
        solver.run<LatticeType<N,IsingModel<N>>>(lattice,sizeof(lattice));
        return lattice.nodes;
    }
#endif

    static void generate_state(std::bitset<N>& state,
                               std::mt19937& engine,
                               std::uniform_int_distribution
                                    <unsigned long long int>& distribution){
        constexpr auto seedSize = 8*sizeof(unsigned long long int);

        state = std::bitset<N>(distribution(engine));
        auto currentSize = seedSize;

        while (currentSize < N){
            state <<= seedSize;
            state  |= std::bitset<N>(distribution(engine));
            currentSize += seedSize;
        }
    }*/
}; // end of class IsingModel
/*
template<size_t N>
double isingEnergy (void* state){
    return IsingModel<N>::energy(static_cast<LatticeType<N,IsingModel<N>>*>(state));
}

template<size_t N>
double isingMeasure(void* stateI, void* stateJ){
    return IsingModel<N>::measure(
            static_cast<LatticeType<N,IsingModel<N>>*>(stateI)->nodes,
            static_cast<LatticeType<N,IsingModel<N>>*>(stateJ)->nodes
           );
}

template<size_t N>
void   isingStep   (const gsl_rng* random __attribute__((unused)), void* state, double step __attribute__((unused))){
    IsingModel<N>::randomize(
            static_cast<LatticeType<N,IsingModel<N>>*>(state)->nodes,
            static_cast<LatticeType<N,IsingModel<N>>*>(state)->randomEngine,
            static_cast<LatticeType<N,IsingModel<N>>*>(state)->uniform,
            static_cast<LatticeType<N,IsingModel<N>>*>(state)->model->get_mask()
           );
}

template<size_t N>
void   isingPrint  (void* state){
    std::cout<<'\t'<<static_cast<LatticeType<N,IsingModel<N>>*>(state)->nodes;
}
*/
} // end namespace metropolis
#endif
