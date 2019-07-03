#ifndef _AUX_H
#define _AUX_H
#include <bitset>
#include <limits>
#include <tuple>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

namespace aux{
template<typename VectorType>
std::pair<double,double> map_geometry(const std::vector<std::tuple<size_t,VectorType,double>>& flippable, double& decayCoeff){
    double minDistance = std::numeric_limits<double>::max();
    double maxDistance = 0.0;
    for(const auto& atom1 : flippable){
        for(const auto& atom2 : flippable){
             auto d = VectorType::norm(std::get<1>(atom1) - std::get<1>(atom2));
            if (d > 1e-7 && d < minDistance) minDistance = d;
            if (d > maxDistance)             maxDistance = d;
        }
    }
    decayCoeff = 24*log10(2)*log(10)/(maxDistance-minDistance);
    return std::make_pair(minDistance,maxDistance);
}

template<typename VectorType>
struct System{
    std::array<VectorType,3> basis;
    std::vector<std::tuple<size_t,VectorType,double>> supercell; 
    std::vector<std::tuple<size_t,VectorType,double>> flippable; 
};
    
template<typename VectorType>
std::vector<std::tuple<unsigned,unsigned,double>> get_interactions(System<VectorType>& system, double decayCoeff){
    std::vector<std::tuple<unsigned,unsigned,double>> d;
    for (const auto& atom1 : system.flippable) {
        for (const auto& atom2 : system.supercell) {
            if(std::fabs(std::get<2>(atom2)) > 1e-9 && std::get<0>(atom1) != std::get<0>(atom2)) 
                d.push_back(std::make_tuple(std::get<0>(atom1),
                                            std::get<0>(atom2),
                                            -exp(decayCoeff*VectorType::norm(
                                                std::get<1>(atom1) - std::get<1>(atom2) ))));
        }
    }
    return d;
}

template<typename State,size_t SITESNUMBER>
void print_state(const State& state, size_t reference, std::bitset<SITESNUMBER> mask){
    for (unsigned i=0; i<reference; ++i){
		if(mask[i])  std::cout<<"\033[32m"<<state[i]<<"\033[39m";
		else         std::cout<<state[i];
    }
    std::cout<<"\033[1m"<<state[reference]<<"\033[0m";
    for (unsigned i=reference+1; i<SITESNUMBER; ++i){
		if(mask[i])  std::cout<<"\033[91m"<<state[i]<<"\033[39m";
		else         std::cout<<state[i];
    }
	    std::cout<<std::endl;
}

template<typename VectorType>
void read_basis(std::array<VectorType,3>& basis,
                char _basisFile[]){ 
    std::size_t pos,buff;
    std::string line;
    std::ifstream directions(_basisFile);
    if(directions.is_open()) {
        size_t id = 0;
        while (getline(directions,line)){
            buff = 0;
            std::array<double,3> position;
            for(int i=0; i<3; ++i){
              position[i] = std::stod(line.substr(buff),&pos);
              buff += pos;
            }
            basis[id++] = position;
        }
        directions.close();
    }
}

template<typename VectorType>
void read_supercell(std::vector<std::tuple<size_t,VectorType,double>>& supercell, char _supercellFile[]){ 
    std::size_t pos,buff;
    std::string line;
    std::ifstream inSupercell(_supercellFile);
    if(inSupercell.is_open()) {
        while (getline(inSupercell,line)){
            size_t idx = std::stoi(line,&pos);
            buff = pos;
            std::array<double,3> position;
            for(int i=0; i<3; ++i){
              position[i] = std::stod(line.substr(buff),&pos);
              buff += pos;
            }
            auto moment = std::stod(line.substr(buff));
            supercell.push_back(std::make_tuple(idx,position,moment));
        }
        inSupercell.close();
    }
}

template<typename VectorType>
void read_flippables(std::vector<std::tuple<size_t,VectorType,double>>& flippable, char _flippableFile[]){ 
    std::size_t pos,buff;
    std::string line;
    std::ifstream inFlippable(_flippableFile);
    if(inFlippable.is_open()) {
        while (getline(inFlippable,line)){
            size_t idx = std::stoi(line,&pos);
            buff = pos;
            std::array<double,3> position;
            for(int i=0; i<3; ++i){
              position[i] = std::stod(line.substr(buff),&pos);
              buff += pos;
            }
            auto moment = std::stod(line.substr(buff));
            flippable.push_back(std::make_tuple(idx,position,moment));
        }
        inFlippable.close();
    }
}
void greetings();
} //end of namespace aux
#endif
