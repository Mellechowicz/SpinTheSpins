#ifndef BILLVECTOR_H
#define BILLVECTOR_H

#include <iostream>
#include <vector>
#include <initializer_list>
#include <stdexcept>
#include <cmath>
#include <type_traits>

namespace rotation{
    // based on https://github.com/Mellechowicz/bill/blob/master/headers/billVector.h
    //     _                       _                                _                      _  _  
    //    | |__  __ __ __ __ _    | |_     ___      _ _   _ _      (_)     ___    _ _     | || | 
    //    | / /  \ V  V // _` |   |  _|   / -_)    | '_| | ' \     | |    / _ \  | ' \     \_, | 
    //    |_\_\   \_/\_/ \__,_|   _\__|   \___|   _|_|_  |_||_|   _|_|_   \___/  |_||_|   _|__/  
    //  _|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_| """"| 
    //  "`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-' 

    template<class Vector, typename T = double>
        class Quaternion{
            // x*i + y*j + z*k + w <=> [s,v] s=w, v=(x,y,z)
            protected:
                T s;
                Vector v;
            public:
                Quaternion():s(1.0),v(){}
                Quaternion(const T&  _s, const Vector&  _v):s(_s),v(_v){}
                Quaternion(const T&& _s, const Vector&& _v):s(_s),v(_v){}
                Quaternion(T _s, std::initializer_list<T> list):s(_s),v(list){}
                Quaternion(Vector&  versor, T&  theta     ):s(cos(theta/2.0)),v(versor.versor()*sin(theta/2.0)){}
                Quaternion(Vector&& versor, T&& theta=M_PI):s(cos(theta/2.0)),v(versor.versor()*sin(theta/2.0)){}
                Quaternion(Vector&  versor, T&& theta=M_PI):s(cos(theta/2.0)),v(versor.versor()*sin(theta/2.0)){}
                Quaternion(const T& x,const T& y,const T& z,const T&& theta):s(cos(theta/2.0)),v(Vector(x,y,z).versor()*sin(theta/2.0)){}

                // Operators
                T  operator[](int n) const{ //element value
                    if (n<0)
                        return s;
                    else
                        return v[n]; 
                }
                T& operator[](int n){ //access elemen
                    if (n<0)
                        return s;
                    else
                        return v[n]; 
                }

                // Arithmetics
                Quaternion operator+(const Quaternion & right) const { // sum
                    return Quaternion(this->s + right.s, this->v + right.v);
                }
                Quaternion operator-(const Quaternion & right) const { // difference 
                    return Quaternion(this->s - right.s, this->v - right.v);
                }
                Quaternion operator*(const Quaternion & right) const { // multiplication (nontrivial)
                    return Quaternion( this->s*right.s - this->v*right.v,
                            (this->s*right.v) + (this->v*right.s) 
                            + (this->v^right.v));
                }

                Quaternion operator/(const T & scalar) const { // division by scalar
                    return Quaternion(s/scalar,v/scalar);
                }
                Quaternion operator*(const T & scalar) const { // multiplication by scalar rhs
                    return Quaternion(s*scalar,v*scalar);
                }
                friend Quaternion operator*(const T & scalar, const Quaternion& right) {
                    // multiplication by scalar lhs
                    return right*scalar;
                }

                Quaternion operator-() const{ // reverse
                    return Quaternion(-s,-v);
                }

                void operator/=(const T & scalar){ // self division
                    s /= scalar;
                    v /= scalar;
                }

                void normalize(){ // normalization
                    T n = norm(*this);
                    s /= n;
                    v /= n;
                }

                Quaternion operator*() const{    // conjugate Quaternion
                    return Quaternion(this->s,-(this->v));
                }

                Quaternion operator!() const{    // reciprocal Quaternion
                    return *(*this)/Quaternion::square_norm(*this);
                }

                static T square_norm(Quaternion q) {
                    return q.s*q.s + Vector::square_norm(q.v);
                }
                static T norm(Quaternion q) {
                    return sqrt(square_norm(q));
                }
                static void rotate(Vector& vec, const Vector& axis, const T& angle){
                    Quaternion q (axis,angle);
                    vec = (q*Quaternion(0.0,vec)*(!q)).get_vector();
                }
                void   rotateMe(Vector& vec) const{
                    vec = ((*this)*Quaternion(0.0,vec)*(!(*this))).get_vector();
                }
                Vector rotate(const Vector& vec) const{
                    return ((*this)*Quaternion(0.0,vec)*(!(*this))).get_vector();
                }
                void set_vector(const Vector&  _v){
                    v = _v;
                }
                void set_vector(const Vector&& _v){
                    v = _v;
                }
                Vector get_vector(){
                    return v;
                }
                void set_scalar(const T&  _s){
                    s = _s;
                }
                void set_scalar(const T&& _s){
                    s = _s;
                }
                T get_scalar(){
                    return s;
                }
                friend std::ostream& operator<<(std::ostream& stream, const Quaternion& q){
                    stream<<q.s<<" | "<<q.v;
                    return stream;
                }
        }; // end class Quaternion

} // end of namespace rotation
#endif
