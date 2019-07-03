#include <iostream>
#include "../../arithmeticvector.h"
#include "../../quaternion.h"

int main(){
    typedef celerium::ArithmeticVector Vector;
    Vector A(1.0,2.0,3.0);
    Vector B(-2.0,-3.0,-1.0);
    std::cout<<"  A = "<< A<<std::endl;
    std::cout<<" -A = "<<-A<<std::endl;
    std::cout<<"  B = "<< B<<std::endl;
    std::cout<<" -B = "<<-B<<std::endl;
    std::cout<<"A+B = "<<A+B<<std::endl;
    std::cout<<"A-B = "<<A-B<<std::endl;
    std::cout<<"A*B = "<<A*B<<std::endl;

    std::cout<<"\n -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n\n"; 

    rotation::Quaternion<Vector> qA(A);
    std::cout<<"  qA= "<< qA<<std::endl;
    std::cout<<" -qA= "<<-qA<<std::endl;
    std::cout<<" *qA= "<<*qA<<std::endl;
    std::cout<<" !qA= "<<!qA<<std::endl;

    std::cout<<"\n -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n\n"; 

    rotation::Quaternion<Vector> q(Vector(2.0,0.0,0.0),M_PI_2);
    std::cout<<"  q = "<< q <<std::endl;
    std::cout<<" -q = "<<-q <<std::endl;
    std::cout<<" *q = "<<*q <<std::endl;
    std::cout<<" !q = "<<!q <<std::endl;

    std::cout<<"\n -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n\n"; 

    rotation::Quaternion<Vector> r(1,0,1,M_PI_2);
    std::cout<<"  r = "<< r <<std::endl;
    std::cout<<" -r = "<<-r <<std::endl;
    std::cout<<" *r = "<<*r <<std::endl;
    std::cout<<" !r = "<<!r <<std::endl;

    std::cout<<"\n -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n\n"; 

    q.rotateMe(A);
    std::cout<<"rot_ q(A) = "<<A<<std::endl;
    (!q).rotateMe(A);
    std::cout<<"rot_!q(A) = "<<A<<std::endl;
    q.rotateMe(B);
    std::cout<<"rot_ q(B) = "<<B<<std::endl;
    (!q).rotateMe(B);
    std::cout<<"rot_!q(B) = "<<B<<std::endl;

    std::cout<<"\n -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n\n"; 

    r.rotateMe(A);
    std::cout<<"rot_ r(A) = "<<A<<std::endl;
    (!r).rotateMe(A);
    std::cout<<"rot_!r(A) = "<<A<<std::endl;
    r.rotateMe(B);
    std::cout<<"rot_ r(B) = "<<B<<std::endl;
    (!r).rotateMe(B);
    std::cout<<"rot_!r(B) = "<<B<<std::endl;

    return 0;
}
