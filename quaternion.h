//
//  quaternion.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 17/1/22.
//

#ifndef quaternion_h
#define quaternion_h

#include <cmath>
#include "vector.h"

template<typename T> //T: specifies the arithmetic type used, e.g. float, double...
class quaternion{
public:
    quaternion();
    quaternion(T angle, vector<T> unitAxis); /*angle is the angle of rotation in degrees, unitAxis is the unit vector
                                              which represents the axis of rotation.*/
    
    T l2norm();
    quaternion<T> conj();
    quaternion<T> logQuat();
    void setQuaternion(T angle, vector<T> unitAxis); /*angle is the angle of rotation in rads, unitAxis is the unit vector
                                                      which represents the axis of rotation.*/
    void set(int i, T t);
    void setToZero();
    T get(int i) const;
    void vector2DtoQuaternion(const vector<T> vec);
    void vector3DtoQuaternion(const vector<T> vec);
    
    quaternion<T> &operator= (const quaternion<T> &right);
    quaternion<T> operator+ (const quaternion<T> &right) const;
    quaternion<T> operator* (const quaternion<T> &right) const;
    template<typename U> friend quaternion<U> operator* (const U &left, quaternion<U> &right);
    quaternion<T> operator* (const T &right) const;
    quaternion<T> operator/ (const T &right) const;
    quaternion<T> operator+= (const quaternion<T> &right);
private:
    T q[4];
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Definitions of functions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
inline quaternion<T>::quaternion(){
    for (int i=0; i<4; i++)
        q[i]=0.;
}

template<typename T>
quaternion<T>::quaternion(T angle, vector<T> unitAxis){
    double sinAngle= sin(angle/2.);
    
    q[0]= cos(angle/2.);
    for (int i=1; i<4; i++)
        q[i]= sinAngle*unitAxis(i-1);
}

template<typename T>
void quaternion<T>::setQuaternion(T angle, vector<T> unitAxis){
    double sinAngle= sin(angle/2.);
    
    q[0]= cos(angle/2.);
    for (int i=1; i<4; i++)
        q[i]= sinAngle*unitAxis(i-1);
    
}

template<typename T>
inline T quaternion<T>::l2norm(){
    T magnitude= sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    return magnitude;
}

template<typename T>
inline void quaternion<T>::set(int i, T t){
    q[i]= t;
}

template<typename T>
inline void quaternion<T>::setToZero(){
    for (int i=0; i<4; i++)
        q[i]= 0.;
}

template<typename T>
inline T quaternion<T>::get(int i) const{
    return q[i];
}

template<typename T>
quaternion<T> quaternion<T>::conj(){
    quaternion<T> conjugate;
    conjugate.q[0]=q[0];
    for (int i=1; i<4; i++)
        conjugate.q[i]= -q[i];
    
    return conjugate;
}

template<typename T>
quaternion<T> quaternion<T>::logQuat(){
    quaternion<T> logQ;
    int angle= 2.*acos(q[0]);
    T sinHAngle= sin(angle/2);
    vector<T> unitAxis(3);
    
    unitAxis.set(0,q[1]/sinHAngle);
    unitAxis.set(1,q[2]/sinHAngle);
    unitAxis.set(2,q[3]/sinHAngle);
    
    logQ.q[0]= 0.;
    
    logQ.q[1]= angle/2.*unitAxis(0);
    logQ.q[2]= angle/2.*unitAxis(1);
    logQ.q[3]= angle/2.*unitAxis(2);
    
    return logQ;
}

template<typename T>
inline void quaternion<T>::vector2DtoQuaternion(const vector<T> vec){
    q[0]=0.;
    q[1]=vec(0);
    q[2]=vec(1);
    q[3]=0.;
}

template<typename T>
inline void quaternion<T>::vector3DtoQuaternion(const vector<T> vec){
    q[0]=0.;
    q[1]=vec(0);
    q[2]=vec(1);
    q[3]=vec(2);
}

template<typename T>
quaternion<T> quaternion<T>::operator+ (const quaternion<T> &right) const{
    quaternion<T> result;
    
    result.q[0]= q[0] + right.q[0];
    result.q[1]= q[1] + right.q[1];
    result.q[2]= q[2] + right.q[2];
    result.q[3]= q[3] + right.q[3];
    
    return result;
}

template<typename T>
quaternion<T> quaternion<T>::operator* (const quaternion<T> &right) const{
    quaternion<T> result;
    result.q[0]= q[0]*right.q[0]- q[1]*right.q[1] - q[2]*right.q[2] - q[3]*right.q[3];
    
    result.q[1]= q[0]*right.q[1] + q[1]*right.q[0] + q[2]*right.q[3] - q[3]*right.q[2];
    result.q[2]= q[0]*right.q[2] + q[2]*right.q[0] + q[3]*right.q[1] - q[1]*right.q[3];
    result.q[3]= q[0]*right.q[3] + q[3]*right.q[0] + q[1]*right.q[2] - q[2]*right.q[1];
    
    return result;
}

template <typename T>
quaternion<T> &quaternion<T>::operator= (const quaternion<T> &right)
{
    if (this==&right)
        return *this;
    for (int i=0;i<4;i++)
        q[i]= right.q[i];
    
    return *this;
}

template <typename T>
quaternion<T> quaternion<T>::operator* (const T &right) const
{
    quaternion<T> result= *this;
    for (int i=0; i<4; i++)
        result.q[i]*=right;

    return result;
}

template<typename T>
quaternion<T> operator* (const T &left, quaternion<T> &right)
{
    quaternion<T> result=right;
    for (int i=0; i<4; i++)
        result.set(i,result.get(i)*left);

    return result;
}

template <typename T>
quaternion<T> quaternion<T>::operator/ (const T &right) const
{
    quaternion<T> result= *this;
    for (int i=0; i<4; i++)
        result.set(i,result.get(i)/right);

    return result;
}

template <typename T>
inline quaternion<T> quaternion<T>::operator+= (const quaternion<T> &right){
    for (int i=0; i<4;i++)
        q[i]+= right.q[i];
    
    return *this;
}

#endif /* quaternion_h */
