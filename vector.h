//
//  vector.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 15/11/21.
//

#ifndef vector_h
#define vector_h

#include <cmath>

template <typename T>
class vector
{
public:
    vector(int N=0,const T &a=0); //default constructor
    vector(int N,int i, const T &a); //constructor for vector with all component zero but one which equals a
    vector(const T &x,const T &y,const T &z=0); //2-d/3-d vector constructor
    vector(const vector<T> &ob);//copy constructor
    ~vector();//destructor
    
    void setDimension(int i);
    int getDimension();
    void set(int i, const T &a) const;
    T l2norm() const;
    vector<T> cross_product(const vector<T> &a, const vector<T> &b);//works only for n<=3
    T findAngle(vector<T> vec); //finds the angle of the vector and the vector vec
    
    T &operator() (int i) const;//reads the i-th component of the vector
    vector<T> operator+ (const vector<T> &right) const;
    vector<T> operator- (const vector<T> &right) const;
    vector<T> operator/ (const T &right) const;
    vector<T> operator* (const T &right) const;
    vector<T> operator+= (const vector<T> &right);
    vector<T> operator-= (const vector<T> &right);
    vector<T> operator*= (const T &right);
    vector<T> operator/= (const T &right);
    template<typename U> friend vector<U> operator* (const U &left, vector<U> &right);
    T operator* (const vector<T> &right) const; //dot product
    vector<T> &operator= (const vector<T> &right);
    template<typename t> friend int operator== (const vector<t> &v1, const vector<t> &v2);
    template <typename U> friend std::ostream &operator<< (std::ostream &left, const vector<U> &right);
protected:
    int dimension;
    T *component;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Definitions of functions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//default constructor
template <typename T>
vector<T>::vector(int N, const T &a)
//N: dimension of vector
//a: value of components of vector
{
    dimension=N;
    if (N == 0)
    {
        component=NULL;
        return;
    }
    
    component= new T [dimension];
    for (int i=0;i<N;i++)
        component[i]=a;
}

//constructor for vector with all component zero but one which equals a
template <typename T>
vector<T>::vector(int N, int i, const T &a)
//N: dimension of vector
//i: i-th component of vector
//a: i-th value of component of vector
{
    dimension=N;
    component= new T [dimension];
    for (int j=0;j<N;j++)
        component[j]=0;
    component[i]=a;
}

//2-d/3-d vector constructor
template <typename T>
vector<T>::vector(const T &x,const T &y,const T &z)
{
    if (z == 0)
        dimension=2;
    else
        dimension=3;
    component= new T [dimension];
    
    component[0]=x;
    component[1]=y;
    if (dimension == 3)
        component[2]=z;
}

//copy constructor
template<typename T>
vector<T>::vector(const vector<T> &ob){
    component= new T[ob.dimension];
    dimension=ob.dimension;
    
    for (int i=0;i<dimension;i++)
        component[i]=ob.component[i];
}

//destructor
template<typename T>
vector<T>::~vector(){
    if (component)
    {
        delete [] component;
        component=NULL;
    }
    dimension=0;
}

template<typename T>
void vector<T>::setDimension(int i){
    if (dimension)
        delete [] component;
    
    dimension=i;
    component= new T[dimension];
    
    for (int i=0; i<dimension; i++)
        component[i]=0.;
}

template<typename T>
inline int vector<T>::getDimension(){
    return dimension;
}

template <typename T>
inline void vector<T>::set(int i, const T &a) const
//sets the i-th component of the vector equal to a
{
    component[i]=a;
}

template<typename T>
T vector<T>::l2norm() const
{
    T dot_product=0;
    for (int i=0; i<dimension;i++)
        dot_product+=component[i]*component[i];
    return sqrt(dot_product);
}

template <typename T>
vector<T> vector<T>::cross_product(const vector<T> &a, const vector<T> &b)
{
    T a2,b2;
    
    if (a.dimension==2)
        a2=0;
    else
        a2=a.component[2];
    
    if (b.dimension==2)
        b2=0;
    else
        b2=b.component[2];
        
    vector<T> c(3,0);

    c.set(0,a.component[1]*b2-a2*b.component[1]);
    c.set(1,a2*b.component[0]-a.component[0]*b2);
    c.set(2,a.component[0]*b.component[1]-a.component[1]*b.component[0]);

    return c;
}

template <typename T>
T vector<T>::findAngle(vector<T> vec){
    T theta, inProd, d;
    vector<T> cProd(3);

    inProd= (*this)*vec;
    d= inProd/(this->l2norm()*vec.l2norm());
    if (fabs(fabs(d)-1) < pow(10,-4) && d > 0 ) // acos can't compute numbers that are very close to +-1
        d=1.;                           //  a bit larger than +-1 and return nan
    else if (fabs(fabs(d)-1) < pow(10,-4) && d < 0 )
        d=-1.;
    theta= acos(d);
    if ( theta >= -M_PI && theta <= -(M_PI - M_PI/15) )
        theta= 2*M_PI+theta;
    
    return theta;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////operators//////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T &vector<T>::operator() (int i) const
//reads the i-th component of the vector
{
    return component[i];
}

template <typename T>
vector<T> vector<T>::operator+ (const vector<T> &right) const
{
    vector<T> result=*this;
    for (int i=0; i<dimension; i++)
        result.component[i]+=right.component[i];

    return result;
}

template <typename T>
vector<T> vector<T>::operator- (const vector<T> &right) const
{
    vector<T> result=*this;
    for (int i=0; i<dimension; i++)
        result.component[i]-=right.component[i];

    return result;
}

template <typename T>
vector<T> vector<T>::operator/ (const T &right) const
{
    vector<T> result= *this;
    for (int i=0; i<dimension; i++)
        result.component[i]/=right;

    return result;
}

//inner product
template <typename T>
vector<T> vector<T>::operator* (const T &right) const
{
    vector<T> result= *this;
    for (int i=0; i<dimension; i++)
        result.component[i]*=right;

    return result;
}

template<typename T>
vector<T> operator* (const T &left, vector<T> &right)
{
    vector<T> result=right;
    for (int i=0; i<result.dimension; i++)
        result.component[i]*=left;

    return result;
}
 
template <typename T>
T vector<T>::operator* (const vector<T> &right) const
{
    T dot_product=0;

    for (int i=0;i<dimension;i++)
        dot_product+=component[i]*right.component[i];

    return dot_product;
}

template <typename T>
inline vector<T> vector<T>::operator+= (const vector<T> &right){
    for (int i=0; i<dimension;i++)
        component[i]+= right.component[i];
    return *this;
}

template <typename T>
inline vector<T> vector<T>::operator-= (const vector<T> &right){
    for (int i=0; i<dimension;i++)
        component[i]-= right.component[i];
    return *this;
}

template <typename T>
inline vector<T> vector<T>::operator/= (const T &right){
    for (int i=0; i<dimension;i++)
        component[i]/= right;
    return *this;
}

template <typename T>
inline vector<T> vector<T>::operator*= (const T &right){
    for (int i=0; i<dimension;i++)
        component[i]*= right;
    return *this;
}

template <typename T>
vector<T> &vector<T>::operator= (const vector<T> &right)
{
    if (this==&right)
        return *this;

    if (component)
        delete [] component;
    component= new T[right.dimension];
    dimension=right.dimension;
    
    for (int i=0;i<dimension;i++)
        component[i]=right.component[i];

    return *this;
}

template<typename T>
int operator== (const vector<T> &v1, const vector<T> &v2){
    if (v1.dimension != v2.dimension)
    {
        std::cout<<"Error: Vectors must have the same dimension in order to check if they are equal.\n";
        exit(0);
    }
    int length=v1.dimension;
    int sOk=0;
    
    
    for (int i=0; i<length; i++)
        if (v1(i) == v2(i))
            sOk+=1;
    
    if (sOk == length)
        return 1;
    else
        return 0;
}

template <typename T>
std::ostream &operator<< (std::ostream &left, const vector<T> &right)
{
    left<<"(";
    for (int i=0;i<right.dimension-1;i++)
        left<<right.component[i]<<",";

    left<<right.component[right.dimension-1]<<")";

    return left;
}
#endif /* vector_h */   
