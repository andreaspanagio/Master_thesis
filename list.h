//
//  list.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 23/11/21.
//

#ifndef list_h
#define list_h

template<typename T>
class list{
public:
    list(int n=0);//default constructor
    list(int n,T t);//constructor
    list(const list<T> &ob);//copy constructor
    ~list();//destructor
    
    void setListSize(int n);//set list size equal to n.
    int size() const;
    void put(int i, T t);// puts an object t at the i-th position of the list
    
    T &operator() (int i);//returns the item stored at the i-th position of the list
    T *operator[] (int i);//returns the i-th pointer of the list
    const list<T> &operator= (const list<T> &right);
    
    template <typename U> friend std::ostream &operator<< (std::ostream &left, const list<U> &right);
protected:
    int number;
    T **item;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Definitions of functions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//default constructor
template<typename T>
list<T>::list(int n){
    number= n;
    
    if (number==0)
        item= NULL;
    else
    {
        item= new T *[number];
        for (int i=0; i<number; i++)
            item[i]= NULL;
    }
}

//constructor
template<typename T>
list<T>::list(int n,T t){
    number= n;
    item= new T *[number];
    for (int i=0; i<number; i++)
        item[i]= new T(t);
}

//copy constructor
template<typename T>
list<T>::list(const list<T> &ob){
    number= ob.number;
    item= new T *[number];
    if (!item)
    
    for (int i=0;i<number;i++)
    {
        if (ob.item[i])
        {
            item[i]= new T(*ob.item[i]);
        }
        else
            item[i]=NULL;
    }
}

//destructor
template<typename T>
list<T>::~list(){
    for (int i=0; i<number;i++)
        delete item[i];
    delete [] item;
}

template<typename T>
void list<T>::setListSize(int n){
    number= n;
    
    if (number==0)
        {
            std::cout<<"Cannot set the size of list to zero.\n";
            exit(0);
        }
    else
    {
        item= new T *[number];
        for (int i=0; i<number; i++)
            item[i]= NULL;
    }
}

template<typename T>
inline int list<T>::size() const{
    return number;
}

template<typename T>
void list<T>::put(int i, T t){
    if (item[i])
        delete item[i];
    
    item[i]= new T(t);
}

//returns the item stored at the i-th position of the list
template<typename T>
T &list<T>::operator() (int i){
    if(item[i])
        return *item[i];
    else
    {
        std::cout<<"Error: This item is empty."<<std::endl;
        exit(0);
    }
}

template<typename T>
T *list<T>::operator[] (int i){
        return item[i];
}

template<typename T>
const list<T> &list<T>::operator= (const list<T> &right){
    if (this == &right)
        return *this;
    
    for (int i=0; i<number; i++)
        delete item[i];
    delete [] item;
    
    number= right.number;
    item= new T * [number];
    
    for (int i=0; i<number; i++)
        item[i]= new T(*right.item[i]);
    
    return *this;
}

template<typename T>
std::ostream &operator<< (std::ostream &left, const list<T> &right){
    for (int i=0;i<right.number;i++)
        left<<"item "<<i<<": "<<*right.item[i]<<std::endl;
    return left;
}

#endif /* list_h */
