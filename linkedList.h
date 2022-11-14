//
//  linkedList.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 31/12/21.
//

#ifndef linkedList_h
#define linkedList_h

template <typename T>
class linkedList{
public:
    linkedList();//default constructor
    linkedList(const T &t,linkedList<T> *N=NULL); //constructor
    linkedList(const linkedList<T> &ob);//copy constructor
    ~linkedList();//destructor
    
    T &operator() (int i);//read/write the i-th item of the linked list
    const T &getItem(int i) const; //read the i-th item of the linked list
    const linkedList<T> *operator[] (int i) const;//returns the next pointer of the i-th element of the linked list
    const linkedList<T> *readNext() const;//returns the pointer next of the first item of the linked list
    const linkedList<T> &first() const; //return the first element of the linked list.
    linkedList<T> &last();//returns the last element of the linked list
    int length() const;//returns the total number of elements in the linked list
    void append(const T &t);//inserts a new item at the end of linked list
    void insertFirstItem(const T &t);//inserts an item at the beginning of the linked list
    void insertNextItem(const T &t);//inserts an item right after the first item
    void dropFirstItem();//deletes the first item of the list
    void dropNextItem();//deletes the item right after the beginning
    void dropThatItem(int i);//deletes the item stored at position i
    
    template<typename t> friend linkedList<t> operator* (const linkedList<t> &left,const linkedList<t> &right);
    const linkedList<T> &operator+= (linkedList<T> &right);
    const linkedList<T> &operator= (const linkedList<T> &right);
    template <typename U> friend std::ostream &operator<< (std::ostream &left, const linkedList<U> &right);
protected:
    T item;
    linkedList<T> *next;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Definitions of functions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//default constructor
template<typename T>
inline linkedList<T>::linkedList(){
    next= NULL;
}

//constructor
template<typename T>
inline linkedList<T>::linkedList(const T &t, linkedList<T> *N): item(t), next(N){
};

//copy constructor
//calls the copy constructor recursively until it bumps
//into the last element for which ob.next=NULL and sets the
//pointer next of the last copied element equal to NULL.
template<typename T>
inline linkedList<T>::linkedList(const linkedList<T> &ob): item(ob.item), next(ob.next ? new linkedList(*ob.next): NULL){
}

//destructor
//once is called for the first object then the destructor
//is called recursively for all the elements till the last
//element of the linked list.
template<typename T>
inline linkedList<T>::~linkedList(){
    if (next)
        delete next;
    next=NULL;
}

template<typename T>
inline const linkedList<T> *linkedList<T>::readNext() const{
    return next;
}

template<typename T>
const T &linkedList<T>::getItem(int i) const{
    if (i == 0)
        return item;
    
    linkedList<T> *that=next;
    for (int j=1; j<i; j++)
        that=that->next;
    return that->item;
}

template<typename T>
inline const linkedList<T> &linkedList<T>::first() const{
   return *this;
}

template<typename T>
linkedList<T> &linkedList<T>::last(){
    if (!next)
        return *this;
    
    linkedList<T> *that=next;
    while (that->next)
        that= that->next;
    
   return *that;
}

template<typename T>
int linkedList<T>::length() const{
    int length=1;
    if (!next)
        return length;
    
    linkedList<T> *that=next;
    while ((that= that->next))
        length+=1;
   return length+1;
}

template<typename T>
inline void linkedList<T>::append(const T &t){
    last().next= new linkedList<T>(t);
    
    return;
}

template<typename T>
inline void linkedList<T>::insertFirstItem(const T &t){
    next = new linkedList<T>(item,next);
    item = t;
}

template<typename T>
inline void linkedList<T>::insertNextItem(const T &t){
    next = new linkedList<T>(t,next);
}

template<typename T>
void linkedList<T>::dropNextItem(){
    if(next)
    {
        if(next->next)
        {
            linkedList<T> *keep = next;
            next = next->next;
            keep->next=NULL; //so that the destructor is not called recursively.
            delete keep;
        }
        else //deletes the last item of the linked list
        {
            delete next;
            next = NULL;
        }
    }
    else
        std::cout<<"Error: Cannot drop non existing next item."<<std::endl;
}

template<typename T>
void linkedList<T>::dropFirstItem()
{
    if(next)
    {
        item = next->item;
        dropNextItem();
    }
    else
        std::cout<<"Error: Cannot drop the only item.\n";
}

/*returns the item of the linked list stored at position i. i=0 is the first item.*/
template<typename T>
void linkedList<T>::dropThatItem(int i){
    if (i == 0)
        dropFirstItem();
    else if (i == 1)
        dropNextItem();
    else
    {
        linkedList<T> *that=next;
        for (int j=1; j<i-1; j++)
            that=that->next;
        that->dropNextItem();
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////operators//////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*returns the item of the linked list stored at position i. i=0 is the first item.*/
template<typename T>
T &linkedList<T>::operator()(int i){
    if (i == 0)
        return item;
    
    linkedList<T> *that=next;
    for (int j=1; j<i; j++)
        that=that->next;
    return that->item;
}

/*returns the next pointer of the linked list stored at position i. i=0 is the first item.*/
template<typename T>
const linkedList<T> *linkedList<T>::operator[](int i) const{
    if (i == 0)
        return next;
    
    linkedList<T> *that=next;
    for (int j=1; j<i; j++)
        that=that->next;
    return that->next;
}

template<typename T>
linkedList<T> operator+ (const linkedList<T> &left, const linkedList<T> &right){
    linkedList<T> result;
    const linkedList<T> *leftNext, *rightNext;
    rightNext= &right.first();
    while (rightNext)
    {
        result.insertNextItem(rightNext->getItem(0));
        rightNext= rightNext->readNext();
    }
    leftNext= &left.first();
    while (leftNext)
    {
        result.insertNextItem(leftNext->getItem(0));
        leftNext= leftNext->readNext();
    }
    result.dropFirstItem();
    return result;
}

template<typename T>
inline const linkedList<T> &linkedList<T>::operator+= (linkedList<T> &right){
    last()= right.first();
    //last().next=&right;
    return *this;
}

template<typename T>
const linkedList<T> &linkedList<T>::operator= (const linkedList<T> &right){
    if (this==&right)
        return *this;
    
    item = right.getItem(0);
    
    if (next)
    {
        if (right.next)
            *next = *right.next;
        else
        {
            delete next;
            next = NULL;
        }
    }
    else
        if (right.next)
            next= new linkedList<T>(*right.next);
    
    return *this;
}

template<typename T>
std::ostream &operator<< (std::ostream &left, const linkedList<T> &right){
    static int i=-1;
    i++;
    left<<"item "<<i<<": "<<std::endl<<right.item<<std::endl;
    
    if (right.readNext())
        std::cout<<*right.readNext();
    else
        i=-1; // resets static variable i
    
    return left;
}


#endif /* linkedList_h */
