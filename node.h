//
//  node.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 2/12/21.
//

#ifndef node_h
#define node_h

#include "vector.h"
#include "linkedList.h"

template<typename T,int N>
//T: specifies the arithmetic type used, e.g. int, float, double...
//N: specifies the dimensionality of the node N can be 1 (1-D), 2 (2-D) or 3 (3-D).
class node: public vector<T>{
public:
    node(const T &loc=0., int ind=-1, int sharing=0, int log=-1);//constructor
    node(const T &locX, const T &locY, const T &locZ=0., int ind=-1, int sharing=0, int log=-1);//2-D/3-D node constructor
    node(const node &n);//copy constructor
    ~node();//constructor
    
    void setIndex(int i);//sets the node's index equal to i
    int getIndex() const;//returns the node's index
    void setLogFr(int i);//sets the node's logFr equal to i
    int getLogFr() const;//returns the node's logFr
    int getSharingCells() const;//returns the number of sharingCells
    void moreSharingCells(void *cellAddress=NULL);//increases the value of sharingCells by 1.
    int lessSharingCells(void *cellAddress=NULL);/*decreases the value of sharingCells by 1.returns 1 if sharingCells==0 and 0
                                                  if sharingCells>0 */
    int noSharingCell() const;//returns 1 (true) if the node is an isolated one (does not belong to any cell).
    void *getSharingCellAddress(int i);
    linkedList<node<T,N> *> getNeighborNodes() const;
    void insertNeighborNode(node<T,N> * const neiNode);
    void dropNeighborNode(int i);
    const linkedList<void *> *getSharingCellNextPointer(int i) const;
    void dropSharingCellAddress(void *cellAddress);
    void dropFirstSharingCellAddress();
    
    void setFlagTrue();
    void setFlagFalse();
    bool getFlag();
    
    const node<T,N> &operator= (const node &right);
    vector<T> operator() () const;//read the location of the node
    void setLocationX(T x);//write the X-component of the location of the node
    void setLocationY(T y);//write the Y-component of the location of the node
    void setLocationZ(T z);//write the Z-component of the location of the node
    T operator[] (int i) const;/*read the i-th component of the node's location vector(i=0 -> x-component,
                          i=1 -> y-component, i=2 -> z-component)*/
    
    template<typename t, int n> friend int operator== (const node<t,n> &n1, const node<t,n> &n2);
    template <typename t, int n> friend std::ostream &operator<< (std::ostream &left, const node<t,n> &right);
private:
    int index;//a unique number used to identify the node when used in a mesh
    int sharingCells;/* indicates the number of cells that share this node.
                       0: the node belongs to no cell.
                      >0: multiple cells share this node, e.g. if sharingCells= 1 it means that
                          this node belongs to only one cell, if 2 it belongs to two cells, etc.*/
    linkedList<void *> sharingCellAddress; /*linked list containing the addresses of the cells that the node is a part of.
                                            Pointers are of type void since cell is higher in the abstract object hierarchy*/
    linkedList<node<T,N> *> neighborNodes; //linked list containing the addresses of the neighboring nodes of the node
    int logFr; /*a number betwen 0 and 4 that indicates if the node belongs to the boundary or not
                 0: node does not belong to the boundary (inner node)
                 1:
                 2:
                 3: node belongs to the inner boundary
                 4: node belongs to the outer boundary */
    bool flagInList; //used to speed up insertion in lists. when a node is put in a list this flag is set to true.
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Definitions of functions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//constructor
template<typename T,int N>
inline node<T,N>::node(const T &loc, int ind, int sharing, int log):
vector<T>(N,loc),index(ind), sharingCells(sharing), logFr(log), flagInList(false),
sharingCellAddress(NULL), neighborNodes(NULL){
}

//2-d/3-D node constructor
template<typename T,int N>
inline node<T,N>::node(const T &locX, const T &locY, const T &locZ, int ind, int sharing, int log):
vector<T>(locX,locY,locZ), index(ind), sharingCells(sharing), logFr(log), flagInList(false),
sharingCellAddress(NULL), neighborNodes(NULL){
}

//copy constructor
template<typename T,int N>
inline node<T,N>::node(const node &n): vector<T>(n()), index(n.index), sharingCells(n.sharingCells),
logFr(n.logFr),sharingCellAddress(n.sharingCellAddress), neighborNodes(n.neighborNodes),flagInList(n.flagInList){
}

//destructor
template<typename T,int N>
inline node<T,N>::~node(){
}

//set index equal to i
template<typename T,int N>
inline void node<T,N>::setIndex(int i){
    index=i;
}

//get index of the node
template<typename T,int N>
inline int node<T,N>::getIndex() const{
    return index;
}

//set logFr equal to i
template<typename T,int N>
inline void node<T,N>::setLogFr(int i){
    logFr=i;
}

//get logFr of the node
template<typename T,int N>
inline int node<T,N>::getLogFr() const{
    return logFr;
}

//get the number of sharingCells
template<typename T,int N>
inline int node<T,N>::getSharingCells() const{
    return sharingCells;
}

//increase the number of sharingCells by 1
template<typename T,int N>
inline void node<T,N>::moreSharingCells(void *cellAddress){
    sharingCellAddress.insertNextItem(cellAddress);
    sharingCells++;
}

//decrease the number of sharingCells by 1
template<typename T,int N>
int node<T,N>::lessSharingCells(void *cellAddress)
/*This function decreases the value of the variable sharingCells by 1 if it is larger than zero otherwise
 it leaves the variable sharingCells equal to zero. Then if the variable sharingCells is larger than zero
 it returns 0 otherwise it returns 1.*/
{
    if (cellAddress && sharingCells>1)
    {
        dropSharingCellAddress(cellAddress);
    }
    return sharingCells ? !(--sharingCells) : 1;
}

template<typename T, int N>
inline void node<T,N>::dropFirstSharingCellAddress(){
    sharingCellAddress.dropFirstItem();
}

template<typename T, int N>
void node<T,N>::dropSharingCellAddress(void *cellAddress){
    int sharingCellLength= sharingCellAddress.length();
    for (int i=0; i<sharingCellLength; i++)
    {
        if (sharingCellAddress(i) == cellAddress)
        {
            sharingCellAddress.dropThatItem(i);
            return;
        }
    }
    std::cout<<"The node does not belong to this cell."<<std::endl;
}

template<typename T, int N>
inline void *node<T,N>::getSharingCellAddress(int i){
    return sharingCellAddress(i);
}

template<typename T, int N>
inline const linkedList<void *> *node<T,N>::getSharingCellNextPointer(int i) const{
    return sharingCellAddress[i];
}

//checks whether the node belongs to no cell
/* This function returns 0 if the node belongs to at least 1 cell (sharingCells>0).
 It returns 1 if the node belongs to no cell (sharingCells==0).*/
template<typename T,int N>
inline int node<T,N>::noSharingCell() const{
    return !sharingCells;
}

// returns the location of the node
template<typename T,int N>
vector<T> node<T,N>::operator() () const{
    vector<T> location(this->dimension);
    for (int i=0; i<this->dimension; i++)
        location.set(i,this->component[i]);
    
    return location;
}

// returns the i-th component of the node's location vector
template<typename T,int N>
inline T node<T,N>::operator[] (int i) const{
    return this->component[i];
}

//write the X-component of the location of the node
template<typename T,int N>
inline void node<T,N>::setLocationX(T x){
    this->component[0]= x;
}

//write the Y-component of the location of the node
template<typename T,int N>
inline void node<T,N>::setLocationY(T y){
    this->component[1]= y;
}

//write the Z-component of the location of the node
template<typename T,int N>
inline void node<T,N>::setLocationZ(T z){
    this->component[2]= z;
}

template<typename T, int N>
inline linkedList<node<T,N> *> node<T,N>::getNeighborNodes() const{
    return neighborNodes;
}

template<typename T, int N>
inline void node<T,N>::insertNeighborNode(node<T,N> * const neiNode){
    neighborNodes.insertNextItem(neiNode);
}

template<typename T, int N>
inline void node<T,N>::dropNeighborNode(int i){
    neighborNodes.dropThatItem(i);
}

template<typename T, int N>
inline void node<T,N>::setFlagTrue(){
    flagInList=true;
}

template<typename T, int N>
inline void node<T,N>::setFlagFalse(){
    flagInList=false;
}

template<typename T, int N>
inline bool node<T,N>::getFlag(){
    return flagInList;
}

//assignment operator
template<typename T,int N>
const node<T,N> &node<T,N>::operator= (const node<T,N> &right){
    if (this == &right)
        return *this;
    
    delete [] this->component;
    this->dimension=right.dimension;
    this->component= new T[this->dimension];
    for (int i=0; i<this->dimension; i++)
        this->component[i]= right.component[i];
    index= right.index;
    sharingCells= right.sharingCells;
    logFr=right.logFr;
    
    sharingCellAddress=right.sharingCellAddress;
    neighborNodes=right.neighborNodes;
    
    return *this;
}

template<typename T, int N>
int operator== (const node<T,N> &n1, const node<T,N> &n2){
    if (n1() == n2())
        return 1;
    else
        return 0;
}

//print node
template<typename T,int N>
std::ostream &operator<< (std::ostream &left, const node<T,N> &right){
    left<<"node "<<right.index<<": "<<"Location "<<right()<<", "
    <<" logFr: "<<right.logFr<<", "
    <<" Sharing Cells: "<<right.sharingCells<<std::endl;
    return left;
}

#endif /* node_h */
