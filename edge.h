//
//  edges.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 14/12/21.
//

#ifndef edge_h
#define edge_h

#include "vector.h"
#include "linkedList.h"
#include "node.h"

template<typename T,int N>
//T: specifies the arithmetic type used, e.g. int, float, double...
//N: specifies the dimensionality of the nodes used in the cell. They can be 1 (1-D), 2 (2-D), 3 (3-D).
class edge{
public:
    edge(node<T,N> *n1=NULL, node<T,N> *n2=NULL);//default constructor
    edge(const edge<T,N> &e); //copy constructor
    ~edge();//destructor
    
    void setIndex(int i);
    void setEdgeLogFr(int i);
    int getIndex() const;
    int getEdgeLogFr() const;
    
    void setNode(int i,node<T,N> *nI);
    node<T,N> *getNode(int i) const;
    
    void setSharingCell(int i, void *cellAddress);
    void *getSharingCell(int i) const;
    const linkedList<void *> *getSharingCellNextPointer(int i) const;
    void dropSharingCell(void *cellAddress);
    void dropFirstSharingCell();
    int getSharingCellsNum() const;
    void moreSharingCells(void *cellAddress=NULL);//increases the value of sharingCells by 1.
    int lessSharingCells(void *cellAddress=NULL);/*decreases the value of sharingCells by 1.returns 1 if sharingCells==0 and 0
                                                  if sharingCells>0 */
    int noSharingCell() const;//returns 1 (true) if the node is an isolated one (does not belong to any cell).
    
    void setFace(int i, void *f1);
    void insertEdgeFace(void *f1);
    void dropFirstEdgeFace();
    void *getFace(int i) const;
    int getEdgeFacesLength() const;
 
    void computeDirectionVector();
    vector<T> getDirectionVector() const;
    T getLength() const;
    T findDirectionAngle();
    
    void setFlagTrue();
    void setFlagFalse();
    bool getFlag() const;
    
    template<typename t, int n> friend int operator== (const edge<t,n> &e1, const edge<t,n> &e2);
    template<typename t, int n> friend std::ostream &operator<< (std::ostream &left, const edge<t,n> &right);
private:
    int index;
    int edgeLogFr;
    int sharingCellsNum;/* indicates the number of cells that share this node.
                       0: the node belongs to no cell.
                      >0: multiple cells share this node, e.g. if sharingCells= 1 it means that
                          this node belongs to only one cell, if 2 it belongs to two cells, etc.*/
    node<T,N> *edgeNode[2];
    linkedList<void *> sharingCell; /*linked list containing the addresses of the cells that the edge is a part of.
                                            Pointers are of type void since cell is higher in the abstract object hierarchy*/
    linkedList<void *> sharingFace; /*linked list containing the addresses of the cells that the edge is a part of.
                                            Pointers are of type void since cell is higher in the abstract object hierarchy*/
    T length; //edge's length
    vector<T> directionVector; //vector that shows the direction of the edge. Its length is equal to edge's length.
    bool flagInList; //used to speed up insertion in lists. when a node is put in a list this flag is set to true.
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Definitions of functions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, int N>
inline edge<T,N>::edge(node<T,N> *n1, node<T,N> *n2): sharingCell(NULL), sharingFace(NULL), index(0), edgeLogFr(-1), sharingCellsNum(0), flagInList(false){
    edgeNode[0]= n1;
    edgeNode[1]= n2;
}

//copy constructor
template<typename T, int N>
inline edge<T,N>::edge(const edge<T,N> &e):sharingCell(e.sharingCell),sharingFace(e.sharingFace), index(e.index), sharingCellsNum(e.sharingCellsNum), edgeLogFr(e.edgeLogFr), flagInList(e.flagInList){
    edgeNode[0]= e.node[0];
    edgeNode[1]= e.node[1];
}

//destructor
template<typename T, int N>
inline edge<T,N>::~edge(){
}

template<typename T, int N>
inline void edge<T,N>::setIndex(int i){
    index= i;
}

template<typename T, int N>
inline void edge<T,N>::setEdgeLogFr(int i){
    edgeLogFr=i;
}

template<typename T, int N>
inline int edge<T,N>::getIndex() const{
    return index;
}

template<typename T, int N>
inline int edge<T,N>::getEdgeLogFr() const{
    return edgeLogFr;
}

template<typename T, int N>
inline void edge<T,N>::setNode(int i, node<T,N> *nI){
    edgeNode[i]= nI;
}

template<typename T, int N>
inline node<T,N> *edge<T,N>::getNode(int i) const{
    return edgeNode[i];
}

template<typename T, int N>
inline void edge<T,N>::insertEdgeFace(void *f1){
    sharingFace.insertNextItem(f1);
}

template<typename T, int N>
inline void edge<T,N>::setFace(int i,void *f1){
    sharingFace(i)= f1;
}

template<typename T, int N>
inline void edge<T,N>::dropFirstEdgeFace(){
    sharingFace.dropFirstItem();
}

template<typename T, int N>
inline int edge<T,N>::getEdgeFacesLength() const{
    return sharingFace.length();
}

template<typename T, int N>
inline void *edge<T,N>::getFace(int i) const{
    return sharingFace.getItem(i);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////SHARING CELLS FUNCTIONS//////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, int N>
inline void edge<T,N>::setSharingCell(int i,void *cellAddress){
    sharingCell(i)= cellAddress;
}

template<typename T, int N>
inline void edge<T,N>::dropFirstSharingCell(){
    sharingCell.dropFirstItem();
}

template<typename T, int N>
void edge<T,N>::dropSharingCell(void *cellAddress){
    int sharingCellLength= sharingCell.length();
    for (int i=0; i<sharingCellLength; i++)
    {
        if (sharingCell(i) == cellAddress)
        {
            sharingCell.dropThatItem(i);
            return;
        }
    }
    std::cout<<"The edge does not belong to this cell."<<std::endl;
}

template<typename T, int N>
inline const linkedList<void *> *edge<T,N>::getSharingCellNextPointer(int i) const{
    return sharingCell[i];
}

template<typename T, int N>
inline void *edge<T,N>::getSharingCell(int i) const{
    return sharingCell.getItem(i);
}

template<typename T, int N>
inline int edge<T,N>::getSharingCellsNum() const{
    return sharingCellsNum;
}

//increase the number of sharingCells by 1
template<typename T,int N>
inline void edge<T,N>::moreSharingCells(void *cellAddress){
    sharingCell.insertNextItem(cellAddress);
    sharingCellsNum++;
}

//decrease the number of sharingCells by 1
template<typename T,int N>
int edge<T,N>::lessSharingCells(void *cellAddress)
/*This function decreases the value of the variable sharingCells by 1 if it is larger than zero otherwise
 it leaves the variable sharingCells equal to zero. Then if the variable sharingCells is larger than zero
 it returns 0 otherwise it returns 1.*/
{
    if (cellAddress && sharingCellsNum > 1)
    {
        dropSharingCell(cellAddress);
    }
    return sharingCellsNum ? !(--sharingCellsNum) : 1;
}

//checks whether the edge belongs to no cell
/* This function returns 0 if the node belongs to at least 1 cell (sharingCells>0).
 It returns 1 if the node belongs to no cell (sharingCells==0).*/
template<typename T,int N>
inline int edge<T,N>::noSharingCell() const{
    return !sharingCellsNum;
}


template<typename T, int N>
std::ostream &operator<< (std::ostream &left, const edge<T,N> &right)
{
    for (int i=0; i<2; i++)
        left<<"node "<<i<<": "<<right.edgeNode[i]->getIndex()<<std::endl;

    return left;
}

template<typename T,int N>
void edge<T,N>::computeDirectionVector(){
    T lengthX, lengthY, lengthZ;
    
    lengthX= edgeNode[1]->operator[](0) - edgeNode[0]->operator[](0);
    lengthY= edgeNode[1]->operator[](1) - edgeNode[0]->operator[](1);
    
    if (N==3)
    {
        lengthZ= edgeNode[1]->operator[](2) - edgeNode[0]->operator[](2);
        directionVector.setDimension(3);
        directionVector.set(2, lengthZ);
    }
    else
    {
        lengthZ= 0;
        directionVector.setDimension(2);
    }
    
    directionVector.set(0, lengthX);
    directionVector.set(1, lengthY);
    length= sqrt(lengthX*lengthX + lengthY*lengthY + lengthZ*lengthZ);
}

template<typename T, int N>
inline T edge<T,N>::getLength() const{
    return length;
}

template<typename T, int N>
inline vector<T> edge<T,N>::getDirectionVector() const{
    return directionVector;
}

template<typename T, int N>
int operator== (const edge<T,N> &e1, const edge<T,N> &e2){
    node<T,N> *cNode;
    int sOk=0;
    
    for (int i=0; i<2; i++)
    {
        cNode= e1.getNode(i);
        for (int j=0; j<2; j++)
            if (*cNode == *(e2.getNode(j)))
            {
                sOk+=1;
                break;
            }
    }
    if (sOk==2) //if sOk==2 it means that there where 2 true ifs, so all the 2 nodes are the same for e1 and e2
        return 1;
    else
        return 0;
}

template<typename T, int N>
T edge<T,N>::findDirectionAngle(){
    double pi=atan(1.)*4.;
    T theta, x, y;
    
    x=directionVector(0);
    y=directionVector(1);

    if (fabs(x) < pow(10,-10) && y > 0)
        theta= pi/2.;
    else if (fabs(x) < pow(10,-10) && y < 0)
        theta= 3*pi/2.;
    else if (fabs(y) < pow(10,-10) && x > 0)
        theta= 0.;
    else if (fabs(y) < pow(10,-10) && x < 0)
        theta= pi;
    else
        theta= atan2(y,x);
    
    return theta;
}

template<typename T, int N>
inline void edge<T,N>::setFlagTrue(){
    flagInList=true;
}

template<typename T, int N>
inline void edge<T,N>::setFlagFalse(){
    flagInList=false;
}

template<typename T, int N>
inline bool edge<T,N>::getFlag() const{
    return flagInList;
}

#endif /* edge_h */
