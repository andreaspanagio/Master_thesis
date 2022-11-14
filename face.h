//
//  face.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 17/8/22.
//

#ifndef face_h
#define face_h

#include <iostream>
#include<cmath>

#include "vector.h"
#include "node.h"
#include "linkedList.h"


template <typename T, int N>
//T: specifies the arithmetic type used, e.g. int, float, double...
//N: specifies the dimensionality of the nodes used in the cell. They can be 1 (1-D), 2 (2-D), 3 (3-D).
class face
{
public:
    face(); //default constructor
    face(node<T,N> *n1, node<T,N> *n2, node<T,N> *n3);
    face(const face<T,N> &f); //copy constructor
    ~face(); //destructor
    
    void insertFaceNode(node<T,N> *n);
    void dropFirstFaceNode();
    void setIndex(int i);
    void setFaceLogFr(int i);
    void setFaceNode(int i, node<T,N> *n);
    
    void findArea();
    void findNormalVector();
    
    int getIndex() const;
    int getFaceLogFr() const;
    int getFaceNodesLength() const;
    bool noBoundaryFace() const; //checks if the face belongs to the boundary. returns true (1) if it is not a boundary face.
    node<T,N> * getFaceNode(int i) const;
    vector<T> getNormalVector() const;
    T getFaceArea() const;
    
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
    
    void setFlagTrue();
    void setFlagFalse();
    bool getFlag() const;
    
    template<typename t, int n> friend int operator== (const face<t,n> &f1, const face<t,n> &f2);
    template<typename t, int n> friend std::ostream &operator<< (std::ostream &left, const face<t,n> &right);
private:
    linkedList<node<T,N> *> faceNode;
    linkedList<void *> sharingCell; /*linked list containing the addresses of the cells that the edge is a part of.
                                      Pointers are of type void since cell is higher in the abstract object hierarchy.
                                      sharingCells can only be 1 (if it is a boundary face) or 2 (if it is an inner face)*/
    int sharingCellsNum;/* indicates the number of cells that share this node.
                       0: the node belongs to no cell.
                      >0: multiple cells share this node, e.g. if sharingCells= 1 it means that
                          this node belongs to only one cell, if 2 it belongs to two cells, etc.*/
    T faceArea;
    int index;
    int faceLogFr;
    vector<T> normalVector;
    bool flagInList; //used to speed up insertion in lists. when a node is put in a list this flag is set to true.
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Definitions of functions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//default constructor
template<typename T, int N>
inline face<T,N>::face(): faceNode(NULL), sharingCell(NULL), faceArea(0.), index(0), normalVector(), faceLogFr(-1), sharingCellsNum(0), flagInList(false){
}

template<>
inline face<double,3>::face(node<double,3> *n1, node<double,3> *n2, node<double,3> *n3): sharingCell(NULL), sharingCellsNum(0), faceArea(0.), index(0), normalVector(), faceLogFr(-1), flagInList(false){
    faceNode.insertNextItem(n1);
    faceNode.insertNextItem(n2);
    faceNode.insertNextItem(n3);
    faceNode.dropFirstItem();
}

//copy constructor
template<typename T, int N>
inline face<T,N>::face(const face<T,N> &f):faceNode(f.faceNode), sharingCell(f.sharingCell), faceArea(f.faceArea), normalVector(f.normalVector), faceLogFr(f.faceLogFr), sharingCellsNum(f.sharingCellsNum), flagInList(f.flagInList){
}

//destructor
template<typename T, int N>
face<T,N>::~face(){
}

template<typename T, int N>
inline void face<T,N>::insertFaceNode(node<T,N> *n){
    faceNode.insertNextItem(n);
}

template<typename T, int N>
inline void face<T,N>::dropFirstFaceNode(){
    faceNode.dropFirstItem();
}

template<typename T, int N>
inline void face<T,N>::setIndex(int i){
    index=i;
}

template<typename T, int N>
inline void face<T,N>::setFaceLogFr(int i){
    faceLogFr=i;
}

template<typename T, int N>
inline void face<T,N>::setFaceNode(int i, node<T,N> *n){
    faceNode(i)=n;
}

//only gives the correct area for a triangular face.
template<typename T, int N>
void face<T,N>::findArea(){
    vector<T> n1(3),a(3),b(3);
    
    a= faceNode(0)->operator()()- faceNode(1)->operator()();
    b= faceNode(0)->operator()()- faceNode(2)->operator()();
    faceArea= 0.5*fabs((n1.cross_product(a,b)).l2norm());
}

template<typename T, int N>
void face<T,N>::findNormalVector(){
    vector<T> a(3),b(3);
    
    a= faceNode(0)->operator()()- faceNode(1)->operator()();
    b= faceNode(0)->operator()()- faceNode(2)->operator()();
    
    normalVector= normalVector.cross_product(a, b);
    normalVector= normalVector/normalVector.l2norm();
    
    if (faceNode(0)->operator()()*normalVector < 0)
        normalVector= normalVector*(-1);
}

template<typename T, int N>
inline int face<T,N>::getIndex() const{
    return index;
}

template<typename T, int N>
inline int face<T,N>::getFaceLogFr() const{
    return faceLogFr;
}

template<typename T, int N>
inline node<T,N> *face<T,N>::getFaceNode(int i) const{
    return faceNode.getItem(i);
}

template<typename T, int N>
inline int face<T,N>::getFaceNodesLength() const{
    return faceNode.length();
}

template<typename T, int N>
inline vector<T> face<T,N>::getNormalVector() const{
    return normalVector;
}

template<typename T, int N>
inline T face<T,N>::getFaceArea() const{
    return faceArea;
}

//checks if the face belongs to the boundary. if it does the face can only belong to one cell
template<typename T, int N>
inline bool face<T,N>::noBoundaryFace() const{
    return sharingCell.readNext();
    //if readNext equals NULL then it is translated to 0 boolean and the face is a boundary one.
    //if readNext returns a memory adress then it is a number and it is translated to boolean 1 and the face
    //is not a boundary one.
    //code is written like this to save time. could as well be written with an if statement but then it would be too slow.
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////SHARING CELLS FUNCTIONS//////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, int N>
inline void face<T,N>::setSharingCell(int i,void *cellAddress){
    sharingCell(i)= cellAddress;
}

template<typename T, int N>
inline void face<T,N>::dropFirstSharingCell(){
    sharingCell.dropFirstItem();
}

template<typename T, int N>
void face<T,N>::dropSharingCell(void *cellAddress){
    int sharingCellLength= sharingCell.length();
    for (int i=0; i<sharingCellLength; i++)
    {
        if (sharingCell(i) == cellAddress)
        {
            sharingCell.dropThatItem(i);
            return;
        }
    }
    std::cout<<"The face does not belong to this cell."<<std::endl;
}

template<typename T, int N>
inline const linkedList<void *> *face<T,N>::getSharingCellNextPointer(int i) const{
    return sharingCell[i];
}

template<typename T, int N>
inline void *face<T,N>::getSharingCell(int i) const{
    return sharingCell.getItem(i);
}

template<typename T, int N>
inline int face<T,N>::getSharingCellsNum() const{
    return sharingCellsNum;
}

//increase the number of sharingCells by 1
template<typename T,int N>
inline void face<T,N>::moreSharingCells(void *cellAddress){
    sharingCell.insertNextItem(cellAddress);
    sharingCellsNum++;
}

//decrease the number of sharingCells by 1
template<typename T,int N>
int face<T,N>::lessSharingCells(void *cellAddress)
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
inline int face<T,N>::noSharingCell() const{
    return !sharingCellsNum;
}


template<typename T, int N>
int operator== (const face<T,N> &f1, const face<T,N> &f2){
    node<T,N> *cNode;
    int sOk=0;
    
    for (int i=0; i<3; i++)
    {
        cNode= f1.getFaceNode(i);
        for (int j=0; j<3; j++)
            if (*cNode == *(f2.getFaceNode(j)))
            {
                sOk+=1;
                break;
            }
    }
    if (sOk==3) //if sOk==3 it means that there where 3 true ifs, so all the 3 nodes are the same for f1 and f2
        return 1;
    else
        return 0;
}

template<typename T, int N>
std::ostream &operator<< (std::ostream &left, const face<T,N> &right){
    left<<"face "<<right.index;
    left<<" is made of the following nodes:\n";
    for (int i=0; i<right.faceNode.length(); i++)
        left<<right.faceNode.getItem(i)->getIndex()<<"\n";
    left<<"It has a surface area of "<<right.faceArea<<"\n";
    left<<"and its normal vector is the following:\n"<<right.normalVector;
    
    return left;
}

template<typename T, int N>
inline void face<T,N>::setFlagTrue(){
    flagInList=true;
}

template<typename T, int N>
inline void face<T,N>::setFlagFalse(){
    flagInList=false;
}

template<typename T, int N>
inline bool face<T,N>::getFlag() const{
    return flagInList;
}

#endif /* face_h */
