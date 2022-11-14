//
//  cell.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 3/12/21.
//

#ifndef cell_h
#define cell_h

#include <cmath>

#include "vector.h"
#include "node.h"
#include "edge.h"
#include "face.h"

template<typename T, int M, int N>
//T: specifies the arithmetic type used, e.g. int, float, double...
//M: specifies the number of nodes to be used to create a cell (e.g. for a triangle use M=3).
//N: specifies the dimensionality of the nodes used in the cell. They can be 1 (1-D), 2 (2-D), 3 (3-D).
class cell{
public:
    
    cell();//constructor
    cell(node<T,N> &A, node<T,N> &B, node<T,N> &C);//triangle constructor
    cell(node<T,N> &A, node<T,N> &B, node<T,N> &C, node<T,N> &D);//tetrahedron constructor
    cell(const cell<T,M,N> &ob);//copy constructor
    ~cell();//destructor
    
    void resetIndices();//reset indices to -1
    void indexing(int &count);//indexing the unindexed vertices
    void setIndex(int i);//set the index of the cell to i
    int getIndex() const;//reads the index of the cell
    void putNode(int i, node<T,N> *n,bool flag1, bool flag2);//puts an existing dynamically allocated node as the i-th vertex of the cell
    void cellArea();//calculates the area of the cell
    void cellVolume();//calculates the volume of the cell (only for 3D cells)
    T getArea() const;//reads the area of the cell
    T getVolume() const;//reads the volume of the cell
    void setCellEdge(int i,edge<T,N> *edgeAdress=NULL);//set the address of the i-th edge
    void insertCellEdge(edge<T,N> *edgeAdress);//set the address of the i-th edge
    void dropFirstCellEdge();
    T cellEdgesLength() const; //return the number of edges in the cell.
    edge<T,N> *getCellEdge(int i) const;//reads the address of the edge
    void setCellFace(int i,face<T,N> *faceAdress=NULL);//set the address of the i-th face
    void insertCellFace(face<T,N> *faceAdress);//set the address of the i-th face
    void dropFirstCellFace();
    T cellFacesLength() const; //return the number of faces in the cell.
    face<T,N> *getCellFace(int i) const;//reads the address of the face
    
    void resetVertexes();
    T cellQuality() const;

    const cell<T,M,N> &operator=(cell<T,M,N> &right);
    node<T,N> *operator() (int i);//read/write i-th vertex adress
    const node<T,N> &operator[] (int i) const;//only read i-th vertex
    template<typename t, int m, int n> friend int operator< (const node<t,n> &nod, const cell<t,m,n> &cel);/*checks
    whether the node n is indeed a vertex in the cell. The symbol < here means belongs. It has been chosen
    because it looks a bit similar to the symbol 'belongs to' of the set theory*/
    template<typename t, int m, int n> friend int operator& (const cell<t,m,n> &e, const cell<t,m,n> &f);/*checks whether
    two cells share an edge*/
    template<typename t, int m, int n>friend std::ostream &operator<< (std::ostream &left, const cell<t,m,n> &right);
private:
    node<T,N> *vertex[M];
    linkedList<edge<T,N> *> cellEdge;
    linkedList<face<T,N> *> cellFace;
    int index;
    T area;
    T volume;
};
using triangle= cell<double,3,2>;
using tetrahedron= cell<double,4,3>;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Definitions of functions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//constructor
template<typename T, int M, int N>
cell<T,M,N>::cell(): index(-1), area(0.), volume(0.), cellEdge(NULL), cellFace(NULL){
    for (int i=0; i<M; i++)
    {
        vertex[i]= new node<T,N>;
        vertex[i]->moreSharingCells((void *)this);
    }
}

//triangle constructor
template<>
inline triangle::cell(node<double,2> &A, node<double,2> &B, node<double,2> &C): index(-1), area(0.), volume(0.), cellEdge(NULL), cellFace(NULL)
/*Works only if the nodes A,B,C that are passed to the constructor are stored dynamically. If so, then their addresses are copied to
 the vertexes of the cell and then the sharingCells variable of the node is incremented by one meaning that the node now is a part
 of one more cell.*/
{
    vertex[0]= &A;
    vertex[1]= &B;
    vertex[2]= &C;
    
    for (int i=0; i<3; i++)
        vertex[i]->moreSharingCells((void *)this);
}

//tetrahedron constructor
template<>
inline tetrahedron::cell(node<double,3> &A, node<double,3> &B, node<double,3> &C, node<double,3> &D): index(-1), area(0.), volume(0.), cellEdge(NULL), cellFace(NULL)
/*Works only if the nodes A,B,C,D that are passed to the constructor are stored dynamically. If so, then their addresses are copied to the vertexes of the cell and then the sharingCells variable of the node is incremented by one meaning that the node now is a part of one more cell.*/
{
    vertex[0]= &A;
    vertex[1]= &B;
    vertex[2]= &C;
    vertex[3]= &D;
    
    for (int i=0; i<4; i++)
        vertex[i]->moreSharingCells((void *)this);
}

//copy constructor
template<typename T, int M, int N>
inline cell<T,M,N>::cell(const cell<T,M,N> &ob):index(ob.index), area(ob.area), volume(ob.volume){
    int cellEdgesNum= ob.cellEdge.length();
    int cellFacesNum= ob.cellEdge.length();
    for (int i=0; i<M; i++)
    {
        vertex[i]= ob.vertex[i];
        vertex[i]->moreSharingCells((void *)this);
    }
    
    cellEdge= ob.cellEdge;
    for (int i=0; i<cellEdgesNum; i++)
        cellEdge(i)->moreSharingCells((void *)this);
    
    cellFace= ob.cellFace;
    for (int i=0; i<cellFacesNum; i++)
        cellFace(i)->moreSharingCells((void *)this);
}

//destructor
template<typename T, int M, int N>
cell<T,M,N>::~cell(){
    int cellEdgesNum= cellEdge.length();
    int cellFacesNum= cellFace.length();
    
    for (int i=0; i<M; i++)
        if (vertex[i]->lessSharingCells((void *)this))
            delete vertex[i];
    for (int i=0; i<cellEdgesNum; i++)
        if (cellEdge(i)->lessSharingCells((void *)this))
            delete cellEdge(i);
    for (int i=0; i<cellFacesNum; i++)
        if (cellFace(i)->lessSharingCells((void *)this))
            delete cellFace(i);
}

//resets the Indices of all the vertexes of the cell and gives them the value -1
template<typename T, int M, int N>
inline void cell<T,M,N>::resetIndices(){
    for (int i=0; i<M;i++)
        vertex[i]->setIndex(-1);
}

/*indexes all the vertexes of the cell starting with the one stored at vertex[0]
 to which it gives the value count*/
template<typename T, int M, int N>
void cell<T,M,N>::indexing(int &count){
    for (int i=0; i<M; i++)
    {
        if (vertex[i]->getIndex()<0)
            vertex[i]->setIndex(count++);
    }
}

//set the index of the cell equal to i
template<typename T, int M, int N>
inline void cell<T,M,N>::setIndex(int i){
    index=i;
}

//returns the index of the cell
template<typename T, int M, int N>
inline int cell<T,M,N>::getIndex() const{
    return index;
}

//puts an existing node as the i-th vertex of the cell
template<typename T, int M, int N>
void cell<T,M,N>::putNode(int i, node<T,N> *n, bool flag1, bool flag2){
    if (vertex[i] == n)
        return;//user tries to set the vertex equal to the node it already points to.
    
    if (flag1)
        if(vertex[i]->lessSharingCells((void *)this) )
            delete vertex[i];
    
    vertex[i]=n;
    if (flag2)
        vertex[i]->moreSharingCells((void *)this);
}

template<typename T, int M, int N>
inline void cell<T,M,N>::setCellEdge(int i, edge<T,N> *edgeAdress){//set the address of the i-th edge
    if (cellEdge(i) == edgeAdress)
        return;//user tries to set the vertex equal to the node it already points to.
    
    if (cellEdge(i))
        if(cellEdge(i)->lessSharingCells((void *)this) )
            delete cellEdge(i);
    
    cellEdge(i)= edgeAdress;
    cellEdge(i)->moreSharingCells((void *)this);
}

template<typename T, int M, int N>
inline void cell<T,M,N>::insertCellEdge(edge<T,N> *edgeAdress){//set the address of the i-th edge
    cellEdge.insertNextItem(edgeAdress);
    cellEdge(1)->moreSharingCells((void *)this);
}

template<typename T, int M, int N>
inline void cell<T,M,N>::dropFirstCellEdge(){
    if (!cellEdge(0))
    {
        cellEdge.dropFirstItem();
        return;
    }
    else
    {
        if (cellEdge(0)->lessSharingCells((void *)this))
            delete cellEdge(0);
        cellEdge.dropFirstItem();
    }
}

template<typename T, int M, int N>
inline edge<T,N> *cell<T,M,N>::getCellEdge(int i) const{//reads the address of the edge
    return cellEdge.getItem(i);
}

template<typename T, int M, int N>
inline T cell<T,M,N>::cellEdgesLength() const{//reads the address of the edge
    return cellEdge.length();
}

template<typename T, int M, int N>
inline void cell<T,M,N>::setCellFace(int i, face<T,N> *faceAdress){//set the address of the i-th edge
    if (cellFace(i) == faceAdress)
        return;//user tries to set the vertex equal to the node it already points to.
    
    if (cellFace(i))
        if(cellFace(i)->lessSharingCells((void *)this) )
            delete cellFace(i);
    
    cellFace(i)= faceAdress;
    cellFace(i)->moreSharingCells((void *)this);
}

template<typename T, int M, int N>
inline void cell<T,M,N>::insertCellFace(face<T,N> *faceAdress){//set the address of the i-th edge
    cellFace.insertNextItem(faceAdress);
    cellFace(1)->moreSharingCells((void *)this);
}

template<typename T, int M, int N>
inline void cell<T,M,N>::dropFirstCellFace(){
    if (!cellFace(0))
    {
        cellFace.dropFirstItem();
        return;
    }
    else
    {
        if (cellFace(0)->lessSharingCells((void *)this))
            delete cellFace(0);
        cellFace.dropFirstItem();
    }
}

template<typename T, int M, int N>
inline face<T,N> *cell<T,M,N>::getCellFace(int i) const{//reads the address of the edge
    return cellFace.getItem(i);
}

template<typename T, int M, int N>
inline T cell<T,M,N>::cellFacesLength() const{//reads the address of the edge
    return cellFace.length();
}

//computes the area of a triangle cell
template<>
void triangle::cellArea(){
    vector<double> posVector0, posVector1, posVector2, crossProduct(3);
    
    posVector0= vertex[0]->operator()();
    posVector1= vertex[1]->operator()();
    posVector2= vertex[2]->operator()();
    crossProduct=crossProduct.cross_product(posVector1-posVector0, posVector2-posVector0);
    
    area= 0.5*crossProduct(2);
}

//computes the surface area of a tetrahedron cell
template<>
void tetrahedron::cellArea(){
    vector<double> posVector0, posVector1, posVector2, posVector3;
    vector<double> crossProduct0(3), crossProduct1(3), crossProduct2(3), crossProduct3(3);
    
    posVector0= vertex[0]->operator()();
    posVector1= vertex[1]->operator()();
    posVector2= vertex[2]->operator()();
    posVector3= vertex[3]->operator()();
    crossProduct0=crossProduct0.cross_product(posVector0-posVector2, posVector0-posVector1);
    crossProduct1=crossProduct1.cross_product(posVector0-posVector1, posVector0-posVector3);
    crossProduct2=crossProduct2.cross_product(posVector1-posVector2, posVector2-posVector3);
    crossProduct3=crossProduct3.cross_product(posVector1-posVector2, posVector0-posVector2);
    
    area= 0.5* ( fabs(crossProduct0.l2norm()) + fabs(crossProduct1.l2norm()) +
                fabs(crossProduct2.l2norm()) + fabs(crossProduct3.l2norm()) );
}

//computes the volume of a tetrahedron cell
template<>
void tetrahedron::cellVolume(){
    vector<double> posVector0, posVector1, posVector2, posVector3;
    vector<double> crossProduct0(3);
    
    posVector0=vertex[0]->operator()();
    posVector1=vertex[1]->operator()();
    posVector2=vertex[2]->operator()();
    posVector3=vertex[3]->operator()();
    crossProduct0=crossProduct0.cross_product(posVector1-posVector3, posVector2-posVector3);
    volume= ((posVector0-posVector3)*crossProduct0)/6.;
}

//returns the area of the cell
template<typename T, int M, int N>
inline T cell<T,M,N>::getArea() const{
    return area;
}

//returns the volume of the cell
template<typename T, int M, int N>
inline T cell<T,M,N>::getVolume() const{
    return volume;
}

//assignment operator
template<typename T, int M, int N>
const cell<T,M,N> &cell<T,M,N>::operator=(cell<T,M,N> &right){
    int cellEdgesNum= cellEdge.length();
    int cellFacesNum= cellFace.length();
    if (this == &right)
        return *this;
    
    index=right.index;
    area=right.area;
    volume= right.volume;
    for (int i=0; i<M; i++)
    {
        if (vertex[i]->lessSharingCells((void *)this))
            delete vertex[i];
    }
    
    for (int i=0; i<M; i++)
    {
        vertex[i]= right.vertex[i];
        vertex[i]->moreSharingCells((void *)this);
    }
    
    for (int i=0; i<cellEdgesNum; i++)
    {
        if (cellEdge(i)->lessSharingCells((void *)this))
            delete cellEdge(i);
    }
    cellEdge= right.cellEdge;
    cellEdgesNum= cellEdge.length();
    for (int i=0; i<cellEdgesNum; i++)
        cellEdge(i)->moreSharingCells((void *)this);
    
    for (int i=0; i<cellFacesNum; i++)
    {
        if (cellFace(i)->lessSharingCells((void *)this))
            delete cellFace(i);
    }
    cellFace= right.cellFace;
    cellFacesNum= cellFace.length();
    for (int i=0; i<cellFacesNum; i++)
        cellFace(i)->moreSharingCells((void *)this);
    
    return *this;
}

template<typename T, int M, int N>
inline node<T,N> *cell<T,M,N>::operator() (int i){
    return vertex[i];
}

template<typename T, int M, int N>
inline const node<T,N> &cell<T,M,N>::operator[] (int i) const{
    return *(vertex[i]);
}

template<typename T, int M, int N>
inline void cell<T,M,N>::resetVertexes(){
    for (int i=0; i<M; i++)
        vertex[i]=NULL;
}

template<>
double triangle::cellQuality() const{
    double metricRinRout, s=0., area, Rin, Rout;
    double a[3];
    int j;
    
    for (int i=0; i<3; i++)
    {
        if (i==2)
            j=0;
        else
            j= i+1;
        
        a[i]= sqrt(pow(vertex[i]->operator()()(0)-vertex[j]->operator()()(0),2)+
                pow(vertex[i]->operator()()(1)-vertex[j]->operator()()(1),2));
        s+=a[i];
    }
    s*= 0.5; //perimeter/2
    area= sqrt(s*(s-a[0])*(s-a[1])*(s-a[2])); //Heron formula
    Rin= area/s; //Radius in
    Rout= 0.25*(a[0]*a[1]*a[2]/area); //Radius out
    metricRinRout= Rout/Rin; //Rout/Rin metric (ideal= 2)
    
    return metricRinRout;
}

template<>
double tetrahedron::cellQuality() const{
    double L,d,V_ideal,V_cell;//,L;
    vector<double> posVec1(3), posVec2(3);
    double metricSkewness;
    
    
    L=0.;
    for (int i=0; i<3; i++)
    {
        posVec1= vertex[i]->operator()();
        for (int j=i+1; j<4; j++)
        {
            posVec2= vertex[j]->operator()();
            d=(posVec1 - posVec2).l2norm();
            if (d > L)
                L= d;
        }
    }
    V_ideal= (sqrt(2.)/12.)*pow(L,3);
    V_cell= fabs(this->volume);
    
    return metricSkewness= (V_ideal-V_cell)/V_ideal;
}

//operator belongs to
template<typename T, int M, int N>
int operator< (const node<T,N> &nod, const cell<T,M,N> &cel)
/* This operator checks if node belongs to a cell. It runs all the vertexes of the cell and checks whether
 their memory address is the same as the one of the node. If the node belongs to the cell then it returns the
 cell's vertex index plus one (so that zero can be used to define the case where the node does not belong to
 cell.*/
{
    for (int i=0;i<M;i++)
    {
        if (&nod == &cel[i])
            return i+1;
    }
    return 0;
}

//operator common edge
template<typename T, int M, int N>
int operator& (const cell<T,M,N> &e, const cell<T,M,N> &f)
/* The operator & checks whether two cells, let them be e and f, have a common edge. The only way to share a common edge
 is to share two nodes. This operator uses the previously declared operator <, which checks if a node belongs to a cell.
 If a node is found to belong to the cell then in the nested loop all the other vertexes, that their indexes are smaller
 than the one found, are checked (this is done to avoid the false case where it finds the same vertex twice and falsly
 returns that the cells indeed share an edge). This way all the vertexes are checked with one another. If the function
 indeed finds two common nodes then it returns 1 (true). Otherwise it returns 0 (false).*/
{
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<i; j++)
        {
            if((e[i]<f) && (e[j]<f))
                return 1;
        }
    }
    return 0;
}

template<typename T, int M, int N>
std::ostream &operator<< (std::ostream &left, const cell<T,M,N> &right){
    left<<"CellIndex: "<<right.index<<"\n";
    left<<"CellArea: "<<right.area<<"\n";
    for (int i=0; i<M; i++)
        left<<"Vertex "<<i<<": "<<*right.vertex[i]<<"\n";
    for(int i=0; i<2*M; i++)
    {
        if (right.cellEdge[i] == NULL)
            continue;
        left<<"edge "<<i<<":"<<" ";
        left<<"node 1: "<<(right.cellEdge[i]->getNode1())->getIndex()<<"\n";
        left<<"        "<<"node 2: "<<(right.cellEdge[i]->getNode2())->getIndex()<<"\n";
        left<<"\n";
    }
    return left;
}
#endif /* cell_h */
