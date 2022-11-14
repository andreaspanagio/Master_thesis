//
//  mesh.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 25/11/21.
//

#ifndef mesh_h
#define mesh_h

#include <string>
#include <cmath>

#include "node.h"
#include "vector.h"
#include "cell.h"
#include "linkedList.h"

template<typename T, int M, int N>
//T: specifies the arithmetic type used, e.g. int, float, double...
//M: specifies the number of nodes to be used to create a cell.
//N: specifies the dimensionality of the nodes used in the cell. They can be 1 (1-D), 2 (2-D), 3 (3-D).
class mesh: public linkedList<cell<T,M,N> *>{
public:
    mesh();
    mesh(T &e);
    ~mesh();
    
    int indexing();//indexing the nodes in the mesh
    int indexingCells();
    void gridQuality() const;
    void setNodesDisplacements(T &xDisplace, T &yDisplace);
private:
    
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Definitions of functions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//constructor
template<typename T, int M, int N>
inline mesh<T,M,N>::mesh(){
}

template<typename T, int M, int N>
inline mesh<T,M,N>::mesh(T &e){
    this->item=e;
}

template<typename T, int M, int N>
mesh<T,M,N>::~mesh(){
    delete this->item;
    //recursive call of the mesh destructor. the mesh destructor is call recursively at the command delete temp
    //until this->next==NULL, which is true only for the last element of the mesh.
    if (this->next)
    {
        mesh<T,M,N> *temp;
        linkedList<cell<T,M,N> *> *tcNext;
        
        tcNext= this->next;
        //necessary to typecast next. otherwise delete temp will call destructor of linkedList and not this one.
        temp= (mesh<T,M,N> *)tcNext;
        delete temp;
    }
    this->next=NULL;
}

template<typename T, int M, int N>
int mesh<T,M,N>::indexing(){
    for (mesh<T,M,N> *runner= this; runner; runner= (mesh<T,M,N> *)runner->next)
        runner->item.resetIndices();/*the loop ends when the runner pointer is set to NULL, because
                                     then the ending criterion (...; runner;...)
                                     is 0 (NULL==0). This happens when the runner points to the last
                                     element of the linked list.*/
    
    int count=0;
    for (mesh<T,M,N> *runner= this; runner; runner= (mesh<T,M,N> *)runner->next)
        runner->item.indexing(count);
    return count;
}

template<typename T, int M, int N>
int mesh<T,M,N>::indexingCells(){
    int count=0;
    for (mesh<T,M,N> *runner= this; runner; runner= (mesh<T,M,N> *)runner->next)
        runner->item.setIndex(count++);
    return count;
}

template<typename T, int M, int N>
void mesh<T,M,N>::gridQuality() const{
    const linkedList<cell<T,M,N> *> *cellNext= this;
    cell<T,M,N> *cCell;
    int meshLength= this->length();
    T *cellQual= new T[meshLength];
    T mean, sdev, min, max;
    int invCells=0;
    
    cellQual[0]= this->item->cellQuality();
    mean= cellQual[0];
    min= cellQual[0];
    max= cellQual[0];
    for (int i=0; i<meshLength-1; i++)
    {
        cellNext= cellNext->readNext();
        cCell= cellNext->getItem(0);
        cellQual[i+1]= cCell->cellQuality();
        mean+= cellQual[i+1];
        if (cCell->getVolume() < 0.)
            invCells+=1;
        if (cellQual[i+1] < min)
            min= cellQual[i+1];
        if (cellQual[i+1] > max)
            max= cellQual[i+1];
        
    }
    mean/= meshLength;
    
    sdev=0.;
    for (int i=0; i<meshLength; i++)
        sdev+= pow(cellQual[i]-mean,2);
    sdev= sqrt(sdev)/sqrt((double) meshLength);
    
    std::cout<<"\n......................................................................\n";
    std::cout<<"Grid Quality     "<<"mean: "<<mean<<"    "<<"sdev: "<<sdev<<"    "<<"    "<<"min: "<<min<<"    "<<
               "max: "<<max<<"    "<<"inverted cells: "<<invCells<<"\n";
    std::cout<<"......................................................................\n\n";
    delete [] cellQual;
}
#endif /* mesh_h */
