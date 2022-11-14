//
//  vtkUnstrMeshTranslator.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 14/7/22.
//
//vtk is the format that ParaView reads for unstructured meshes
#ifndef vtkUnstrMeshTranslator_h
#define vtkUnstrMeshTranslator_h

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

template<typename T, int M, int N>
//T: specifies the arithmetic type used, e.g. int, float, double...
//M: T specifies the number of nodes in the cell (e.g. for a triangle use M=3).
//N: specifies the dimensionality of the nodes used in the cell. They can be 1 (1-D), 2 (2-D), 3 (3-D).
class vtkUnstrMeshTranslator{
public:
    vtkUnstrMeshTranslator(const createMesh<T,M,N> &m, const std::string filename, const std::string pathname=""); //constructor
private:
    
};

template<typename T, int M, int N>
vtkUnstrMeshTranslator<T,M,N>::vtkUnstrMeshTranslator(const createMesh<T,M,N> &m, const std::string filename, const std::string pathname){
    std::ifstream inData;
    std::fstream outData;
    std::string filename1= pathname + filename + ".vtk";
    node<T,N> *cNode;
    tetrahedron *cTetra;
    std::string keyword;
    int nodesNumber= m.getNodesNumber();
    int cellsNumber= m.getCellsNumber();
    const linkedList<cell<T,M,N> *> *cellNext;
    
    //write the .vtk file
    //open filename1 to write the necessary data
    outData.open(filename1,std::ios::out);
    if (!outData)
    {
        std::cout<<"\nError: Could not create file"<<filename1<<".nod.\n";
        exit(0);
    }
    outData.precision(17);
    
    outData<<"# vtk DataFile Version 2.0\n";
    outData<<"squares unstructured grid\n";
    outData<<"ASCII\n";
    outData<<"DATASET UNSTRUCTURED_GRID\n\n";
    outData<<"POINTS "<<nodesNumber<<" float\n";
    
    cellNext= &m().first();
    for (int i=0; i<nodesNumber; i++)
    {
        cTetra= cellNext->getItem(0);
        for (int j=0; j<M; j++)
        {
            cNode= cTetra->operator()(j);
            if (cNode->getFlag())
                continue;
            outData<<cNode->operator[](0)<<" "<<cNode->operator[](1)<<" "<<cNode->operator[](2)<<"\n";
            cNode->setFlagTrue();
        }
        cellNext= cellNext->readNext();
    }
    
    cellNext= &m().first();
    for (int i=0; i<cellsNumber; i++)
    {
        cTetra= cellNext->getItem(0);
        cellNext= cellNext->readNext();
        
\
        for (int j=0; j<M; j++)
        {
            cNode= cTetra->operator()(j);
            if (!cNode->getFlag())
                continue;
            cNode->setFlagFalse();
        }
    }

    outData<<"\n\nCELLS ";
    outData<<cellsNumber<<" ";
    outData<<5*cellsNumber<<"\n";
    cellNext= &m().first();
    for (int i=0; i<cellsNumber; i++)
    {
        cTetra= cellNext->getItem(0);
        outData<<4<<" ";
        outData<<cTetra->operator[](0).getIndex()-1<<" ";
        outData<<cTetra->operator[](1).getIndex()-1<<" ";
        outData<<cTetra->operator[](2).getIndex()-1<<" ";
        outData<<cTetra->operator[](3).getIndex()-1<<"\n";
        cellNext= cellNext->readNext();
    }
    outData<<"\nCELL_TYPES ";
    outData<<cellsNumber<<"\n";
    for (int i=0; i<cellsNumber; i++)
    {
        //outData<<5<<"\n";//this number indicates the cell type of each cell. 5 is the code number for triangle
        outData<<10<<"\n";//this number indicates the cell type of each cell. 10 is the code number for tetrahedron
    }
    outData<<"CELL_DATA "<<cellsNumber<<"\n";
    outData<<"SCALARS cell_Quality float "<<1<<"\n";
    outData<<"LOOKUP_TABLE default\n";

    cellNext= &m().first();
    for (int i=0; i<cellsNumber; i++)
    {
        cTetra= cellNext->getItem(0);
        outData<<cTetra->cellQuality()<<"\n";
        cellNext= cellNext->readNext();
    }
    
    outData.close();
}

#endif /* vtkUnstrMeshTranslator_h */
