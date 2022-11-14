//
//  createMesh.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 1/12/21.
//

#ifndef createMesh_h
#define createMesh_h

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include "omp.h"

#include "list.h"
#include "linkedList.h"
#include "node.h"
#include "edge.h"
#include "cell.h"
#include "mesh.h"
#include "face.h"

template<typename T, int M, int N>
//T: specifies the arithmetic type used, e.g. int, float, double...
//M: T specifies the number of nodes in the cell (e.g. for a triangle use M=3).
//N: specifies the dimensionality of the nodes used in the cell. They can be 1 (1-D), 2 (2-D), 3 (3-D).
class createMesh{
public:
    createMesh(const std::string filename,const std::string pathname=""); //constructor
    ~createMesh();//destructor
    
    const mesh<T,M,N> &operator() () const; //returns the mesh that has been initialized from the files
    void findEdgesLogFr();
    void findFacesLogFr();
    
    int getNodesNumber() const;
    int getInnerBoundaryNodesNumber() const;
    int getOuterBoundaryNodesNumber() const;
    int getEdgesNumber() const;
    int getInnerBoundaryEdgesNumber() const;
    int getOuterBoundaryEdgesNumber() const;
    int getFacesNumber() const;
    int getInnerBoundaryFacesNumber() const;
    int getOuterBoundaryFacesNumber() const;
    int getCellsNumber() const;

    template<typename t, int m, int n> friend class meshDisplacement;
protected:
    mesh<T,M,N> cMesh;
    
    int nodesNumber; //total number of nodes in the grid
    int edgesNumber; //total number of edges in the grid
    int facesNumber; //total number of faces in the grid
    int cellsNumber; //total number of cells in the grid
    int iBNodesNumber; //total number of inner boundary edges in the grid
    int oBNodesNumber; //total number of outer boundary edges in the grid
    int iBEdgesNumber; //total number of inner boundary edges in the grid
    int oBEdgesNumber; //total number of outer boundary edges in the grid
    int iBFacesNumber; //total number of inner boundary faces in the grid
    int oBFacesNumber; //total number of outer boundary faces in the grid
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Definitions of functions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//constructor for 3-D mesh consisting of tetrahedron cells
template<>
createMesh<double,4,3>::createMesh(const std::string filename,const std::string pathname)
{
    std::ifstream inData;
    std::string filename1= pathname + filename + ".nod";
    std::string filename2= pathname + filename + ".ele";
    int *fileIntElement= new int[4]{-1}; /*variables to temporarily store each integer element of the files.
                                          4 variables are needed to store each vertex of a tetrahedron.*/
    int inverted=0, shCellsLength, cFacesLength, cEdgesLength, faceNo, edgeNo, count=-1;;
    double fileDoubleElement;
    vector<double> edgeNormalVector(3),node1NormalVector(3),node2NormalVector(3);
    node<double,3> *cNode[3]{NULL},*nod=NULL;
    edge<double,3> *edgeElement[6]{NULL}, *cEdge[6]{NULL};
    face<double,3> *faceElement[4]{NULL}, *cFace[6]{NULL};
    tetrahedron *cTetra=NULL, *cell2=NULL;
    const linkedList<tetrahedron *> *cTetraNext;
    node<double,3> **nodes;
    bool flag1,flag2, condition1, condition2;
    
    omp_set_max_active_levels(2);
    //open file .nod to construct node objects.
    inData.open(filename1,std::ios::in);
    if (!inData)
    {
        std::cout<<"Error: No such file"<<std::endl;
        std::cout<<"Is this filepath and filename correct?\n"<<filename1<<"\n";
        exit(0);
    }
    inData.seekg(0,std::ios::beg); //set runner to read form the start of the open file
    
    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();
    
    std::cout<<"Initializing nodes...\n";
    inData>>nodesNumber;
    nodes= new node<double,3> *[nodesNumber];
    facesNumber=0;
    edgesNumber=0;
    iBNodesNumber=0;
    oBNodesNumber=0;
    while (inData >> fileDoubleElement)
    /* current while-loop reads the elements of each line of the .nod file one by one and initializes the nodes. It stops once it reaches the end of the file.
     Format of .nod file is:
        - First Element is an integer number,lets call it N, which is the total number of nodes in the grid (it has been extracted
          and saved at variable nodesNumber above).
        - The next N elements correspond to the logFr value (which is an integer) of each node.
        - The next N elements correspond to the x-component (which is a double)  of each node.
        - The last N elements correspond to the y-component (which is a double) of each node.
     The index each node takes corresponds to the position its elements are appearing in the .nod file, e.g. node with index 1 is the node that its logFr, x-component and y-component appear at the first position relative to the logFr, x-component and y-component of the rest of the nodes. That node is saved in the node[0] position of the nodes array.
     NOTE: Each node is saved in the list position [index-1], e.g node with index 1 is saved in list's position 0 (nodesList(0)), node with index 2 is saved in list's position 1 (nodesList(1)). This happens just because in C/C++, arrays are numbered starting from 0, while in FORTRAN (default language for using .NOD files) the default numbering starts from 1. */
    {
        count+=1;
        if (count < nodesNumber)
        {
            nod= new node<double,3>;
            nod->setIndex(count+1);
            nod->setLogFr((int)fileDoubleElement);
            if ((int)fileDoubleElement == 3)
                iBNodesNumber++;
            else if ((int)fileDoubleElement == 4)
                oBNodesNumber++;
            nodes[count]= nod;
        }
        else if (count < 2*nodesNumber) //read x-component from file and set node's x component equal to it
            nodes[count-nodesNumber]->set(0,fileDoubleElement);
        else if (count < 3*nodesNumber) //read y-component from file and set node's x component equal to it
            nodes[count-2*nodesNumber]->set(1,fileDoubleElement);
        else if (count < 4*nodesNumber) //read z-component from file and set node's x component equal to it
            nodes[count-3*nodesNumber]->set(2,fileDoubleElement);
    }
    inData.close();
    //end of reading .nod file
    
    std::cout<<"Creating cells...\n";
    //open .ele file to construct cell objects//
    inData.open(filename2,std::ios::in);
    if (!inData)
    {
        std::cout<<"Error: No such file"<<std::endl;
        std::cout<<"Is this filepath and filename correct?\n"<<filename2<<"\n";
        exit(0);
    }
    inData.seekg(0,std::ios::beg); //set runner to read form the start of the open file
    inData>>cellsNumber;
    count=0;
    while (inData>>fileIntElement[0])
    /*
     Current while loop reads the elements of the .ele file which has the index of the three nodes that
     correspond to the three vertexes of each triangle. Each triangle is constructed and then is inserted
     into the grid which is a mesh object (linked list that has objects of type cell). At the same time a
     list containing the indexes of the triangles, that each node is a part of is created. Also a linked
     list containing all the edges of the mesh, a linked list containing all the boundary edges and a linked
     list containing all the boundary nodes are created.
     */
    {
        count+=1;
        //starts with the creation of a triangle as it reads it from .ele file
        inData>>fileIntElement[1];
        inData>>fileIntElement[2];
        inData>>fileIntElement[3];
        cTetra= new tetrahedron;
        cTetra->putNode(0,nodes[fileIntElement[0]-1],1,1);
        cTetra->putNode(1,nodes[fileIntElement[2]-1],1,1);
        cTetra->putNode(2,nodes[fileIntElement[1]-1],1,1);
        cTetra->putNode(3,nodes[fileIntElement[3]-1],1,1);
        cTetra->setIndex(count);
        cTetra->cellVolume();
        if (cTetra->getVolume() < 0)/*checks if area of tetrahedron <0. if so, then it interchanges the vertexes 1,2
                                     of the tetrahedron and by doing so the area is >0.*/
        {
            cTetra->putNode(1,nodes[fileIntElement[1]-1],0,0);
            cTetra->putNode(2,nodes[fileIntElement[2]-1],0,0);
            inverted+=1;
            cTetra->cellVolume();//recalculates the area.
        }
        if (cTetra->getVolume() < pow(10,-10))
        {
            std::cout<<"Error: Vanished triangle. Volume of tetrahedron "<<cTetra->getIndex()<<" is too small."<<std::endl;
            exit(0);
        }
        //end of creation of the triangle and isertion of it into the mesh.
        cMesh.insertNextItem(cTetra);
                                            //complementary analysis//
        // saves all the faces of the cell. a face of a tetrahedron is described by three nodes.
        faceNo=0;
        for (int vertexIndex=0; vertexIndex<2; vertexIndex++)
            for (int vertexIndex2=vertexIndex+1; vertexIndex2<3; vertexIndex2++)
                for (int vertexIndex3=vertexIndex2+1; vertexIndex3<4; vertexIndex3++) //forLoop1
                {
                    cNode[0]= nodes[fileIntElement[vertexIndex]-1];
                    cNode[1]= nodes[fileIntElement[vertexIndex2]-1];
                    cNode[2]= nodes[fileIntElement[vertexIndex3]-1];
                    faceElement[faceNo]= new face<double,3>(cNode[0],cNode[1],cNode[2]);
                    faceNo++;
                }
        #pragma omp parallel for num_threads(faceNo) default(none) private(flag1,flag2,cNode,shCellsLength,cell2,cFacesLength) shared(cTetra,faceElement,cFace,facesNumber)
        {
            for (int iF=0; iF<4; iF++)
            {
                flag1=1;
                flag2=0;

                cNode[0]= faceElement[iF]->getFaceNode(0);
                cNode[1]= faceElement[iF]->getFaceNode(1);
                cNode[2]= faceElement[iF]->getFaceNode(2);
                //search each node's neighboring cells to find if any of its faces is the same with the faceElement.
                //A face can only belong to two cells, so if a second cell is found the loop breaks.
                for (int i=0; i<3; i++) //forLoop2
                {
                    if (cNode[i]->getSharingCellNextPointer(0) == NULL)
                        continue;
                    shCellsLength= cNode[i]->getSharingCells();
                    for (int j=1; j<shCellsLength+1; j++) //forLoop3
                    {
                        cell2= (tetrahedron *)cNode[i]->getSharingCellAddress(j);
                        if (cell2->getIndex() == cTetra->getIndex())
                            continue;
                        cFace[iF]= cell2->getCellFace(1);
                        if (cell2 != NULL && cFace[iF] != NULL)
                        {
                            cFacesLength= cell2->cellFacesLength();
                            for (int k=1; k<cFacesLength; k++)
                            {
                                cFace[iF]= cell2->getCellFace(k);
                                if (*cFace[iF] == *faceElement[iF])
                                {
                                    #pragma omp critical
                                    {
                                        cTetra->insertCellFace(cFace[iF]);
                                    }
                                    flag1=0;
                                    flag2=1;
                                    break;
                                }
                            }
                        }
                    } //end of forLoop3
                    if (flag2)
                        break;
                } //end of forLoop2
                if (flag1)
                {
                    faceElement[iF]->findArea();
                    faceElement[iF]->findNormalVector();
                    faceElement[iF]->setFaceLogFr(0); //initialize it to 0, a.k.a. inner face. if boundary face, gonna be rewritten.
                    #pragma omp critical
                    {
                        cTetra->insertCellFace(faceElement[iF]);
                        facesNumber++;
                        faceElement[iF]->setIndex(facesNumber);
                    }
                    faceElement[iF]=NULL;
                }
            }
        }
        for (int i=0; i<faceNo; i++)
            delete faceElement[i];

        // saves all the edges of the cell. an edge is described by the two nodes that create it.
        edgeNo=0;
        for (int vertexIndex=0; vertexIndex<3; vertexIndex++)
            for (int vertexIndex2=vertexIndex+1; vertexIndex2<4; vertexIndex2++)
            {
                cNode[0]= nodes[fileIntElement[vertexIndex]-1];
                cNode[1]= nodes[fileIntElement[vertexIndex2]-1];
                edgeElement[edgeNo]= new edge<double,3>(cNode[0],cNode[1]);
                edgeNo++;
            }
        #pragma omp parallel for num_threads(edgeNo) default(none) private(flag1, cNode, cell2, shCellsLength, condition1, condition2, cEdgesLength) shared(cTetra, cEdge, cFace, edgeElement, edgesNumber)
        { //start of parallel region
            for (int iE=0; iE<6; iE++)
            { //forLoop1
                cNode[0]= edgeElement[iE]->getNode(0);
                cNode[1]= edgeElement[iE]->getNode(1);
                //search each node's neighboring cells to find if any of its edges is the same with the edgeElement.
                flag1=1;
                for (int i=0; i<2; i++)
                { //forLoop2
                    if (cNode[i]->getSharingCellNextPointer(0) == NULL)
                        continue;
                    shCellsLength= cNode[i]->getSharingCells();
                    for (int j=1; j<shCellsLength+1; j++)
                    { //forLoop3
                        cell2= (tetrahedron *)cNode[i]->getSharingCellAddress(j);
                        if (cell2->getIndex() == cTetra->getIndex())
                            continue;
                        cEdge[iE]= cell2->getCellEdge(1);
                        if (cell2 != NULL && cEdge[iE] != NULL)
                        {
                            cEdgesLength= cell2->cellEdgesLength();
                            for (int k=1; k<cEdgesLength; k++)
                            {
                                cEdge[iE]= cell2->getCellEdge(k);
                                if (*cEdge[iE] == *edgeElement[iE])
                                {
                                    #pragma omp critical
                                    {
                                        cTetra->insertCellEdge(cEdge[iE]);
                                    }
                                    //find this edge at which of the 4 faces of the cell belongs to
                                    //if number of cores available >= edgeNo*faceNo then set num_threads equal to faceNo. otherwise, possible slow down
                                    //if so, declare faceNo as shared variable
                                    #pragma omp parallel for num_threads(1) default(none) private(condition1, condition2, cFace) shared(cTetra, cEdge, iE)
                                    {
                                        for (int k=1; k<5; k++)
                                        {
                                            cFace[iE]= cTetra->getCellFace(k);
                                            condition1= cEdge[iE]->getNode(0) == cFace[iE]->getFaceNode(0) || cEdge[iE]->getNode(0) == cFace[iE]->getFaceNode(1) || cEdge[iE]->getNode(0) == cFace[iE]->getFaceNode(2);
                                            condition2= cEdge[iE]->getNode(1) == cFace[iE]->getFaceNode(0) || cEdge[iE]->getNode(1) == cFace[iE]->getFaceNode(1) || cEdge[iE]->getNode(1) == cFace[iE]->getFaceNode(2);
                                            if ( condition1 && condition2)
                                            {
                                                #pragma omp critical
                                                {
                                                    cEdge[iE]->insertEdgeFace(cFace[iE]);
                                                }
                                            }
                                        }
                                    }
                                    flag1=0;
                                    break;
                                }
                            }
                        }
                        if (!flag1)
                            break;
                    } //end of forLoop3
                    if (!flag1)
                        break;
                } //end of forLoop2
                if (flag1)
                {
                    edgeElement[iE]->computeDirectionVector();
                    #pragma omp critical
                    {
                        edgesNumber++;
                        edgeElement[iE]->setIndex(edgesNumber);
                        cTetra->insertCellEdge(edgeElement[iE]);
                        cNode[0]->insertNeighborNode(cNode[1]);
                        cNode[1]->insertNeighborNode(cNode[0]);
                    }
                    //find at which of the 4 faces of the cell the edge belongs to
                    //if number of cores available >= edgeNo*faceNo then set num_threads equal to faceNo. otherwise, possible slow down
                    //if so, declare faceNo as shared variable
                    #pragma omp parallel for num_threads(1) default(none) private(condition1, condition2, cFace) shared(cTetra, edgeElement, iE)
                    {
                        for (int k=1; k<5; k++)
                        {
                            cFace[iE]= cTetra->getCellFace(k);
                            condition1= edgeElement[iE]->getNode(0) == cFace[iE]->getFaceNode(0) || edgeElement[iE]->getNode(0) == cFace[iE]->getFaceNode(1) || edgeElement[iE]->getNode(0) == cFace[iE]->getFaceNode(2);
                            condition2= edgeElement[iE]->getNode(1) == cFace[iE]->getFaceNode(0) || edgeElement[iE]->getNode(1) == cFace[iE]->getFaceNode(1) || edgeElement[iE]->getNode(1) == cFace[iE]->getFaceNode(2);
                            
                            if ( condition1 && condition2)
                                edgeElement[iE]->insertEdgeFace(cFace[iE]);
                        }
                    }
                    edgeElement[iE]->dropFirstEdgeFace();
                    edgeElement[iE]= NULL;
                }
            } //forLoop1
        } //end of parallel region
        for (int i=0; i<edgeNo; i++)
            delete edgeElement[i];
    }
    inData.close();
    //end of reading file .ele
    cMesh.dropFirstItem();
    //drop first item of linked lists sharingCellsAddress,neighborNode CellFace, CellEdge because its a dummy entry.
    //Could have been dropped when the list was initialized before, but we would need an if statement, that would cost time.
    #pragma omp parallel for default(none) private(cNode) shared(nodesNumber, nodes)
    {
        for (int i=0; i<nodesNumber; i++)
        {
            cNode[0]= nodes[i];
            cNode[0]->dropFirstSharingCellAddress();
            cNode[0]->dropNeighborNode(0);
        }
    }
    cTetraNext= &cMesh.first();
    for (int i=0; i<cellsNumber; i++)
    {
        cTetra= cTetraNext->getItem(0);
        cTetra->dropFirstCellEdge();
        cTetra->dropFirstCellFace();
        cTetraNext= cTetraNext->readNext();
    }
    findFacesLogFr();
    findEdgesLogFr();
    
    #pragma omp parallel for default(none) shared(nodesNumber, nodes)
    {
        for (int i=0; i<nodesNumber; i++)
            nodes[i]=NULL;
    }
    delete [] nodes;
    
    if (inverted>0)
        std::cout<<"On the initial grid existed "<<inverted<<" inverted cells.\n";

    std::cout<<"\n\n\n";
    std::cout<<"Number of Nodes:                   "<<nodesNumber<<"\n";
    std::cout<<"Number of Edges:                   "<<edgesNumber<<"\n";
    std::cout<<"Number of faces:                   "<<facesNumber<<"\n";
    std::cout<<"Number of Cells:                   "<<cellsNumber<<"\n";
    std::cout<<"Number of Boundary Nodes:          "<<iBNodesNumber+oBNodesNumber<<"\n";
    std::cout<<"Number of Inner Boundary Nodes:    "<<iBNodesNumber<<"\n";
    std::cout<<"Number of Outer Boundary Nodes:    "<<oBNodesNumber<<"\n";
    std::cout<<"Number of Boundary Edges:          "<<iBEdgesNumber+oBEdgesNumber<<"\n";
    std::cout<<"Number of Inner Boundary Edges:    "<<iBEdgesNumber<<"\n";
    std::cout<<"Number of Outer Boundary Edges:    "<<oBEdgesNumber<<"\n";
    std::cout<<"Number of Boundary Faces:          "<<iBFacesNumber+oBFacesNumber<<"\n";
    std::cout<<"Number of Inner Boundary Faces:    "<<iBFacesNumber<<"\n";
    std::cout<<"Number of Outer Boundary Faces:    "<<oBFacesNumber<<"\n";
    std::cout<<"\n\n\n";

    /*
     int EulerType= nodesNumber - edgesNumber + facesNumber - cellsNumber
                    - 0.5*(iBNodesNumber + oBNodesNumber) + 0.25*(iBFacesNumber + oBFacesNumber);
    if (EulerType != 0)
    {
        std::cout<<"Euler's type = "<<EulerType<<"\n";
        std::cout<<"Some edges or faces are missing from your list(s). Logical debug necessary.\n";
        exit(0);
    }
     */
    // Stop measuring time and calculate the elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
}


//destructor
template<typename T, int M, int N>
createMesh<T,M,N>::~createMesh(){
}

//returns the mesh that has been initialized by the files
template<typename T,int M,int N>
inline const mesh<T,M,N> &createMesh<T,M,N>::operator() () const{
    return cMesh;
}

template<typename T,int M,int N>
inline int createMesh<T,M,N>::getNodesNumber() const{
    return nodesNumber;
}

template<typename T,int M,int N>
inline int createMesh<T,M,N>::getInnerBoundaryNodesNumber() const{
    return iBNodesNumber;
}

template<typename T,int M,int N>
inline int createMesh<T,M,N>::getOuterBoundaryNodesNumber() const{
    return oBNodesNumber;
}

template<typename T,int M,int N>
inline int createMesh<T,M,N>::getEdgesNumber() const{
    return edgesNumber;
}

template<typename T,int M,int N>
inline int createMesh<T,M,N>::getInnerBoundaryEdgesNumber() const{
    return iBEdgesNumber;
}

template<typename T,int M,int N>
inline int createMesh<T,M,N>::getOuterBoundaryEdgesNumber() const{
    return oBEdgesNumber;
}

template<typename T,int M,int N>
inline int createMesh<T,M,N>::getFacesNumber() const{
    return facesNumber;
}

template<typename T,int M,int N>
inline int createMesh<T,M,N>::getInnerBoundaryFacesNumber() const{
    return iBFacesNumber;
}

template<typename T,int M,int N>
inline int createMesh<T,M,N>::getOuterBoundaryFacesNumber() const{
    return oBFacesNumber;
}


template<typename T,int M,int N>
inline int createMesh<T,M,N>::getCellsNumber() const{
    return cellsNumber;
}

template<typename T, int M, int N>
void createMesh<T,M,N>::findFacesLogFr(){
    cell<T,M,N> *cCell;
    face<T,N> *cFace;
    const linkedList<cell<T,M,N> *> *cellNext;
    int fNum;
    
    //check boundary faces to drop faces that have been saved in list but are not in the boundary
    //only way for a face to belong to the boundary is if it belongs to only one cell.
    std::cout<<"Find boundary faces logFr...\n";
    iBFacesNumber=0;
    oBFacesNumber=0;
    cellNext= &cMesh.first();
    for (int i=0; i<cellsNumber; i++)
    {
        cCell= cellNext->getItem(0);
        cellNext= cellNext->readNext();
        fNum= cCell->cellFacesLength();
        #pragma omp parallel for num_threads(fNum) default(none) private(cFace) shared(fNum, cCell)
        {
            for (int j=0; j<fNum; j++)
            {
                cFace=cCell->getCellFace(j);
                if (cFace->getFlag())
                    continue;
                
                cFace->dropFirstSharingCell();
                
                if (cFace->noBoundaryFace())
                {
                    cFace->setFaceLogFr(0);
                    cFace->setFlagTrue();
                }
                else
                    if (cFace->getFaceNode(0)->getLogFr() == 3)
                    {
                        cFace->setFaceLogFr(3);
                        cFace->setFlagTrue();
                        iBFacesNumber++;
                    }
                    else
                    {
                        cFace->setFaceLogFr(4);
                        cFace->setFlagTrue();
                        oBFacesNumber++;
                    }
            }
        }
    }
    cellNext= &cMesh.first();
    for (int i=0; i<cellsNumber; i++)
    {
        cCell= cellNext->getItem(0);
        cellNext= cellNext->readNext();
        fNum= cCell->cellFacesLength();
        #pragma omp parallel for num_threads(fNum) default(none) private(cFace) shared(fNum, cCell)
        {
            for (int j=0; j<fNum; j++)
            {
                cFace=cCell->getCellFace(j);
                if (!cFace->getFlag())
                    continue;
                
                cFace->setFlagFalse();
            }
        }
    }
}

template<typename T, int M, int N>
void createMesh<T,M,N>::findEdgesLogFr(){
    cell<T,M,N> *cCell;
    edge<T,N> *cEdge;
    face<T,N> *cFace;
    const linkedList<cell<T,M,N> *> *cellNext;
    int eNum, cFacesNum;
    
    std::cout<<"Find boundary edges logFr...\n";
    iBEdgesNumber=0;
    oBEdgesNumber=0;
    cellNext= &cMesh.first();
    for (int i=0; i<cellsNumber; i++)
    {
        cCell= cellNext->getItem(0);
        cellNext= cellNext->readNext();
        eNum= cCell->cellEdgesLength();
        #pragma omp parallel for num_threads(eNum) default(none) private(cEdge, cFace, cFacesNum) shared(eNum, cCell)
        {
            for (int j=0; j<eNum; j++)
            {
                cEdge=cCell->getCellEdge(j);
                if (cEdge->getFlag())
                    continue;
                
                cEdge->dropFirstSharingCell();
                
                cEdge->setEdgeLogFr(0);
                cEdge->setFlagTrue();
                cFacesNum= cEdge->getEdgeFacesLength();
                for (int k=0; k<cFacesNum; k++)
                {
                    cFace= (face<double,3> *)cEdge->getFace(k);
                    if (cFace->getFaceLogFr() == 3)
                    {
                        cEdge->setEdgeLogFr(3);
                        iBEdgesNumber++;
                        break;
                    }
                    else if (cFace->getFaceLogFr() == 4)
                    {
                        cEdge->setEdgeLogFr(4);
                        oBEdgesNumber++;
                        break;
                    }
                }
            }
        }
    }
    
    cellNext= &cMesh.first();
    for (int i=0; i<cellsNumber; i++)
    {
        cCell= cellNext->getItem(0);
        cellNext= cellNext->readNext();
        eNum= cCell->cellEdgesLength();
        #pragma omp parallel for num_threads(eNum) default(none) private(cEdge, cFace, cFacesNum) shared(eNum, cCell)
        {
            for (int j=0; j<eNum; j++)
            {
                cEdge=cCell->getCellEdge(j);
                if (!cEdge->getFlag())
                    continue;
                
                cEdge->setFlagFalse();
            }
        }
    }
}
#endif /* createMesh_h */
