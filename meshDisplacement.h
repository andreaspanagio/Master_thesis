//
//  meshDisplacement.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 1/1/22.
//

#ifndef meshDisplacement_h
#define meshDisplacement_h

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "omp.h"

#include "vector.h"
#include "node.h"
#include "cell.h"
#include "face.h"
#include "createMesh.h"
#include "linkedList.h"
#include "quaternion.h"

template<typename T, int M, int N>
//T: specifies the arithmetic type used, e.g. int, float, double...
//M: T specifies the number of nodes in the cell (e.g. for a triangle use M=3).
//N: specifies the dimensionality of the nodes used in the cell. They can be 1 (1-D), 2 (2-D), 3 (3-D).
class meshDisplacement{
public:
    meshDisplacement(createMesh<T,M,N> &m,const std::string pathname="");
    ~meshDisplacement();
    void calculateDisplacement(node<T,N> **nodes, int nodesNum);
    
    //adaptive mesh deformation techniques available
    void IDW(bool slNodes);
    void IDW(node<T,N> **inNodes, node<T,N> **bNodes, int inNodesNum, int bNodesNum, T cLength, T expo1, T expo2, T acap, T alfa);
    void IDWquaternion(bool slNodes);
    void IDWquaternion(node<T,N> **inNodes, node<T,N> **bNodes, face<T,N> **bFaces, int inNodesNum, int bNodesNum, int bFacesNum, T cLength, T expo1, T expo2, T acap, T alfa);
    
    void findQuaternions(node<T,N> **bNodes, face<T,N> **bFaces, int bNodesNum, int bFacesNum);
    void findTranslations(node<T,N> **bNodes, int bNodesNum);
    T domainlength(node<T,N> **bNodes, int bNodesNum);
    T findAlfa(node<T,N> **bNodes, int bNodesNum, T cLength, T eta);
    T findAlfa(node<T,N> **iBNodes, node<T,N> **oBNodes, int iBNodesNum, int oBNodesNum, T cLength, T eta);
    void extractNodes(createMesh<T,M,N> &m);
    void extractFaces(createMesh<T,M,N> &m);
    void updateMesh(createMesh<T,M,N> &m); //updates the positions of the nodes, the cells volume etc.
    vector<T> *getDisplace() const;
    
    void slidingNodes(T cLength, T expo1, T expo2, T acap, T alfa,int method);
    vector<T> findCrossPoint(face<T,N> &f, vector<T> &vec);
    void findMaxCoords(node<T,N> **bNodes, int bNodesNum, T &xMin,T &xMax,T &yMin,T &yMax,T &zMin,T &zMax);
    void findBoundaryTurningPoints(face<T,N> **bFaces, int bFacesNum);
    
    void writeDisplacedMesh(createMesh<T,M,N> &m, const std::string filename, const std::string pathname="");
    void writeNodFile(createMesh<T,M,N> &m, const std::string filename);
private:
    node<T,N> **inNodes; //innerNodes
    node<T,N> **iBNodes; //innerBoundaryNodes
    node<T,N> **oBNodes; //outerBoundaryNodes
    face<T,N> **iBFaces; //innerBoundaryFaces
    face<T,N> **oBFaces; //outerBoundaryFaces
    int inNodesNum;
    int iBNodesNum;
    int oBNodesNum;
    int iBFacesNum;
    int oBFacesNum;
    
    linkedList<node<T,N> *> oBTNodes; //linked list of the boundary nodes where the orientation of the outer boundary changes
    list<linkedList<face<T,N> *>> oBTFaces; //list containing the adresses of the faces that turning points belong to
    int oBTNodesNum; //total number of outer boundary turning points in the grid
    
    vector<T> *displace;
    vector<T> *translation; //used in the IDWquaternion method
    quaternion<T> *rotQ; //used in the IDWquaternion method
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Definitions of functions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//constructor for 3D mesh consisting of tetrahedrons
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
meshDisplacement<double,4,3>::meshDisplacement(createMesh<double,4,3> &m, const std::string pathname):translation(NULL), rotQ(NULL), oBTNodes(NULL), oBTFaces(){
    int moveIt;
    
    inNodesNum= m.nodesNumber - m.iBNodesNumber - m.oBNodesNumber;
    iBNodesNum= m.iBNodesNumber;
    oBNodesNum= m.oBNodesNumber;
    extractNodes(m);
    iBFacesNum= m.iBFacesNumber;
    oBFacesNum= m.oBFacesNumber;
    extractFaces(m);
    
    displace= new vector<double>[m.nodesNumber];
    //initialization of translation
    #pragma omp parallel for default(none) shared(translation,m)
    {
        for (int i=0; i<m.nodesNumber; i++)
            displace[i].setDimension(3);
    }
    
    while(1)
    {
        std::cout<<"Enter logFr of nodes to be displaced (or -1 to stop displacing nodes): ";
        std::cin>>moveIt;
    
        std::cout<<"\nFor nodes with logFr = "<<moveIt;
        std::cout<<"\n------------------------\n";
        if (moveIt == 3)
            calculateDisplacement(iBNodes,iBNodesNum);
        else if (moveIt == 4)
            calculateDisplacement(oBNodes,oBNodesNum);
        else
            return;
    } //end of while loop
}

//destructor
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
meshDisplacement<T,M,N>::~meshDisplacement(){
    delete [] displace;
    delete [] inNodes;
    delete [] iBNodes;
    delete [] oBNodes;
    delete [] iBFaces;
    delete [] oBFaces;
    if (translation)
        delete [] translation;
    if (rotQ)
        delete [] rotQ;
}

//IDW
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
void meshDisplacement<T,M,N>::IDW(bool slNodes){
    T cLength, expo1, expo2, alfa, acap, eta, sumW, dist, term1, term2, wme;
    vector<T> sumD(3), dr(3);
    node<T,N> *node1, *node2;
    
    //find characteristic length
    cLength= domainlength(oBNodes, oBNodesNum);
    
    expo1=3;
    expo2=5;
    acap=1.;
    eta=5.;
    alfa= findAlfa(iBNodes,iBNodesNum,cLength,eta);
    
    if (slNodes)
    {
        slidingNodes(cLength,expo1,expo2,acap,alfa,1);
        //recalculate alpha cause now outer boundary nodes have moved as well
        alfa=findAlfa(iBNodes,oBNodes,iBNodesNum,oBNodesNum,cLength,eta);
    }
    
    std::cout<<"Displacing inner nodes...\n";
    //find displace for inNodes
    #pragma omp parallel for default(none) private(node1, node2, sumW, dist, term1, term2, wme)  firstprivate(sumD, dr) shared(displace, inNodesNum, iBNodesNum, oBNodesNum, cLength, alfa, acap, expo1, expo2, inNodes, iBNodes, oBNodes)
    { //start of parallel region
        for (int i=0; i<inNodesNum; i++)
        {
            node1= inNodes[i];
            sumD.set(1,0.);
            sumD.set(2,0.);
            sumD.set(3,0.);
            sumW=0.;
            for (int j=0; j<iBNodesNum; j++)
            {
                node2= iBNodes[j];
                dr= node2->operator()()- node1->operator()();
                dist= dr.l2norm();
                term1= cLength/dist;
                term2= alfa*term1;
                wme= acap * (pow(term1,expo1) + pow(term2,expo2));
                sumD+= wme*displace[node2->getIndex()-1];
                sumW+= wme;
            }
            for (int j=0; j<oBNodesNum; j++)
            {
                node2= oBNodes[j];
                dr= node2->operator()()- node1->operator()();
                dist= dr.l2norm();
                term1= cLength/dist;
                term2= alfa*term1;
                wme= acap * (pow(term1,expo1) + pow(term2,expo2));
                sumD+= wme*displace[node2->getIndex()-1];
                sumW+= wme;
            }
            displace[node1->getIndex()-1]= sumD/sumW;
        }
    } //end of parallel region
}

template<typename T, int M, int N>
void meshDisplacement<T,M,N>::IDW(node<T,N> **inNodes, node<T,N> **bNodes,int inNodesNum,int bNodesNum, T cLength, T expo1, T expo2, T acap, T alfa){
    T sumW, dist, term1, term2, wme;
    vector<T> sumD(3), dr(3);
    node<T,N> *node1, *node2;
    
    //find displace for inNodes
    #pragma omp parallel for default(none) private(node1, node2, sumW, dist, term1, term2, wme)  firstprivate(sumD, dr) shared(displace, inNodesNum, bNodesNum, cLength, alfa, acap, expo1, expo2, inNodes, bNodes)
    { //start of parallel region
        for (int i=0; i<inNodesNum; i++)
        {
            node1= inNodes[i];
            sumD.set(1,0.);
            sumD.set(2,0.);
            sumD.set(3,0.);
            sumW=0.;
            for (int j=0; j<bNodesNum; j++)
            {
                node2= bNodes[j];
                dr= node2->operator()()- node1->operator()();
                dist= dr.l2norm();
                term1= cLength/dist;
                term2= alfa*term1;
                wme= acap * (pow(term1,expo1) + pow(term2,expo2));
                sumD+= wme*displace[node2->getIndex()-1];
                sumW+= wme;
            }
            displace[node1->getIndex()-1]= sumD/sumW;
        }
    } //end of parallel region
}

//IDW-Quaternion
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
void meshDisplacement<T,M,N>::IDWquaternion(bool slNodes){
    T cLength, expo1, expo2, alfa, acap, eta, sumW, dist, term1, term2, wme;
    node<T,N> *node1, *node2;
    quaternion<T> sumR, R1;
    vector<T> dr(3), sumD(3), r1(3);
    int oBSLNodesNum=0, nodesNum= inNodesNum + iBNodesNum + oBNodesNum;
    
    rotQ= new quaternion<double>[nodesNum];
    translation= new vector<double>[nodesNum];
    //initialization of translation
    for (int i=0; i<nodesNum; i++)
        translation[i].setDimension(3);

    //find characteristic length
    cLength= domainlength(oBNodes, oBNodesNum);
    //enter values for IDW parameters
    expo1=3;
    expo2=5;
    acap=1.;
    eta=5.;
    alfa= findAlfa(iBNodes,iBNodesNum,cLength,eta);
    //analyze displace into a rotation and a translation
    findQuaternions(iBNodes, iBFaces, iBNodesNum, iBFacesNum);
    findTranslations(iBNodes, iBNodesNum);
    if (slNodes)
    {
        slidingNodes(cLength,expo1,expo2,acap,alfa,2);
        oBSLNodesNum= 1;
        findQuaternions(oBNodes,oBFaces, oBNodesNum, oBFacesNum);
        findTranslations(oBNodes,oBNodesNum);
        //recalculate alpha cause now outer boundary nodes have moved as well
        alfa=findAlfa(iBNodes,oBNodes,iBNodesNum,oBNodesNum,cLength,eta);
    }
    
    std::cout<<"Displacing inner nodes...\n";
    //IDW-quaternion for displacing internal nodes.
    #pragma omp parallel for default(none) private(node1, node2, sumR, sumW, dist, term1, term2, wme, R1)  firstprivate(sumD, dr, r1) shared(displace, inNodesNum, iBNodesNum, oBNodesNum, cLength, alfa, acap, expo1, expo2, inNodes, iBNodes, oBNodes, translation, rotQ, oBSLNodesNum)
    { //start of parallel region
        for (int i=0; i<inNodesNum; i++)
        {
            node1= inNodes[i];
            sumD.set(0,0.);
            sumD.set(1,0.);
            sumD.set(2,0.);
            sumR.setToZero();
            sumW=0.;
            for (int j=0; j<iBNodesNum; j++)
            {
                node2= iBNodes[j];
                dr= node2->operator()() - node1->operator()();
                dist= dr.l2norm();
                term1= cLength/dist;
                term2= alfa*term1;
                wme= acap * (pow(term1,expo1) + pow(term2,expo2));
                sumD+= wme*translation[node2->getIndex()-1];
                sumR+= wme*rotQ[node2->getIndex()-1];
                sumW+= wme;
            }
            for (int j=0; j<oBNodesNum; j++)
            {
                node2= oBNodes[j];
                dr= node2->operator()() - node1->operator()();
                dist= dr.l2norm();
                term1= cLength/dist;
                term2= alfa*term1;
                wme= acap * (pow(term1,expo1) + pow(term2,expo2));
                for (int k=0; k<oBSLNodesNum; k++)
                {
                    //in case no sliding nodes are used oBSLNodesNum==0, so this loop is ignored. This is done
                    //to save time cause both oBTranslation and oBRotQ are zero. In case sliding nodes are used
                    //oBSLNodesNum==1 so the computations are performed one time.
                    sumD+= wme*translation[node2->getIndex()-1];
                    sumR+= wme*rotQ[node2->getIndex()-1];
                }
                sumW+= wme;
            }
            //Rotation
            sumR= sumR/sumW;
            sumR= sumR/sumR.l2norm(); //R is not a unit quaternion, so it has to be normalized
            R1.set(0,0.);
            R1.set(1,node1->operator[](0));
            R1.set(2,node1->operator[](1));
            R1.set(3,node1->operator[](2));
            R1= sumR*R1*sumR.conj();
            r1.set(0,R1.get(1));
            r1.set(1,R1.get(2));
            r1.set(2,R1.get(3));
            displace[node1->getIndex()-1]= sumD/sumW  + r1 - node1->operator()();
        }
    } //end of parallel region
}

template<typename T, int M, int N>
void meshDisplacement<T,M,N>::IDWquaternion(node<T,N> **inNodes,node<T,N> **bNodes, face<T,N> **bFaces, int inNodesNum,int bNodesNum, int bFacesNum, T cLength, T expo1, T expo2, T acap, T alfa){
    T sumW, dist, term1, term2, wme;
    node<T,N> *node1, *node2;
    quaternion<T> sumR, R1;
    vector<T> dr(3), sumD(3), r1(3);
    
    //IDW-quaternion for displacing internal nodes.
    #pragma omp parallel for default(none) private(node1, node2, sumR, sumW, dist, term1, term2, wme, R1)  firstprivate(sumD, dr, r1) shared(displace, inNodesNum, bNodesNum, cLength, alfa, acap, expo1, expo2, inNodes, bNodes, translation, rotQ)
    { //start of parallel region
        for (int i=0; i<inNodesNum; i++)
        {
            node1= inNodes[i];
            sumD.set(0,0.);
            sumD.set(1,0.);
            sumD.set(2,0.);
            sumR.setToZero();
            sumW=0.;
            for (int j=0; j<bNodesNum; j++)
            {
                node2= bNodes[j];
                dr= node2->operator()()-node1->operator()();
                dist= dr.l2norm();
                term1= cLength/dist;
                term2= alfa*term1;
                wme= acap * (pow(term1,expo1) + pow(term2,expo2));
                sumD+= wme*translation[node2->getIndex()-1];
                sumR+= wme*rotQ[node2->getIndex()-1];
                sumW+= wme;
            }
            //Rotation
            sumR= sumR/sumW;
            sumR= sumR/sumR.l2norm(); //R is not a unit quaternion, so it has to be normalized
            R1.set(1,node1->operator[](0));
            R1.set(2,node1->operator[](1));
            R1.set(3,node1->operator[](2));
            R1= sumR*R1*sumR.conj();
            r1.set(0,R1.get(1));
            r1.set(1,R1.get(2));
            r1.set(2,R1.get(3));
            displace[node1->getIndex()-1]= r1 + sumD/sumW - node1->operator()();
        }
    }//end of parallel region
}

//Sliding Nodes Technique
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
void meshDisplacement<T,M,N>::slidingNodes(T cLength, T expo1, T expo2, T acap, T alfa, int method){
    int *turnPointsIndex=NULL, flag, count, bannedNodeIndLength, FNodesLength, neiCellsLength, invertedCells, neiFacesLength, bannedFaceIndLength=0, flagBan;
    T pDist, dMin, nodesDist, displ, xMin, xMax, yMin, yMax, zMin, zMax;
    quaternion<T> R, R1;
    vector<T> n(3), cPoint(3), cPointOld(3), dVec(3), dVecN(3), pF(3), nodeVecDist(3), tempVec(3), fProj(3), lastVec(3), cProd(3);
    node<T,N> *cNode=NULL, *tNode=NULL, *cFNode[3]{NULL}, *minNode[3]{NULL};
    cell<T,M,N> *cCell;
    face<T,N> *cFace=NULL, *minFace=NULL;
    linkedList<int> bannedNodeInd, bannedFaceInd;
    
    std::cout<<"Using sliding nodes..\n";
    if (method == 1) //IDW
        IDW(oBNodes, iBNodes, oBNodesNum, iBNodesNum, cLength, expo1, expo2, acap, alfa);
    else if (method == 2) //IDW-quaternions
        IDWquaternion(oBNodes, iBNodes, iBFaces, oBNodesNum, iBNodesNum, iBFacesNum, cLength, expo1, expo2, acap, alfa);
    //end of outer/inner boundary nodes displacement
    
    //find maximum coords of outer boundary
    findMaxCoords(oBNodes,oBNodesNum,xMin,xMax,yMin,yMax,zMin,zMax);
    
    //project nodes back to the boundary.
    #pragma omp parallel for default(none) private(cNode, dVec, displ, flag, dMin, cFace, flagBan, bannedFaceIndLength, bannedFaceInd, pDist, FNodesLength, pF, minFace, cPoint) shared(oBNodes, displace, oBFacesNum, oBFaces, xMin, xMax, yMin, yMax, zMin, zMax)
    {
        for (int i=0; i<oBNodesNum; i++)
        { //for_cNode
            cNode= oBNodes[i];
            dVec= displace[cNode->getIndex()-1];
            displ= dVec.l2norm();
            if (displ < pow(10,-6))
                continue;
            
            //find closest face to displaced node
            dVec= cNode->operator()()+dVec;
            
            flag=1;
            while (flag)
            {
                flag=0;
                dMin=pow(10,10);
                for (int j=0; j<oBFacesNum; j++)
                {
                    cFace= oBFaces[j];
                    flagBan=0;
                    bannedFaceIndLength= bannedFaceInd.length();
                    if (bannedFaceIndLength > 1)
                    {
                        for (int k=1; k<bannedFaceIndLength; k++)
                            if (cFace->getIndex() == bannedFaceInd(k))
                            {
                                flagBan=1;
                                break;
                            }
                    }
                    if (flagBan)
                        continue;
                    
                    //find distance from diplaced node to nodes in plane.
                    pDist=0.;
                    FNodesLength= cFace->getFaceNodesLength();
                    for (int k=0; k<FNodesLength; k++)
                    {
                        pF= cFace->getFaceNode(k)->operator()();
                        pDist+= (dVec - pF).l2norm();
                    }
                    if (pDist<dMin)
                    {
                        dMin= pDist;
                        minFace= cFace;
                    }
                }
                cFace=NULL;
                
                //find cross point between minFace and displaced cNode
                cPoint= findCrossPoint(*minFace,dVec);
                
                //check if cross point lies out of the boundary. if so ban this face
                //and find the next one that is closest to the node.
                if ( (cPoint(0) < xMin*1.1) &&  (cPoint(0) > xMax*1.1) &&
                    (cPoint(1)  < yMin*1.1) &&  (cPoint(1) > yMax*1.1) &&
                    (cPoint(2)  < zMin*1.1) &&  (cPoint(2) > zMax*1.1))
                {
                    bannedFaceInd.insertNextItem(minFace->getIndex());
                    flag=1;
                }
            }
            if (bannedFaceIndLength >1)
                for (int k=1; k<bannedFaceIndLength; k++)
                    bannedFaceInd.dropNextItem();
            
            dVec= cPoint - cNode->operator()();
            displace[cNode->getIndex()-1]= dVec;
        } //end of for_cNode
    }
    
    findBoundaryTurningPoints(oBFaces, oBFacesNum);
    //if turning points exist in outer boundary find closest nodes to them and put them at these locations
    if (oBTNodesNum > 1)
    {
        turnPointsIndex= new int[oBTNodesNum];
        for (int i=0; i<oBTNodesNum; i++)
        {//forLoop1
            tNode= oBTNodes(i);
            neiFacesLength= oBTFaces(i).length();
            invertedCells=1;
            while (invertedCells)
            {
                dMin=pow(10,10);
                for (int j=0; j<oBNodesNum; j++)
                {//forLoop2 (specific cNode)
                    cNode= oBNodes[j];
                    
                    flag=0;
                    bannedNodeIndLength= bannedNodeInd.length();
                    if (bannedNodeIndLength > 1)
                    {
                        for (int k=1; k<bannedNodeIndLength; k++)
                            if (cNode->getIndex() == bannedNodeInd(k))
                            {
                                flag=1;
                                break;
                            }
                    }
                    if (flag)
                        continue;
                    
                    nodesDist=0.;
                    count=0;//counts how many times a face is found
                    //find distance of all the neighboring faces nodes to the turing point.
                    neiCellsLength= cNode->getSharingCells();
                    for (int k=0; k<neiCellsLength; k++)
                    {//forLoop3 (specific cCell)
                        cCell= (cell<T,M,N> *)cNode->getSharingCellAddress(k);
                        for (int k2=0; k2<4; k2++)
                        {//forLoop4 (specific cFace)
                            cFace= cCell->getCellFace(k2);
                            //check if face is indeed an outer boundary face
                            if (cFace->getFaceLogFr() == 0)
                                continue;
                            cFNode[0]= cFace->getFaceNode(0);
                            cFNode[1]= cFace->getFaceNode(1);
                            cFNode[2]= cFace->getFaceNode(2);
                            for (int k3=0; k3<3; k3++)
                            {//forLoop5
                                dVec= displace[cFNode[k3]->getIndex()-1];
                                dVec= cFNode[k3]->operator()()+dVec;
                                //find distance from diplaced node and node in turning point;
                                nodesDist+=( dVec - tNode->operator()() ).l2norm();
                                count++;
                            }//end of forLoop5
                        }//end of forLoop4 (specific cFace)
                    }//end of forLoop3 (specific cCell)
                    //nodeDist must be normalized by number of nodes that were found.
                    nodesDist= nodesDist/count;
                    if (nodesDist<dMin)
                    {
                        dMin= nodesDist;
                        minNode[0]= cNode;
                    }
                }//end of forLoop2 (specific cNode)
                neiCellsLength= minNode[0]->getSharingCells();
                invertedCells=0;
                
                //search all the neighboring cells to check if anyone is inverted. if yes then this is not the right cell.
                for (int j=0; j<neiCellsLength; j++)
                { //forLoop specific cCell
                    cCell= (tetrahedron *)minNode[0]->getSharingCellAddress(j);
                    if (cCell->getVolume() < 0)
                    {
                        invertedCells++;
                        bannedNodeInd.insertNextItem(minNode[0]->getIndex());
                        break;
                    }
                }
            }// end of whileLoop
            turnPointsIndex[i]=minNode[0]->getIndex();
            dVec= tNode->operator()() - minNode[0]->operator()();
            displace[minNode[0]->getIndex()-1]= dVec;
            
            //delete banned nodes cause they are specific to the node checked.
            bannedNodeIndLength= bannedNodeInd.length()-1;
            for (int j=0; j<bannedNodeIndLength; j++)
                bannedNodeInd.dropNextItem();
        }//end forLoop1
        delete [] turnPointsIndex;
    }//end ifLoop1
}
//Find Quaternions of Boundary Nodes
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
void meshDisplacement<T,M,N>::findQuaternions(node<T,N> **bNodes, face<T,N> **bFaces, int bNodesNum, int bFacesNum){
    T theta1, theta2, theta3;
    quaternion<T> rotQ1, rotQ2, rotQ3, semiNewPosVecQ, perpPVecNQ;
    vector<T> semiNewPosVec(3),perpPVecN(3),perpPVecNNew(3), RotationAxis1(3), RotationAxis2(3), RotationAxis3(3), newPosVec(3), displaceVec(3), originVec(3), n(3), nNew(3);
    node<T,N> *cNode[2]{NULL}, originNode, nodeNew[2];
    face<T,N> *cFace, faceOld, faceNew;
    //check if the whole inner object has moved.
    //Initially the object was positioned, so that its "center of mass" was at the origin of the axes (0,0,0).
    //if its "center of mass" remains at the origin, then the object has not moved, but only rotated around the origin.
    //find "center of mass"
    if (bNodes[0]->getIndex() ==3)//if logFr != 3 I am finding the rotation of outer boundary which is not translated. so originVec=0.
    {
        #pragma omp parallel for default(none) private(cNode,displaceVec,newPosVec) shared(bNodesNum,bNodes) reduction(+:originVec)
        {
            for (int i=0; i<bNodesNum; i++)
            {
                cNode[0]= bNodes[i];//current node
                displaceVec= displace[cNode[0]->getIndex()-1];
                newPosVec= cNode[0]->operator()() + displaceVec;
                originVec+= newPosVec;
            }
        }
        originVec= originVec/bNodesNum;
    }
    //else originVec is initialized to zero. Originvec is zero when we find quaternions for outer boundary nodes.
    
    //initialize faceOld and faceNew
    faceOld.setFaceNode(0,&originNode);
    faceNew.setFaceNode(0,&originNode);
    for (int i=0; i<2; i++)
    {
        faceOld.insertFaceNode(&originNode);
        faceNew.insertFaceNode(&originNode);
    }
    
    //find rotation axis
    #pragma omp parallel for default(none) private(cFace, cNode, nodeNew, theta1, theta2, theta3, rotQ1, rotQ2, rotQ3, perpPVecNQ, semiNewPosVecQ, n, nNew, displaceVec, newPosVec) shared(bFacesNum, bFaces, displace, originVec, rotQ) firstprivate(faceOld, faceNew, perpPVecNNew, perpPVecN, semiNewPosVec, RotationAxis1, RotationAxis2, RotationAxis3)
    {
        for (int i=0; i<bFacesNum; i++)
        {
            cFace= bFaces[i];
            //create a face that is in between the position vectors of the boundary nodes.
            for (int j=0; j<3; j++)
            {
                cNode[0]=cFace->getFaceNode(j);
                faceOld.setFaceNode(1,cNode[0]);
                if (j == 2)
                    cNode[1]=cFace->getFaceNode(0);
                else
                    cNode[1]=cFace->getFaceNode(j+1);
                faceOld.setFaceNode(2,cNode[1]);
                
                faceOld.findNormalVector();
                n= faceOld.getNormalVector();
                //find position of boundary nodes after displacement
                for (int k=0; k<2; k++)
                {
                    displaceVec= displace[cNode[k]->getIndex()-1];
                    newPosVec= cNode[k]->operator()() + displaceVec-originVec;
                    nodeNew[k].setLocationX(newPosVec(0));
                    nodeNew[k].setLocationY(newPosVec(1));
                    nodeNew[k].setLocationZ(newPosVec(2));
                    faceNew.setFaceNode(k+1,&nodeNew[k]);
                }
                faceNew.findNormalVector();
                nNew= faceNew.getNormalVector();
                //first rotation axis. need two perpendicular axis to specify 3d rotation axis.
                RotationAxis1= RotationAxis1.cross_product(n, nNew);
                if (RotationAxis1.l2norm() >pow(10,-10))
                    RotationAxis1=RotationAxis1/RotationAxis1.l2norm();
                theta1= n.findAngle(nNew);
                rotQ1.setQuaternion(theta1,RotationAxis1);
                //find 2nd rotation axis. use vector normal to n and perpendicular to posVec. It lies in face plane.
                perpPVecNNew= perpPVecNNew.cross_product(newPosVec, nNew); //newPosVec is the new position of cNode[1]
                perpPVecN= perpPVecN.cross_product(cNode[1]->operator()(), n);
                perpPVecNQ.vector3DtoQuaternion(perpPVecN);
                perpPVecNQ= rotQ1*perpPVecNQ*rotQ1.conj();
                perpPVecN(0)= perpPVecNQ.get(1);
                perpPVecN(1)= perpPVecNQ.get(2);
                perpPVecN(2)= perpPVecNQ.get(3);
                RotationAxis2= RotationAxis2.cross_product(perpPVecN, perpPVecNNew);
                if (RotationAxis2.l2norm() >pow(10,-10))
                    RotationAxis2=RotationAxis2/RotationAxis2.l2norm();
                theta2= perpPVecN.findAngle(perpPVecNNew);
                rotQ2.setQuaternion(theta2,RotationAxis2);
                //finally find rotation of posVec and use it as 3rd rotation axis. when the face is being deformed and does not
                //keep its original shape, then its posVec can move independently so must find a 3rd rot axis.
                for (int k=0; k<2; k++)
                {
                    semiNewPosVecQ.vector3DtoQuaternion(cNode[k]->operator()());
                    semiNewPosVecQ= rotQ1*semiNewPosVecQ*rotQ1.conj();
                    semiNewPosVecQ= rotQ2*semiNewPosVecQ*rotQ2.conj();
                    semiNewPosVec(0)=semiNewPosVecQ.get(1);
                    semiNewPosVec(1)=semiNewPosVecQ.get(2);
                    semiNewPosVec(2)=semiNewPosVecQ.get(3);
                    
                    displaceVec= displace[cNode[k]->getIndex()-1];
                    newPosVec= cNode[k]->operator()() + displaceVec-originVec;
                    
                    RotationAxis3= RotationAxis3.cross_product(semiNewPosVec, newPosVec);
                    if (RotationAxis3.l2norm() >pow(10,-10))
                        RotationAxis3=RotationAxis3/RotationAxis3.l2norm();
                    theta3= semiNewPosVec.findAngle(newPosVec);
                    rotQ3.setQuaternion(theta3,RotationAxis3);
                    rotQ[cNode[k]->getIndex()-1]=rotQ3*rotQ2*rotQ1;
                }
            }
        }
    }
}

//Find translation of Boundary Nodes (used in the IDWquaternion method)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
void meshDisplacement<T,M,N>::findTranslations(node<T,N> **bNodes, int bNodesNum){
    quaternion<T> R;
    vector<T> r(3);
    node<T,N> *cNode;
    
    #pragma omp parallel for default(none) private(cNode,R) firstprivate(r) shared(bNodesNum,bNodes,translation,displace)
    {
        for (int i=0; i<bNodesNum; i++)
        {
            cNode= bNodes[i];//current node
            R.vector3DtoQuaternion(cNode->operator()());
            R= rotQ[cNode->getIndex()-1]*R*rotQ[cNode->getIndex()-1].conj();
            r.set(0,R.get(1));
            r.set(1,R.get(2));
            r.set(2,R.get(3));
            translation[cNode->getIndex()-1]= cNode->operator()() + displace[cNode->getIndex()-1] -r;
        }
    }
}

//find characteristic length
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
T meshDisplacement<T,M,N>::domainlength(node<T,N> **bNodes, int bNodesNum){
    T dLength=0., L;
    node<T,N> *cNode;
    vector<T> posVec(3);
    
   #pragma omp parallel for default(none) private(cNode, L, posVec) shared(bNodesNum, bNodes) reduction(max:dLength)
    {
        for (int i=0; i<bNodesNum; i++)
        {
            cNode= bNodes[i];
            posVec= cNode->operator()();
            L=posVec.l2norm();
            if (dLength < L)
                dLength= L;
        }
    }
    std::cout<<"\nInverse Distance Weighting: Characteristic Length= "<<dLength<<"\n";
    
    return dLength;
}

//write displaced grid to file
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
void meshDisplacement<T,M,N>::writeDisplacedMesh(createMesh<T,M,N> &m, const std::string filename, const std::string pathname){
    std::fstream outData;
    std::string filename1;
    int inval=0;
    T newVolume;
    cell<T,M,N> *cCell;
    const linkedList<cell<T,M,N> *> *cellNext;
    
    filename1= pathname + filename +".nod";
    writeNodFile(m,filename1);
    
    filename1= pathname + "list_invalid_cells.dat";
    outData.open(filename1,std::ios::out);
    if (!outData)
    {
        std::cout<<"\nError: Could not create file list_invalid_cells.dat.\n";
        exit(0);
    }
    outData.precision(17);
    
    outData<<"List of invalid cells:\n";
    outData<<"Cell ID "<<" New volume\n";
    
    cellNext= &m.cMesh.first();
    for (int i=0; i<m.cellsNumber; i++)
    {
        cCell= cellNext->getItem(0);
        newVolume= cCell->getVolume();
        if ( newVolume <= 0.)
        {
            outData<<cCell->getIndex()<<" "<<newVolume<<"\n";
            inval+=1;
        }
        cellNext= cellNext->readNext();
    }
    outData<<"\n"<<inval<<" invalid cells, out of "<<m.cellsNumber;
    std::cout<<"\n"<<inval<<" invalid cells, out of "<<m.cellsNumber<<"\n";
}

//write .nod file
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
void meshDisplacement<T,M,N>::writeNodFile(createMesh<T,M,N> &m, const std::string filename){
    std::fstream outData;
    node<T,N> **nodes;
    int nodesNum= inNodesNum + iBNodesNum + oBNodesNum;

    nodes= new node<T,N> *[nodesNum];
    #pragma omp parallel for default(none) shared(inNodesNum, nodes, inNodes)
    {
        for (int i=0; i<inNodesNum; i++)
            nodes[inNodes[i]->getIndex()-1]= inNodes[i];
    }
    #pragma omp parallel for default(none) shared(iBNodesNum, nodes, iBNodes)
    {
        for (int i=0; i<iBNodesNum; i++)
            nodes[iBNodes[i]->getIndex()-1]= iBNodes[i];
    }
    #pragma omp parallel for default(none) shared(oBNodesNum, nodes, oBNodes)
    {
        for (int i=0; i<oBNodesNum; i++)
            nodes[oBNodes[i]->getIndex()-1]= oBNodes[i];
    }
    
    outData.open(filename,std::ios::out);
    if (!outData)
    {
        std::cout<<"\nError: Could not create file"<<filename<<".nod.\n";
        exit(0);
    }
    outData.precision(17);
    
    outData<<m.getNodesNumber()<<"\n";
    for (int i=0; i<nodesNum; i++)
        outData<<nodes[i]->getLogFr()<<" ";
    outData<<"\n";
    for (int i=0; i<nodesNum; i++)
        outData<<nodes[i]->operator[](0)<<" ";
    outData<<"\n";
    for (int i=0; i<nodesNum; i++)
        outData<<nodes[i]->operator[](1)<<" ";
    outData<<"\n";
    for (int i=0; i<nodesNum; i++)
        outData<<nodes[i]->operator[](2)<<" ";
    outData.close();
    delete [] nodes;
}

//return displace
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
vector<T> *meshDisplacement<T,M,N>::getDisplace() const{
    return displace;
}

//find alfa
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
T meshDisplacement<T,M,N>::findAlfa(node<T,N> **bNodes, int bNodesNum, T cLength, T eta){
    T s[bNodesNum], sMean=0., maxS=0., alfa, temp=0.;
    
    #pragma omp parallel for default(none) shared(bNodesNum, s, displace, bNodes) reduction(+:sMean)
    {
        for (int i=0; i<bNodesNum; i++)
        {
            s[i]= displace[bNodes[i]->getIndex()-1].l2norm();
            sMean+=s[i];
        }
    }
    sMean= sMean/bNodesNum;
    #pragma omp parallel for default(none) private(temp) shared(bNodesNum, s, sMean) reduction(max:maxS)
    {
        for (int i=0; i<bNodesNum; i++)
        {
            temp= s[i] - sMean;
            if (temp > maxS)
                maxS= temp;
        }
    }
    return alfa=eta/cLength*maxS;
}

template<typename T, int M, int N>
T meshDisplacement<T,M,N>::findAlfa(node<T,N> **iBNodes, node<T,N> **oBNodes, int iBNodesNum, int oBNodesNum, T cLength, T eta){
    T sI[iBNodesNum], sO[oBNodesNum], sMean, sMeanI=0., sMeanO=0., maxS, maxSI=0., maxSO=0., alfa, temp;
    
    #pragma omp parallel for default(none) shared(iBNodesNum, sI, displace, iBNodes) reduction(+:sMeanI)
    {
        for (int i=0; i<iBNodesNum; i++)
        {
            sI[i]= displace[iBNodes[i]->getIndex()-1].l2norm();
            sMeanI+=sI[i];
        }
    }
    #pragma omp parallel for default(none) shared(oBNodesNum, sO, displace,oBNodes) reduction(+:sMeanO)
    {
        for (int i=0; i<oBNodesNum; i++)
        {
            sO[i]= displace[oBNodes[i]->getIndex()-1].l2norm();
            sMeanO+=sO[i];
        }
    }
    sMean= (sMeanI+sMeanO)/(iBNodesNum + oBNodesNum);
    
    #pragma omp parallel for default(none) private(temp) shared(iBNodesNum, sI, sMean) reduction(max:maxSI)
    {
        for (int i=0; i<iBNodesNum; i++)
        {
            temp= sI[i] - sMean;
            if (temp > maxSI)
                maxSI= temp;
        }
    }
    
    #pragma omp parallel for default(none) private(temp) shared(oBNodesNum, sO, sMean) reduction(max:maxSO)
    {
        for (int i=0; i<oBNodesNum; i++)
        {
            temp= sO[i] - sMean;
            if (temp > maxSO)
                maxSO= temp;
        }
    }
    if (maxSO > maxSI)
        maxS= maxSO;
    else
        maxS= maxSI;
    return alfa=eta/cLength*maxS;
}

//update mesh
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
void meshDisplacement<T,M,N>::updateMesh(createMesh<T,M,N> &m){
    node<T,N> *cNode;
    edge<T,N> *cEdge;
    face<T,N> *cFace;
    cell<T,M,N> *cCell;
    const linkedList<cell<T,M,N> *> *cellNext;
    int cEdgesNumber, cFacesNumber;
    
    //find new position of nodes
    #pragma omp parallel for default(none) private(cNode) shared(inNodes,displace)
    {
        for (int i=0; i<inNodesNum; i++)
        {
            cNode= inNodes[i];
            cNode->setLocationX(cNode->operator[](0) +displace[cNode->getIndex()-1](0));
            cNode->setLocationY(cNode->operator[](1) +displace[cNode->getIndex()-1](1));
            cNode->setLocationZ(cNode->operator[](2) +displace[cNode->getIndex()-1](2));
        }
    }
    #pragma omp parallel for default(none) private(cNode) shared(iBNodes,displace)
    {
        for (int i=0; i<iBNodesNum; i++)
        {
            cNode= iBNodes[i];
            cNode->setLocationX(cNode->operator[](0) +displace[cNode->getIndex()-1](0));
            cNode->setLocationY(cNode->operator[](1) +displace[cNode->getIndex()-1](1));
            cNode->setLocationZ(cNode->operator[](2) +displace[cNode->getIndex()-1](2));
        }
    }
    #pragma omp parallel for default(none) private(cNode) shared(oBNodes,displace)
    {
        for (int i=0; i<oBNodesNum; i++)
        {
            cNode= oBNodes[i];
            cNode->setLocationX(cNode->operator[](0) +displace[cNode->getIndex()-1](0));
            cNode->setLocationY(cNode->operator[](1) +displace[cNode->getIndex()-1](1));
            cNode->setLocationZ(cNode->operator[](2) +displace[cNode->getIndex()-1](2));
        }
    }


    //compute volume and area of the each cell, area and normal vector of each face and direction vector of each edge
    cellNext= &m.cMesh.first();
    for (int i=0; i<m.cellsNumber; i++)
    {
        cCell= cellNext->getItem(0);
        cCell->cellVolume();
        cCell->cellArea();
        cellNext= cellNext->readNext();
        
        cEdgesNumber= cCell->cellEdgesLength();
        #pragma omp parallel for num_threads(cEdgesNumber) default(none) private(cEdge) shared(cCell, cEdgesNumber)
        {
            for (int j=0; j<cEdgesNumber; j++)
            {
                cEdge= cCell->getCellEdge(j);
                if (cEdge->getFlag())
                    continue;
                
                cEdge->computeDirectionVector();
                cEdge->setFlagTrue();
            }
        }
        
        cFacesNumber= cCell->cellFacesLength();
        
        #pragma omp parallel for num_threads(cFacesNumber) default(none) private(cFace) shared(cCell, cFacesNumber)
        {
            for (int j=0; j<cFacesNumber; j++)
            {
                cFace= cCell->getCellFace(j);
                if (cFace->getFlag())
                    continue;
                cFace->findArea();
                cFace->findNormalVector();
                cFace->setFlagTrue();
            }
        }
    }
    
    //set flags to zero in case grid is updated again
    cellNext= &m.cMesh.first();
    for (int i=0; i<m.cellsNumber; i++)
    {
        cCell= cellNext->getItem(0);
        cellNext= cellNext->readNext();
        
        cEdgesNumber= cCell->cellEdgesLength();
        for (int j=0; j<cEdgesNumber; j++)
        {
            cEdge= cCell->getCellEdge(j);
            if (!cEdge->getFlag())
                continue;
            cEdge->setFlagFalse();
        }
        
        cFacesNumber= cCell->cellFacesLength();
        for (int j=0; j<cFacesNumber; j++)
        {
            cFace= cCell->getCellFace(j);
            if (!cFace->getFlag())
                continue;
            cFace->setFlagFalse();
        }
    }
}

//find cross point between face and line in 3D space
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
vector<T> meshDisplacement<T,M,N>::findCrossPoint(face<T,N> &f, vector<T> &vec)
{
    vector<T> n(3), cPoint(3);
    T nX, nY, nZ, xP, yP, zP, a, b, c, lamda;
    
    n=f.getNormalVector();
    vec= vec/vec.l2norm();
    
    nX= n(0);
    nY= n(1);
    nZ= n(2);
    xP= f.getFaceNode(0)->operator[](0);
    yP= f.getFaceNode(0)->operator[](1);
    zP= f.getFaceNode(0)->operator[](2);
    a= vec(0);
    b= vec(1);
    c= vec(2);
    
    lamda= (nX*xP + nY*yP + nZ*zP)/(nX*a + nY*b + nZ*c);
    
    cPoint(0)= a*lamda;
    cPoint(1)= b*lamda;
    cPoint(2)= c*lamda;
    
    return cPoint;
}

//find maximum coordinates of boundary
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, int M, int N>
void meshDisplacement<T,M,N>::findMaxCoords(node<T,N> **bNodes,int bNodesNum, T &xMin,T &xMax,T &yMin,T &yMax,T &zMin,T &zMax)
{
    node<T,N> *cNode;
    vector<T> cVec;
    
    xMin=  pow(10,10);
    xMax= -pow(10,10);
    yMin=  pow(10,10);
    yMax= -pow(10,10);
    zMin=  pow(10,10);
    zMax= -pow(10,10);
    
    #pragma omp parallel for default(none) private(cNode,cVec) shared(bNodes, bNodesNum) reduction(min: xMin, yMin, zMin) reduction(max: xMax, yMax, zMax)
    {
        for (int i=0; i<bNodesNum; i++)
        {
            cNode= bNodes[i];
            cVec= cNode->operator()();
            if (cVec(0) < xMin)
                xMin= cVec(0);
            if (cVec(0) > xMax)
                xMax= cVec(0);
            if (cVec(1) < yMin)
                yMin= cVec(1);
            if (cVec(1) > yMax)
                yMax= cVec(1);
            if (cVec(2) < zMin)
                zMin= cVec(2);
            if (cVec(2) > zMax)
                zMax= cVec(2);
        }
    }
}

template<typename T, int M, int N>
void meshDisplacement<T,M,N>::extractNodes(createMesh<T,M,N> &m){
    cell<T,M,N> *cCell;
    node<T,N> *cNode;
    const linkedList<cell<T,M,N> *> *cellNext;
    int countIn=0, countInB=0, countOutB=0;
    
    inNodes= new node<T,N> *[inNodesNum];
    iBNodes= new node<T,N> *[iBNodesNum];
    oBNodes= new node<T,N> *[oBNodesNum];
    
    cellNext= &m.cMesh.first();
    for (int i=0; i<m.cellsNumber; i++)
    {
        cCell= cellNext->getItem(0);
        for (int j=0; j<M; j++)
        {
            cNode= cCell->operator()(j);
            if (!cNode->getFlag())
            {
                if (cNode->getLogFr() == 0)
                {
                    inNodes[countIn]=cNode;
                    countIn++;
                    cNode->setFlagTrue();
                }
                else if (cNode->getLogFr() == 3)
                {
                    iBNodes[countInB]=cNode;
                    countInB++;
                    cNode->setFlagTrue();
                }
                else
                {
                    oBNodes[countOutB]=cNode;
                    countOutB++;
                    cNode->setFlagTrue();
                }
            }
        }
        cellNext= cellNext->readNext();
    }
    //set flagInList to false, so it can be used again
    #pragma omp parallel for default(none) shared(inNodesNum, inNodes)
    {
        for (int i=0; i<inNodesNum; i++)
            inNodes[i]->setFlagFalse();
    }
    #pragma omp parallel for default(none) shared(iBNodesNum, iBNodes)
    {
        for (int i=0; i<iBNodesNum; i++)
            iBNodes[i]->setFlagFalse();
    }
    #pragma omp parallel for default(none) shared(oBNodesNum, oBNodes)
    {
        for (int i=0; i<oBNodesNum; i++)
            oBNodes[i]->setFlagFalse();
    }
}

template<typename T, int M, int N>
void meshDisplacement<T,M,N>::extractFaces(createMesh<T,M,N> &m){
    cell<T,M,N> *cCell;
    face<T,N> *cFace;
    int cFacesNum;
    const linkedList<cell<T,M,N> *> *cellNext;
    int countInB=0, countOutB=0;
    
    iBFaces= new face<T,N> *[iBFacesNum];
    oBFaces= new face<T,N> *[oBFacesNum];
    
    cellNext= &m.cMesh.first();
    for (int i=0; i<m.cellsNumber; i++)
    {
        cCell= cellNext->getItem(0);
        cFacesNum= cCell->cellFacesLength();
        for (int j=0; j<cFacesNum; j++)
        {
            cFace= cCell->getCellFace(j);
            if (cFace->getFlag())
                continue;
            
            if (cFace->getFaceLogFr() == 3)
            {
                iBFaces[countInB]=cFace;
                countInB++;
                cFace->setFlagTrue();
            }
            else if (cFace->getFaceLogFr() == 4)
            {
                oBFaces[countOutB]=cFace;
                countOutB++;
                cFace->setFlagTrue();
            }
        }
        cellNext= cellNext->readNext();
    }
    //set flagInList to false, so it can be used again
    #pragma omp parallel for default(none) shared(iBFacesNum, iBFaces)
    {
        for (int i=0; i<iBFacesNum; i++)
            iBFaces[i]->setFlagFalse();
    }
    #pragma omp parallel for default(none) shared(oBFacesNum, oBFaces)
    {
        for (int i=0; i<oBFacesNum; i++)
            oBFaces[i]->setFlagFalse();
    }
}


template<typename T,int M, int N>
void meshDisplacement<T,M,N>::calculateDisplacement(node<T,N> **nodes, int nodesNum){
    double xDisp, yDisp, zDisp, x0, y0, z0, rotationAngle, xComp, yComp, zComp;
    vector<double> rotationAxis(3), displ(3), dr(3), drNew(3), rcNew(3);
    node<T,N> *cNode;
    quaternion<double> inQ, drNewQ; //inQ: quaternion to store the initial rotation that will be applied to the object
                                   //drNewQ: quaternion to store the transformed vector drNew
    
    std::cout<<"Enter Displacement in x-axis: ";
    std::cin>>xDisp;
    std::cout<<"\nEnter Displacement in y-axis: ";
    std::cin>>yDisp;
    std::cout<<"\nEnter Displacement in z-axis: ";
    std::cin>>zDisp;
    std::cout<<"\nEnter new x-center: ";
    std::cin>>x0;
    std::cout<<"\nEnter new y-center: ";
    std::cin>>y0;
    std::cout<<"\nEnter new z-center: ";
    std::cin>>z0;
    std::cout<<"\nEnter axis of rotation:\n";
    std::cout<<"x-component: ";
    std::cin>>xComp;
    std::cout<<"\ny-component: ";
    std::cin>>yComp;
    std::cout<<"\nz-component: ";
    std::cin>>zComp;
    std::cout<<"\nEnter rotation angle (in degrees): ";
    std::cin>>rotationAngle;
    
    displ.set(0,xDisp);
    displ.set(1,yDisp);
    displ.set(2,zDisp);
    rcNew.set(0,x0);
    rcNew.set(1,y0);
    rcNew.set(2,z0);
    rotationAxis.set(0,xComp);
    rotationAxis.set(1,yComp);
    rotationAxis.set(2,zComp);
    rotationAxis= rotationAxis/rotationAxis.l2norm();
    rotationAngle=rotationAngle*M_PI/180.;
    
    inQ.setQuaternion(rotationAngle,rotationAxis);
    
    #pragma omp parallel for default(none) private(cNode, dr, drNew, drNewQ) shared(nodesNum, nodes, displace, displ, rcNew, inQ)
    {
        for (int i=0; i<nodesNum; i++)
        {
            cNode= nodes[i];
            dr= cNode->operator()() - rcNew;
            drNew= dr;
            drNewQ.vector3DtoQuaternion(drNew);
            drNewQ= inQ*drNewQ*inQ.conj();
            drNew.set(0,drNewQ.get(1));
            drNew.set(1,drNewQ.get(2));
            drNew.set(2,drNewQ.get(3));
            
            displace[cNode->getIndex()-1]= rcNew + drNew - cNode->operator()() + displ;
        }
    }
}

template<typename T,int M,int N>
void meshDisplacement<T,M,N>::findBoundaryTurningPoints(face<T,N> **bFaces, int bFacesNum){
    face<T,N> *cFace=NULL, *neiFace=NULL;
    node<T,N> *cNode=NULL;
    cell<T,M,N> *cCell=NULL;
    int flag1=0, flag2=0, fNodeInd[3], cNodeInd, count;
    T theta;
    int FNodesLength1, FNodesLength2, shNodesLength, cFacesLength;
    
   #pragma omp parallel for default(none) private(cFace, FNodesLength1, flag1, cNode, shNodesLength, cCell, cFacesLength, flag2, neiFace, FNodesLength2, theta) shared(bFacesNum, bFaces, oBTNodes, oBTNodesNum)
    {
        for (int i=0; i<bFacesNum; i++)
        {//forLoop1
            cFace= bFaces[i];
            FNodesLength1= cFace->getFaceNodesLength();
            for (int j=0; j<FNodesLength1; j++)
            {//forLoop2
                flag1=0;
                cNode= cFace->getFaceNode(j);
                oBTNodesNum= oBTNodes.length();
                for (int m1=1; m1<oBTNodesNum; m1++)
                    if (*cNode == *oBTNodes(m1))
                    {
                        flag1=1;
                        break;
                    }
                if (flag1)
                    continue;
                
                shNodesLength= cNode->getSharingCells();
                for (int k=0; k<shNodesLength; k++)
                {//forLoop3
                    cCell= (tetrahedron *)cNode->getSharingCellAddress(k);
                    //find the boundary face of this cell
                    cFacesLength= cCell->cellFacesLength();
                    for (int m1=0; m1<cFacesLength; m1++)
                    {
                        flag1=0;
                        flag2=0;
                        neiFace= cCell->getCellFace(m1);
                        FNodesLength2= neiFace->getFaceNodesLength();
                        for (int m2=0; m2<FNodesLength2; m2++)
                            if (neiFace->getFaceNode(m2)->getLogFr() != 4)
                            {
                                flag1=1;
                                break;
                            }
                        if (flag1 && m1 == cFacesLength-1)
                        {
                            flag2=1;
                            break;
                        }
                        else if (flag1)
                            continue;
                        else
                            break;
                    }
                    if (flag2)
                        continue;
                    theta= acos(cFace->getNormalVector()*neiFace->getNormalVector())/
                    (cFace->getNormalVector().l2norm()*neiFace->getNormalVector().l2norm());
                    /*subjective judgement of the programmer: if the change in orientation of the two boundary faces is
                     larger than pi/8. (=22.5 degrees), then this is a big change and the node in between the two segments must
                     be kept in position in order to keep the geometry of the boundary unchanged (this node is gonna be needed
                     when the sliding boundary technique is used.*/
                    if (fabs(theta) > M_PI/8.)
                    {
                        #pragma omp critical
                        {
                            if (oBTNodesNum==1)
                            {
                                oBTNodes.insertNextItem(cNode);
                                oBTNodesNum++;
                            }
                            else
                            {
                                flag2=1;
                                for (int m1=1; m1<oBTNodesNum; m1++)
                                {
                                    if (*cNode == *oBTNodes(m1))
                                    {
                                        flag2=0;
                                        break;
                                    }
                                }
                                if (flag2)
                                {
                                    oBTNodes.insertNextItem(cNode);
                                    oBTNodesNum++;
                                }
                            }
                        }
                    }
                }//end forLoop3
            }//end forLoop2
        }//end forLoop1
    }

    if (oBTNodesNum > 1)
    {
        oBTNodes.dropFirstItem();
        oBTNodesNum--;
        
        oBTFaces.setListSize(oBTNodesNum);
        #pragma omp parallel for default(none) private(cNode, cNodeInd, cFace, fNodeInd, count) shared(oBTNodesNum, oBTNodes, bFaces, oBTFaces)
        {//start of parallel region
            for (int i=0; i<oBTNodesNum; i++)
            {
                cNode= oBTNodes(i);
                cNodeInd= cNode->getIndex();
                count=0;
                for (int j=0; j<oBFacesNum; j++)
                {
                    cFace= bFaces[j];
                    for (int k=0; k<3; k++)
                        fNodeInd[k]= cFace->getFaceNode(k)->getIndex();
                    if ( (fNodeInd[0] == cNodeInd) || (fNodeInd[1] == cNodeInd) || (fNodeInd[2] == cNodeInd) )
                    {
                            if (count == 0)
                                oBTFaces.put(i,cFace);
                            else
                                oBTFaces(i).insertNextItem(cFace);
                            count++;
                    }
                }
            }
        }//end of parallel region
    }
}
#endif /* meshDisplacement_h */
