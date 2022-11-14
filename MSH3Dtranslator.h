//
//  MSH3Dtranslator.h
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 24/7/22.
//

/*
 This is a class that takes as an input an msh file and returns a .nod and . ele file.
 The msh file should be a mesh that has been created around a surface imported in GMSH as a
 stl file.
 */
#ifndef MSH3Dtranslator_h
#define MSH3Dtranslator_h


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

class MSH3Dtranslator{
public:
    MSH3Dtranslator(const std::string filename, const std::string pathname=""); //constructor
    MSH3Dtranslator(const std::string filename, const std::string filenameNew, const std::string pathname); //constructor
    
    double getXc();
    double getYc();
    double getZc();
private:
    double Xc; // x-coordinate centre of inner boudary surface.
    double Yc; // y-coordinate centre of inner boudary surface.
    double Zc; // z-coordinate centre of inner boudary surface
};

//write .nod and .ele file from .msh file
MSH3Dtranslator::MSH3Dtranslator(const std::string filename, const std::string pathname)
{
    std::ifstream inData;
    std::fstream outData;
    std::string filename1= pathname + filename + ".msh";
    std::string filename2= pathname + filename + ".nod";
    std::string filename3= pathname + filename + ".ele";
    std::string keyword;
    int fileIntElement, ElementType, numberofTetrahedrons, numberofElements, kk1, kk2;
    int nodesNumber, nodesBlocksNumber, blNodesNum, cLogFr;
    double fileDoubleElement, sumX, sumY, sumZ;
    int *logFr;
    double *coords[3];
    
    //open filename1.msh to read its data
    inData.open(filename1,std::ios::in);
    if (!inData)
    {
        std::cout<<"Error: No such file\n";
        std::cout<<"Is this filepath and filename correct?\n"<<filename1<<"\n";
        exit(0);
    }
    inData.seekg(0,std::ios::beg); //set runner to read form the start of the open file

    //search for $Nodes section
    while(std::getline(inData, keyword))
    {
        if (keyword == "$Nodes")
            break;
    }
    inData>>nodesBlocksNumber; //number of nodes Blocks. gmsh partitions the mesh in these blocks for some reason.
    inData>>nodesNumber; // number of nodes in the mesh
    inData>>fileIntElement; // skip
    inData>>fileIntElement; // skip

    //initialize arrays
    logFr= new int[nodesNumber];
    for (int i=0;i<3;i++)
        coords[i]= new double[nodesNumber]{0.};
    
    sumX=0.;
    sumY=0.;
    sumZ=0.;
    Xc=0.;
    Yc=0.;
    Zc=0.;
    kk1=0;
    kk2=0;
    for (int k=0; k<nodesBlocksNumber; k++)
    {
        for (int i=0; i<3; i++)
            inData>>fileIntElement; //move the cursor 3 elements. These 3 elements contain no useful info
        inData>>blNodesNum; //number of nodes in block k
        //write the logFr for inner boundary nodes
        if (k==0)
            cLogFr= 3;
        else if (k==1)
            cLogFr= 4;
        else
            cLogFr= 0;
        for (int i=0; i<blNodesNum; i++)
        {
            logFr[kk1]= cLogFr; //currentLogFr;
            kk1+=1;
        }
        //skip nodes index
        for (int i=0; i<blNodesNum; i++)
            inData>>fileIntElement;
        //save coords of each node
        if (k == 0)
        {
            for (int i=0; i<blNodesNum; i++)
            {
                for (int j=0; j<3; j++)
                {
                    inData>>fileDoubleElement;
                    coords[j][kk2]= fileDoubleElement;
                }
                sumX+= coords[0][kk2];
                sumY+= coords[1][kk2];
                sumZ+= coords[2][kk2];
                kk2+=1;
            }
            Xc= sumX/blNodesNum;
            Yc= sumY/blNodesNum;
            Zc= sumZ/blNodesNum;
        }
        else
        {
            for (int i=0; i<blNodesNum; i++)
            {
                for (int j=0; j<3; j++)
                {
                    inData>>fileDoubleElement;
                    coords[j][kk2]= fileDoubleElement;
                }
                kk2+=1;
            }
        }
    }
    //write the .ele file
    //open filename3 to write the necessary data from filename1
    outData.open(filename3,std::ios::out);
    if (!outData)
    {
        std::cout<<"\nError: Could not create file"<<filename<<".nod.\n";
        exit(0);
    }
    outData.precision(17);
    //search for $Elements section
    while(std::getline(inData, keyword))
    {
        if (keyword == "$Elements")
            break;
    }
    for (int i=0; i<4; i++)
        inData>>fileIntElement; //skip
    
    while (1)
    {
        inData>>fileIntElement; //skip this entry
        inData>>fileIntElement; //skip this entry
        inData>>ElementType;
        if (ElementType == 4)
            break;
        inData>>numberofElements;
        if (ElementType == 15)
        {
            for (int i=0; i<numberofElements; i++)
            {
                inData>>fileIntElement; //skip this entry
                inData>>fileIntElement; //skip this entry
            }
        }
        else if (ElementType == 1)
        {
            for (int i=0; i<numberofElements; i++)
            {
                inData>>fileIntElement; //skip this entry
                inData>>fileIntElement; //skip this entry
                inData>>fileIntElement; //skip this entry
            }
        }
        else if (ElementType == 2)
        {
            for (int i=0; i<numberofElements; i++)
            {
                inData>>fileIntElement; //skip this entry
                inData>>fileIntElement; //skip this entry
                inData>>fileIntElement; //skip this entry
                inData>>fileIntElement; //skip this entry
            }
        }
    }
    inData>>numberofTetrahedrons;
    outData<<numberofTetrahedrons;// store number of tetrahedrons
    outData<<"\n";
    for (int i=0; i<numberofTetrahedrons; i++)
    {
        inData>>fileIntElement; //skip
        for (int j=0; j<4; j++)
        {
            inData>>fileIntElement;
            outData<<fileIntElement<<" ";
        }
        outData<<"\n";
    }
    outData.close();
    //the .ele file is complete
    inData.close();
    //done reading .msh file
    
    //write the .nod file
    //open filename2 to write the necessary data
    outData.open(filename2,std::ios::out);
    if (!outData)
    {
        std::cout<<"\nError: Could not create file"<<filename2<<".nod.\n";
        exit(0);
    }
    outData.precision(17);
    
    outData<<nodesNumber; //write number of Nodes to .nod file
    outData<<"\n"; //change line
    for (int i=0; i<nodesNumber; i++)
        outData<<logFr[i]<<" ";
    outData<<"\n";
    
    for (int i=0; i<nodesNumber; i++)
        outData<<coords[0][i]-Xc<<" ";
        outData<<"\n";
    for (int i=0; i<nodesNumber; i++)
        outData<<coords[1][i]-Yc<<" ";
        outData<<"\n";
    for (int i=0; i<nodesNumber; i++)
        outData<<coords[2][i]-Zc<<" ";
    outData.close();
    //the .nod file is complete
    delete [] logFr;
    for (int i=0; i<3; i++)
        delete [] coords[i];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Write .msh file from .nod file
MSH3Dtranslator::MSH3Dtranslator(const std::string filename, const std::string filenameNew, const std::string pathname)
{
    std::ifstream inData;
    std::fstream outData;
    std::string filename1= pathname + filename + ".msh";
    std::string filename2= pathname + filenameNew + ".msh";
    std::string filename3= pathname + filenameNew + ".nod";
    std::string keyword;
    int fileIntElement, blNodesNumber, nodesNumber, kk;
    double fileDoubleElement;
    double *coords[3];
    
    //open filename3.nod to read its data
    inData.open(filename3,std::ios::in);
    if (!inData)
    {
        std::cout<<"Error: No such file\n";
        std::cout<<"Is this filepath and filename correct?\n"<<filename3<<"\n";
        exit(0);
    }
    inData.seekg(0,std::ios::beg); //set runner to read form the start of the open file
    
    inData>>nodesNumber;
    for (int i=0; i<3; i++)
        coords[i]= new double[nodesNumber]{0.};
    //set file pointer to point after the logFr values at the x-coordinate of the first node
    for (int i=0; i<nodesNumber; i++)
        inData>>fileIntElement;
    for (int i=0; i<3; i++)
        for (int j=0; j<nodesNumber; j++)
            inData>>coords[i][j];
    inData.close();
    
    //open filename1.msh to copy its data to filename2.msh
    inData.open(filename1,std::ios::in);
    if (!inData)
    {
        std::cout<<"Error: No such file\n";
        std::cout<<"Is this filepath and filename correct?\n"<<filename1<<"\n";
        exit(0);
    }
    inData.seekg(0,std::ios::beg); //set runner to read form the start of the open file
    
    //write elements in the new .msh file
    //open filename2.msh to write data from filename1.msh and from coords array
    outData.open( filename2,std::ios::out);
    if (!outData)
    {
        std::cout<<"\nError: Could not create file"<<filename2<<".\n";
        exit(0);
    }
    outData.precision(17);
    
    //copy all the entries until you find $Entities. Entry $Entities is not copied.
    while (std::getline(inData, keyword))
    {
        if (keyword == "$Entities")
            break;
        outData<<keyword<<"\n";
    }
    //skip all the entries until you find the entry $Nodes. when you find it, copy it to the file and break.
    while (std::getline(inData, keyword))
    {
        if (keyword == "$Nodes")
        {
            outData<<keyword<<"\n";
            break;
        }
    }
    //write the new coords to filename3.msh
    inData>>fileIntElement;
    outData<<fileIntElement<<" ";
    inData>>nodesNumber;
    outData<<nodesNumber<<" ";
    inData>>fileDoubleElement;
    outData<<fileDoubleElement<<" ";
    inData>>fileDoubleElement;
    outData<<fileDoubleElement<<" ";
    outData<<"\n";
    kk=0;
    //copy blocks of nodes (number of blocks= 3)
    for (int k=0; k<3; k++)
    {
        for (int i=0; i<3; i++)
        {
            inData>>fileIntElement;
            outData<<fileIntElement<<" ";
        }
        inData>>blNodesNumber; //number of nodes in Block k
        outData<<blNodesNumber;
        outData<<"\n";
        for (int i=0; i<blNodesNumber; i++)
        {
            inData>>fileIntElement; //index of Nodes;
            outData<<fileIntElement;
            outData<<"\n";
        }
        for (int i=0; i<blNodesNumber; i++)
        {
            inData>>fileDoubleElement; //skip this value
            outData<<coords[0][kk]<<" ";
            inData>>fileDoubleElement; //skip this value
            outData<<coords[1][kk]<<" ";
            inData>>fileDoubleElement; //skip this value
            outData<<coords[2][kk]<<" ";
            outData<<"\n";
            kk+=1;
        }
    }
    //copy the rest of the file as it is.
    while (std::getline(inData, keyword))
    {
        outData<<keyword<<"\n";
        if (keyword == "$EndElements")
            break;
    }
    inData.close();
    outData.close();
    
    for (int i=0; i<3; i++)
        delete [] coords[i];
}

double MSH3Dtranslator::getXc(){
    return Xc;
}

double MSH3Dtranslator::getYc(){
    return Yc;
}

double MSH3Dtranslator::getZc(){
    return Zc;
}
#endif /* MSH3Dtranslator_h */
