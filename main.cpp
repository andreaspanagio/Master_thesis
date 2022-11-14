//
//  main.cpp
//  Master_Thesis_C++
//
//  Created by Andreas Panagiotopoulos on 15/11/21.
//

#include <iostream>
#include <filesystem>
#include <string>

#include "createMesh.h"
#include "meshDisplacement.h"
#include "MSH3Dtranslator.h"
#include "vtkUnstrMeshTranslator.h"

int main()
{
    std::cout.precision(10);
    std::string pathname="";
    std::string filename="";
    std::string meshRes,filename1,filename2;
    int method, flag=0;
    bool slNodes;
    
    if (filename.find(".msh") != std::string::npos)
    {
        filename.resize(filename.size() - 4);
        MSH3Dtranslator(filename,pathname);
        flag=1;
    }
    std::cout<<"Reading grid...\n";
    createMesh<double,4,3> mesh1(filename,pathname); //tetrahedron 3D mesh
    
    mesh1().gridQuality(); //calculate grid quality
    vtkUnstrMeshTranslator(mesh1,filename,pathname);
    meshDisplacement<double,4,3> displ(mesh1,pathname); //ask user to input displacements
    std::cout<<"Select Solution Method/Model:\n";
    std::cout<<"-----------------------------\n";
    std::cout<<"<1> IDW-Plain (Inverse Distance Weighting)\n";
    std::cout<<"<2> IDW with Rotations-Translations (Quaternions)\n";
    std::cin>>method;
    std::cout<<"Do you want to use the sliding nodes technique? (enter 0 for NO or 1 for YES)\n";
    std::cin>>slNodes;
    
    switch (method) {
        case 1:
            displ.IDW(slNodes);
            break;
        case 2:
            displ.IDWquaternion(slNodes);
            break;
        default:
            std::cout<<"Please insert 1 (for IDW-Plain) or 2 (for IDW with Rotations-Translations).\n";
            exit(0);
    }
    displ.updateMesh(mesh1); //update the grid after the displacement
    mesh1().gridQuality(); //calculate new grid quality
    
    
    std::cout<<"\nEnter the RootFileName of the resulting grid: ";
    std::cin>>meshRes;
    displ.writeDisplacedMesh(mesh1,meshRes,pathname); //write the resulted grid into a file
    
    filename1= pathname + filename + ".ele";
    filename2= pathname + meshRes +".ele";
    std::filesystem::copy_file(filename1, filename2);
    //translate .nod and .ele files to .msh file
    if (flag == 1)
        MSH3Dtranslator(filename,meshRes,pathname);
    vtkUnstrMeshTranslator(mesh1,meshRes,pathname);
    std::cout<<"All files have been writen\n";
    return 0;
}
