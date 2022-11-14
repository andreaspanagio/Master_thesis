adaptDisplace: main.o createMesh.o meshDisplacement.o MSH3Dtranslator.o vtkUnstrMeshTranslator.o list.o linkedList.o quaternion.o vector.o node.o edge.o face.o cell.o mesh.o
	g++-11 -fopenmp main.o createMesh.o meshDisplacement.o MSH3Dtranslator.o vtkUnstrMeshTranslator.o list.o linkedList.o quaternion.o vector.o node.o edge.o face.o cell.o mesh.o -o adaptDisplace

main.o: main.cpp
	g++-11 main.cpp
	
createMesh.o: createMesh.cpp
	g++-11 -fopenmp createMesh.h
	
meshDisplacement.o: meshDisplacement.h
	g++-11 -fopenmp meshDisplacement.h
	
MSH3Dtranslator.o: MSH3Dtranslator.h
	g++-11 MSH3Dtranslator.h
	
vtkUnstrMeshTranslator.o: vtkUnstrMeshTranslator.h
	g++-11 -Xpreprocessor vtkUnstrMeshTranslator.h
	
list.o: list.h
	g++-11 -Xpreprocessor list.h
	
linkedList.o: linkedList.h
	g++-11 -Xpreprocessor linkedList.h
	
quaternion.o: quaternion.h
	g++-11 -Xpreprocessor quaternion.h
	
vector.o: vector.h
	g++-11 -Xpreprocessor vector.h
	
node.o: node.h
	g++-11 -Xpreprocessor node.h
	
edge.o: edge.h
	g++-11 -Xpreprocessor edge.h
	
face.o: face.h
	g++-11 -Xpreprocessor face.h
	
cell.o: cell.h
	g++-11 -Xpreprocessor cell.h
	
mesh.o: mesh.h
	g++-11 -Xpreprocessor mesh.h
	
#
# Clean out object files and the executable.
#
clean:
	rm *.o adaptDisplace
