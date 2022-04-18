// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2022 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>
//



#ifndef ZSPACE_MODULES_FACTORY_H
#define ZSPACE_MODULES_FACTORY_H



#pragma once
#include <vector>
#include <string>
using namespace std;

#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>
#include <headers/zInterface/functionsets/zFnParticle.h>

#include <headers/zModules/base/zSpace_Modules.h>
#include <headers/zModules/base/zSpace_ModulesEnumerators.h>
#include <headers/zModules/base/zSpace_ComputeMesh.h>

namespace  zSpace
{

	//--------------------------
	//---- GLOBAL VARIABLES
	//--------------------------

	//--------------------------
	//----  CREATE METHODS
	//--------------------------

	/*! \brief This method creates a zObjMesh from connectivity and position containers.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	*  	\param	[out]	outMesh					- output mesh object.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void createMeshOBJ(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, zObjMesh& outMesh);

	/*! \brief This method creates a zComputeMesh from connectivity and position containers.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	*  	\param	[in]	makeDynamic				- input boolean initialise the particle of the mesh.
	*  	\param	[out]	outMesh					- output compute mesh object.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void createComputeMesh(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, bool makeDynamic, zComputeMesh& outMesh);

	//--------------------------
	//----  SET METHODS
	//--------------------------

	/*! \brief This method sets the triangles of input zComputeMesh from input containers.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	* 	\param	[in]	_triCounts				- input container of number of triangles per polygon of the mesh.
	*	\param	[in]	_triConnects			- input container of triangle connectivity. Collapsed 1D array of size numFaces * (numtriangles per face * 3).
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void setTriangles(zComputeMesh& inMesh, int numFaces, int* _triCounts, int* _triConnects);

	//--------------------------
	//----  UPDATE METHODS
	//--------------------------

	/*! \brief This method updates the V matrix of input zComputeMesh from the current vertex positions.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void updateMatrixV(zComputeMesh& inMesh);

	/*! \brief This method updates the FTris matrix of input zComputeMesh from the current triangle connectivity.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void updateMatrixFTris(zComputeMesh& inMesh);

	//--------------------------
	//----  COMPUTE METHODS
	//--------------------------

		/*! \brief This utility method computes the polyhedral volumes and centers of each face of the mesh.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	*  	\param	[out]	fCenters				- output container of face center positions.
	*  	\param	[out]	fVolumes				- output container of face volumes.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void computeFaceVolumes(zComputeMesh& inMesh, zPointArray& fCenters, zDoubleArray& fVolumes);

	/*! \brief This utility method computes the face normals of the mesh.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	*  	\param	[out]	fNormals				- output container of face normals.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void computeFaceNormals(zComputeMesh& inMesh, zVectorArray& fNormals);

	/*! \brief This utility method computes the face normals of the mesh.
	*
	*	\param	[in]	inMesh					- input compute mesh object.	
	*  	\param	[out]	fDeviations				- output container of face planarity deviations.
	*	\warning only works on quad mesh
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void computeQuadPlanarityDeviation(zComputeMesh& inMesh, zDoubleArray& fDeviations);

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/base/zSpace_MeshUtilities.cpp>
#endif

#endif
