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



#ifndef ZSPACE_MODULES_PLANARISE_H
#define ZSPACE_MODULES_PLANARISE_H



#pragma once
#include <vector>
#include <string>
using namespace std;

#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>
#include <headers/zInterface/functionsets/zFnParticle.h>

#include <headers/zSpace_Modules.h>

namespace  zSpace
{

	//--------------------------
	//---- GLOBAL VARIABLES
	//--------------------------

	/*!	\brief container of  particle objects  */
	vector<zObjParticle> o_Particles;

	/*!	\brief container of particle function set  */
	vector<zFnParticle> fnParticles;

	/*!	\brief container of vertex positions  */
	zPointArray vertexPositions;

	/*!	\brief 2D container of polygon vertex connectivity  */
	zInt2DArray polygons;

	/*!	\brief 2D container of triangle vertex connectivity */
	zInt2DArray triangles;

	//--------------------------
	//---- UTILITY METHOD
	//--------------------------

	/*! \brief This utility method computes the polyhedral volumes and centers of each face of the mesh.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polygons				- input 2D container of polygons.
	*	\param	[in]	_triangles				- input 2D container of triangles. 
	* 	\param	[in]	_core					- input object of zUtilsCore.
	*  	\param	[out]	fCenters				- output container of face center positions. 
	*  	\param	[out]	fVolumes				- output container of face volumes.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void computeFaceVolumes(zPointArray& _vertexPositions, zInt2DArray& _polygons, zInt2DArray& _triangles, zUtilsCore& _core, zPointArray& fCenters, zFloatArray& fVolumes);

	/*! \brief This utility method computes the face normals of the mesh.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polygons				- input 2D container of polygons.
	*	\param	[in]	_triangles				- input 2D container of triangles.
	* 	\param	[in]	_core					- input object of zUtilsCore.
	*  	\param	[out]	fNormals				- output container of face normals.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void computeFaceNormals(zPointArray& _vertexPositions, zInt2DArray& _polygons, zInt2DArray& _triangles, zUtilsCore& _core, zVectorArray& fNormals);

	//--------------------------
	//---- EXTERNAL METHOD
	//--------------------------

	/*! \brief This method planarises the input quad mesh using the closest distance between diagonals as the tolerance test.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	* 	\param	[in]	updatetopology			- input boolean to update topology.
	*  	\param	[in]	numIterations			- input number of iterations to run the solver.
	*  	\param	[in]	tolerance				- input planarity tolerance.
	*  	\param	[out]	outVertexPositions		- output container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*  	\param	[out]	outDeviations			- output container of planarity deviations per face/polygon.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void quadPlanarise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, bool updatetopology, int numIterations, double tolerance, double* outVertexPositions, double* outDeviations);

	/*! \brief This method planarises the input mesh using the polyhedral volume of the face/polygon as the tolerance test.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	_triCounts				- input container of number of triangles per polygon of the mesh.
	*	\param	[in]	_triConnects			- input container of triangle connectivity. Collapsed 1D array of size numFaces * (numtriangles per face * 3).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	* 	\param	[in]	updatetopology			- input boolean to update topology.
	*  	\param	[in]	numIterations			- input number of iterations to run the solver.
	*  	\param	[in]	tolerance				- input planarity tolerance.
	*  	\param	[out]	outVertexPositions		- output container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*  	\param	[out]	outDeviations			- output container of planarity deviations per face/polygon.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void planarise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int* _fixedVertices, int numVerts, int numFaces, int numFixed, bool updatetopology, int numIterations, double tolerance, double* outVertexPositions, double* outDeviations);

	/*! \brief This method returns the input mesh planarity deviation using the polyhedral volume of the face/polygon as the tolerance test.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	_triCounts				- input container of number of triangles per polygon of the mesh.
	*	\param	[in]	_triConnects			- input container of triangle connectivity. Collapsed 1D array of size numFaces * (numtriangles per face * 3).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	* 	\param	[in]	updatetopology			- input boolean to update topology.
	*  	\param	[out]	outDeviations			- output container of planarity deviations per face/polygon.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void planarityDeviation(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, bool updatetopology, double* outDeviations);

	/*! \brief This method returns the input quad mesh planarity deviation using the closest distance between diagonals as the tolerance test.
*
*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
* 	\param	[in]	numVerts				- input number of vertices in the mesh.
*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
* 	\param	[in]	updatetopology			- input boolean to update topology.
*  	\param	[out]	outDeviations			- output container of planarity deviations per face/polygon.
*	\since version 0.0.4
*/
	extern "C" ZSPACE_MODULES void quadPlanarityDeviation(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, bool updatetopology, double* outDeviations);

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zSpace_Planarise.cpp>
#endif

#endif
