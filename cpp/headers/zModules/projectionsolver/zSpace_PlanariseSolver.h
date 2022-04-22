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

#ifndef ZSPACE_MODULES_PLANARISESOLVER_H
#define ZSPACE_MODULES_PLANARISESOLVER_H

#pragma once
#include <vector>
#include <string>
using namespace std;

#include <headers/zInterface/functionsets/zFnParticle.h>

#include <headers/zModules/base/zSpace_Modules.h>
#include <headers/zModules/base/zSpace_ComputeMesh.h>
#include <headers/zModules/base/zSpace_MeshUtilities.h>

#include <headers/zModules/projectionsolver/zSpace_ProjectionForces.h>

namespace  zSpace
{
	//--------------------------
	//---- GLOBAL VARIABLES
	//--------------------------
	zComputeMesh planariseMesh;
	zPlanarType planarisationType;

	//--------------------------
	//---- EXTERNAL METHODS FOR PLANARISATION
	//--------------------------

	/*! \brief This method initialises the solver for planarisation.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	_triCounts				- input container of number of triangles per polygon of the mesh.
	*	\param	[in]	_triConnects			- input container of triangle connectivity. Collapsed 1D array of size numFaces * (numtriangles per face * 3).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	*  	\param	[in]	volPlanarise			- input boolean to set planarity solver type - quad diagonal distance or volume.
	*  	\param	[out]	outDeviations			- output container of planarity deviations per face/polygon.
	*	\return			int					- output boolean - true if setup is successful.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES int planariseSolver_initialise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, bool volPlanarise, double* outDeviations);

	/*! \brief This method planarises the mesh.
	*
	*  	\param	[in]	numIterations			- input number of iterations to run the solver.
	*  	\param	[in]	tolerance				- input planarity tolerance.
	*  	\param	[in]	volPlanarise			- input boolean to set planarity solver type - quad diagonal distance or volume.
	*  	\param	[out]	outVertexPositions		- output container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*  	\param	[out]	outDeviations			- output container of planarity deviations per face/polygon.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void planariseSolver_compute(int numIterations, double tolerance, double* outVertexPositions, double* outDeviations);

	//--------------------------
	//---- EXTERNAL METHODS FOR CONSTRAINTS
	//--------------------------

	/*! \brief This method makes the vertices specified by the input contatiner fixed.
	*
	*	\param	[in]	_fixedVertices			- input container of anchor point indicies.
	*	\param	[in]	numFixed				- number of fixed vertices.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void planariseSolver_setFixed(int* _fixedVertices, int numFixed);

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/projectionsolver/zSpace_PlanariseSolver.cpp>
#endif

#endif
