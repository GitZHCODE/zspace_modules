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

#ifndef ZSPACE_MODULES_CCFSOLVER_H
#define ZSPACE_MODULES_CCFSOLVER_H

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
	zComputeMesh ccfMesh;
	zPlanarSolverType ccf_planarisationType;

	//--------------------------
	//---- EXTERNAL METHOD FOR CCF
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
	*  	\param	[out]	outPlanarityDeviations	- output container of planarity deviations per face/polygon.
	*  	\param	[out]	outGaussianCurvatures	- output container of gaussian curvatures per vertex.
	*	\return			bool					- output boolean - true if setup is successful.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES bool ccfSolver_initialise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, double* outPlanarityDeviations, double* outGaussianCurvatures);

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
	extern "C" ZSPACE_MODULES void ccfSolver_compute(int numIterations, double tolerance, double* outVertexPositions, double* outPlanarityDeviations, double* outGaussianCurvatures);

	
	//--------------------------
	//---- EXTERNAL METHODS FOR CONSTRAINTS
	//--------------------------

	/*! \brief This method makes the vertices specified by the input contatiner fixed.
	*
	*	\param	[in]	_fixedVertices			- input container of anchor point indicies.
	*	\param	[in]	numFixed				- number of fixed vertices.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void ccfSolver_setFixed(int* _fixedVertices, int numFixed);
}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/projectionsolver/zSpace_CCFSolver.cpp>
#endif

#endif
