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

#include <headers/zSpace_Modules.h>

namespace  zSpace
{

	//--------------------------
	//---- GLOBAL VARIABLES
	//--------------------------

	//--------------------------
	//---- METHOD
	//--------------------------

	/*! \brief This method constructs a zMesh from from float arrays.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	*  	\param	[in]	out_mesh				- output z mesh.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void ConstructTopology(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, zObjMesh& out_mesh);

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zSpace_Factory.cpp>
#endif

#endif
