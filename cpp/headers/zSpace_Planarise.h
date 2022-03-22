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

	zPointArray vertexPositions;

	zInt2DArray polygons;

	//--------------------------
	//---- METHOD
	//--------------------------
	
	ZSPACE_MODULES double quadPlanarityDeviation(int index, zPointArray &_vertexPositions, zInt2DArray& _polygons, zUtilsCore& _core);
	
	extern "C" ZSPACE_MODULES void quadPlanarise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, bool updatetopology, int numIterations, double tolerance, double* outVertexPositions,  double* outDeviations);

	extern "C" ZSPACE_MODULES void planarise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, bool updatetopology, int numIterations, double tolerance, double* outVertexPositions, double* outDeviations);

}

#if defined(ZSPACE_STATIC_LIBRARY)  || defined(ZSPACE_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
//#include<zSpace_Planarise.cpp>
#endif

#endif
