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
#include <headers/zModules/base/zSpace_SolverUtilities.h>

#include <headers/zModules/projectionsolver/zSpace_ProjectionForces.h>

namespace  zSpace
{

	//--------------------------
	//---- EXTERNAL METHODS FOR PLANARISATION
	//--------------------------

	/*! \brief This method planarises the global compute mesh.
	*
	*  	\param	[in]	numIterations			- input number of iterations to run the solver.
	*  	\param	[in]	tolerance				- input planarity tolerance.
	*  	\param	[in]	volPlanarise			- input boolean to set planarity solver type - quad diagonal distance or volume.
	*  	\param	[out]	outVertexPositions		- output container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*  	\param	[out]	outDeviations			- output container of planarity deviations per face/polygon.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void computeMesh_planarise(int numIterations, bool volPlanarise, double tolerance, double* outVertexPositions, double* outDeviations);

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/projectionsolver/zSpace_PlanariseSolver.cpp>
#endif

#endif
