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

#ifndef ZSPACE_MODULES_FDM_H
#define ZSPACE_MODULES_FDM_H

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
	//---- METHODS
	//--------------------------

	/*! \brief This method computes the Edge Node Matrix for the input mesh.
	*
	*	\param		[in]	numRows								- number of rows in the out matrix.
	*	\return				zSparseMatrix						- edge node matrix.
	*	\since version 0.0.4
	*/	
	zSparseMatrix getEdgeNodeMatrix(int numRows);

	/*! \brief This method computes the sub Matrix of a sparse matrix.
	*
	*	\param		[in]	C									- input sparse matrix.
	*	\param		[in]	nodes								- container of integers.
	*	\return				zSparseMatrix								- sub matrix.
	*	\since version 0.0.4
	*/
	zSparseMatrix subMatrix(zSparseMatrix& C, vector<int>& nodes);

	/*! \brief This method computes the sub Matrix of a matrix.
	*
	*	\param		[in]	C									- input sparse matrix.
	*	\param		[in]	nodes								- container of integers.
	*	\return				MatrixXd							- sub matrix.
	*	\since version 0.0.4
	*/
	MatrixXd subMatrix(MatrixXd& X, vector<int>& nodes);

	//--------------------------
	//---- EXTERNAL METHODS FOR FDM
	//--------------------------

	/*! \brief This method computes the result based on the force density method.
	*
	*	\details Based on Schek, H-J. "The force density method for form finding and computation of general networks." Computer methods in applied mechanics and engineering 3.1 (1974): 115-134. (https://www.sciencedirect.com/science/article/pii/0045782574900450)
				and Linkwitz, K. (2014). Force density method. Shell Structures for Architecture: Form Finding and Optimization, Routledge, Oxon and New York, 59-71.
	*  	\param	[in]	forceDensities			- input container of force densities per vertex.
	*  	\param	[in]	tolerance				- input container of vertex mass per vertex.
	*  	\param	[out]	outVertexPositions		- output container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\return			int						- output boolean - true if fdm is successful.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES int computeMesh_fdm(double *vForceDensities, double* vMass, double* outVertexPositions);

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/statics/zSpace_FDM.cpp>
#endif

#endif
