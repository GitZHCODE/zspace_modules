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



#ifndef ZSPACE_MODULES_SOLVERUTILS_H
#define ZSPACE_MODULES_SOLVERUTILS_H



#pragma once
#include <vector>
#include <string>
using namespace std;


#include <headers/zModules/base/zSpace_MeshUtilities.h>


namespace  zSpace
{

	//--------------------------
	//---- GLOBAL VARIABLES
	//--------------------------
	
	/*!	\brief global solver variable. A container of booleans to indicate if the vertex is fixed or not. */
	extern zBoolArray vFixed;

	/*!	\brief global solver variable. Timestep variable of the solver. */
	extern float dT;

	/*!	\brief global solver variable. Integration type of the solver. */
	extern zIntergrationType intType;

	//--------------------------
	//----  EXTERNAL METHODS
	//--------------------------

	/*! \brief This method makes the global compute mesh vertices specified by the input contatiner fixed.
	*
	*	\param	[in]	_fixedVertices			- input container of anchor point indicies.
	*	\param	[in]	numFixed				- number of fixed vertices.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void computeMesh_setFixed(int* _fixedVertices, int numFixed);

	/*! \brief This method sets the solvers time step variable.
	*
	*	\param	[in]	_dT						- input timestep variable. Should be between 0 & 1.
	*	\param	[in]	numFixed				- number of fixed vertices.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void computeMesh_setDT(double _dT);

	/*! \brief This method sets the mass of the global compute mesh vertices specified by the input contatiner.
	*
	*	\param	[in]	_vMass					- input container of vertex mass.
	*	\param	[in]	numVerts				- number of vertices. Should match the number of vertices of the mesh.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void computeMesh_setMass(int* _vMass, int numVerts);

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/base/zSpace_SolverUtilities.cpp>
#endif

#endif
