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

#ifndef ZSPACE_MODULES_SPATIALBIN_H
#define ZSPACE_MODULES_SPATIALBIN_H

#pragma once
#include <vector>
#include <string>
using namespace std;

#include <headers/zModules/base/zSpace_SolverUtilities.h>

namespace  zSpace
{
	//--------------------------
	//---- GLOBAL VARIABLES
	//--------------------------
	
	/*!	\brief global variable for binning. container of bin positions */
	extern zPointArray bin_positions;

	/*!	\brief global variable for binning. container of vertex indicies inside the bin per object */
	extern zInt2DArray bin_vIds;
	
	/*!	\brief global variable for binning. Resolution of Spatial Bin */
	extern int bin_RES;

	/*!	\brief global variable for binning. domain of bin bounds */
	extern zDomainVector bin_bounds;

	/*!	\brief global variable for binning. unit vector of bin */
	extern zVector bin_unitVec;

	//--------------------------
	//----  METHODS
	//--------------------------
	
	/*! \brief This method sets the mass of the global compute mesh vertices specified by the input contatiner.
	*
	*	\param	[in]	_vMass					- input container of vertex mass.
	*	\param	[in]	numVerts				- number of vertices. Should match the number of vertices of the mesh.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void spatialBin(zComputeMesh &inMesh, int _res);

	/*! \brief This method sets the mass of the global compute mesh vertices specified by the input contatiner.
	*
	*	\param	[in]	_vMass					- input container of vertex mass.
	*	\param	[in]	numVerts				- number of vertices. Should match the number of vertices of the mesh.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES bool spatialBin_getID(zPoint& inPoint, int &outID);

	ZSPACE_MODULES void spatialBin_getNeighbourRing(int binID, zIntArray& ringNeighbours);
}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/base/zSpace_SpatialBinUtilities.cpp>
#endif

#endif
