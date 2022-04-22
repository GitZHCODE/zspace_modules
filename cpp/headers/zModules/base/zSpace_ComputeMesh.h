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



#ifndef ZSPACE_MODULES_MESHCONTAINERS_H
#define ZSPACE_MODULES_MESHCONTAINERS_H



#pragma once
#include <vector>
#include <string>
using namespace std;

#include <headers/zInterface/functionsets/zFnParticle.h>
#include <headers/zModules/base/zSpace_Modules.h>

namespace  zSpace
{
	struct zComputeMesh
	{
		/*!	\brief container of  particle objects  */
		vector<zObjParticle> o_Particles;

		/*!	\brief container of particle function set  */
		vector<zFnParticle> fnParticles;

		/*!	\brief container of vertex positions  */
		zPointArray vertexPositions;

		/*!	\brief 2D container of polygon vertex connectivity  */
		zInt2DArray polygons;

		/*!	\brief 2D container of edge vertex connectivity  */
		zInt2DArray edges;			

		/*!	\brief 2D container of triangle vertex connectivity */
		zInt2DArray triangles;	

		/*!	\brief matrix of vertex positions. This is to be set to call IGL methods */
		MatrixXd V;

		/*!	\brief matrix of Tri-Face connectivity. This is to be set to call IGL methods */
		MatrixXi FTris;

		/*!	\brief number of vertices */
		int nV;

		/*!	\brief number of edges */
		int nE;

		/*!	\brief number of polygons */
		int nF;

		/*!	\brief number of triangles */
		int nT;

	};

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
// All defined OK so do nothing
#endif

#endif
