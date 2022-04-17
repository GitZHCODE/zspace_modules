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

		/*!	\brief 2D container of triangle vertex connectivity */
		zInt2DArray triangles;		

	};

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
// All defined OK so do nothing
#endif

#endif
