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

#ifndef ZSPACE_MODULES_PROJECTIONFORCES_H
#define ZSPACE_MODULES_PROJECTIONFORCES_H

#pragma once

#include <headers/zModules/base/zSpace_MeshUtilities.h>

namespace  zSpace
{

	//--------------------------
	//---- PROJECTION FORCE METHODS
	//--------------------------

	ZSPACE_MODULES void addPlanarityForces(zComputeMesh &inMesh, zPlanarType type, double &tolerance, bool &exit , zDoubleArray &planarityDeviations, zPointArray &targetCenters, zVectorArray &targetNormals);
	
	//ZSPACE_MODULES void addGaussianForces();

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/projectionsolver/zSpace_ProjectionForces.cpp>
#endif

#endif
