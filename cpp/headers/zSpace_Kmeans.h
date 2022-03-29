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



#ifndef ZSPACE_MODULES_KMEANS_H
#define ZSPACE_MODULES_KMEANS_H



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


	extern "C" ZSPACE_MODULES void kMeansClustering(double* _data, int datastride, int numClusters, int numIterations, double* outClusters, double* outClusterCentroids);

}


#endif
