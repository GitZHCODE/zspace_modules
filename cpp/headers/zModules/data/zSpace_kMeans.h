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

#include <headers/zModules/base/zSpace_Modules.h>

namespace  zSpace
{

	//--------------------------
	//---- GLOBAL VARIABLES
	//--------------------------

	//--------------------------
	//---- METHOD
	//--------------------------

	/*! \brief This method computes & classifies the input data into input number of clusters using the K-Means Algorithm.
	*
	*	\param	[in]	_data					- input container of data. Collapsed 1D array of size numDataPoints * datastride. 
	*	\param	[in]	numDataPoints			- input number of data points.
	*	\param	[in]	datastride				- input stride of data ie number of features in a datapoint.
	* 	\param	[in]	numClusters				- input number of clusters.
	*  	\param	[in]	numIterations			- input number of iterations.
	*  	\param	[out]	outClusters				- output cluster index per data point.
	*  	\param	[out]	outClusterCentroids		- output cluster centroids/means.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void kMeansClustering(double* _data, int numDataPoints, int datastride, int numClusters, int numIterations, int* outClusters, double* outClusterCentroids);

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/data/zSpace_kMeans.cpp>
#endif

#endif
