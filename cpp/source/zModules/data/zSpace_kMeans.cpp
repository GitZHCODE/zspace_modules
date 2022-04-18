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

#include <headers/zModules/data/zSpace_kMeans.h>

namespace  zSpace
{
	ZSPACE_MODULES_INLINE void kMeansClustering(double* _data, int numDataPoints, int datastride, int numClusters, int numIterations, int* outClusters, double* outClusterCentroids)
	{
		zUtilsCore core;

		// get min max value of the datapoints
		
		vector<zDomainFloat> minMax;
		minMax.assign(datastride, zDomainFloat());

		for (int j = 0; j < datastride; j++)
		{
			minMax[j].min = std::numeric_limits<float>::max();
			minMax[j].max = std::numeric_limits<float>::min();
		}

		for (int i = 0; i < numDataPoints; i++)
		{
			for (int j = 0; j < datastride; j++)
			{
				if (_data[i * datastride + j] < minMax[j].min) minMax[j].min = _data[i * datastride + j];

				if (_data[i * datastride + j] > minMax[j].max) minMax[j].max = _data[i * datastride + j];
			}
			
		}

		// Initialise means
		int meanSize = numClusters * datastride;
		double*  means = new double[meanSize];
		srand(time(NULL));
		for (int i = 0; i < numClusters; i++)
		{				
			for (int j = 0; j < datastride; j++)
				means[i* datastride + j] = core.randomNumber_double(minMax[j].min, minMax[j].max);
			
		}
				
		// Initialise container to store tempMeans, cluster counts and  clusterID
		
		double* tempMeans = new double[meanSize];
		int* clusterCounts = new int[numClusters];	
		int* clusters = new int[numDataPoints];		
		
		
		// Compute means
		for (int k = 0; k < numIterations; k++)
		{
			// reset cluster counts
			for (int i = 0; i < numClusters; i++) clusterCounts[i] = 0;

			for (int i = 0; i < numDataPoints; i++)
			{
				// get cluster index
				double minDist = 10000000;
				int clusterID = -1;

				for (int j = 0; j < numClusters; j++)
				{
					double tempDist = 0;
					for (int l = 0; l < datastride; l++)
					{
						int dataID = (i * datastride) + l;
						int meanID = (j * datastride) + l;

						tempDist += pow(_data[dataID] - means[meanID], 2);
					}

					tempDist = sqrt(tempDist);
										

					if (tempDist < minDist)
					{
						minDist = tempDist;
						clusterID = j;
					}
				}
								
				// store and update cluster, clustercounts
				clusters[i] = clusterID;
				clusterCounts[clusterID]+= 1;

				// update temp mean
				for (int l = 0; l < datastride; l++)
				{
					int dataID = (i * datastride) + l;
					int meanID = (clusterID * datastride) + l;
					tempMeans[meanID] += _data[dataID];
				}
			}

			// update mean
			
			for (int j = 0; j < numClusters; j++)
			{	
				for (int l = 0; l < datastride; l++)
				{
					if(clusterCounts[j] > 0) means[j * datastride + l] = tempMeans[j * datastride + l] / clusterCounts[j];
					tempMeans[j * datastride + l] = 0.0;
				}
			}

		}

		std::copy(clusters, clusters + numDataPoints, outClusters);
		std::copy(means, means + (meanSize), outClusterCentroids);

		// delete dynamic memory
		delete[] means;
		delete[] tempMeans;
		delete[] clusterCounts;
		delete[] clusters;
	
	}
	
}