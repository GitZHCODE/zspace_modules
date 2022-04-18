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

#include <headers/zModules/projectionsolver/zSpace_MinimalSurfaceSolver.h>

namespace  zSpace
{

	//---- EXTERNAL METHODS FOR Minimal Surface

	ZSPACE_MODULES_INLINE bool msSolver_initialise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, double* outMeanCurvature)
	{
		bool out = false;

		if (_vertexPositions && _polyCounts && _polyConnects && _triCounts && _triConnects && outMeanCurvature)
		{
			createComputeMesh(_vertexPositions, _polyCounts, _polyConnects, numVerts, numFaces, true, msMesh);

			// set triangles and matrices for igl method call
			setTriangles(msMesh, numFaces, _triCounts, _triConnects);
			updateMatrixV(msMesh);
			updateMatrixFTris(msMesh);
			
			//compute mean curvature
			VectorXd vMeanCurvature;
			MatrixXd vMeanNormal;
			computeMeanCurvature(msMesh, vMeanCurvature, vMeanNormal);

			// update deviations
			for (int i = 0; i < msMesh.nV; i++)
			{
				outMeanCurvature[i] = vMeanCurvature[i];
			}


			out = true;
		}

		return out;
	}

	//---- EXTERNAL METHODS FOR CONSTRAINTS

	ZSPACE_MODULES_INLINE void msSolver_setFixed(int* _fixedVertices, int numFixed)
	{
		setFixed(msMesh, _fixedVertices, numFixed);
	}

}

