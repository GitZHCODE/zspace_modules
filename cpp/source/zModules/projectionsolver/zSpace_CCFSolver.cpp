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

#include <headers/zModules/projectionsolver/zSpace_CCFSolver.h>

namespace  zSpace
{

	//---- EXTERNAL METHODS FOR CCF

	ZSPACE_INLINE bool ccfSolver_initialise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, double* outPlanarityDeviations, double* outGaussianCurvatures)
	{
		bool out = false;

		if (_vertexPositions && _polyCounts && _polyConnects && outPlanarityDeviations && outGaussianCurvatures)
		{
			createComputeMesh(_vertexPositions, _polyCounts, _polyConnects, numVerts, numFaces, true, ccfMesh);
		
			// set planarity type
			ccf_planarisationType = zQuadPlanar;

			// compute planarity deviations
			zDoubleArray fDeviations;
			computeQuadPlanarityDeviation(ccfMesh, fDeviations);

			//compute gaussian curvature
			zDoubleArray vGauss;

			zObjMesh oMesh;
			zFnMesh fnMesh(oMesh);
			
			

			//computeGaussianCurvatures(compMesh,)

			// update deviations
			for (int i = 0; i < fDeviations.size(); i++)
			{
				outPlanarityDeviations[i] = fDeviations[i];
			}
			

			out = true;
		}

		return out;
	}

	//---- EXTERNAL METHODS FOR CONSTRAINTS

	ZSPACE_INLINE void ccfSolver_setFixed(int* _fixedVertices, int numFixed)
	{
		if (numFixed >= ccfMesh.vertexPositions.size()) throw std::invalid_argument(" error: number of fixed vertices greater than number of vertices.");

		// set fixed
		if (_fixedVertices)
		{			
			for (int i = 0; i < numFixed; i++)
			{
				ccfMesh.fnParticles[_fixedVertices[i]].setFixed(true);
			}
		}
	}


}

