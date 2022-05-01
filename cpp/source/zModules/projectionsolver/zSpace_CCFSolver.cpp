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

	ZSPACE_MODULES_INLINE void computeMesh_ccf(int numIterations, double tolerance, double* outVertexPositions, double* outPlanarityDeviations, double* outGaussianCurvatures)
	{
		zDoubleArray fDeviations;
		VectorXd vGauss;

		zPointArray fCenters;
		zVectorArray fNormals;

		bool exit = false;

		for (int i = 0; i < numIterations; i++)
		{
			// compute deviations
			computeQuadPlanarityDeviation(compMesh, fDeviations);

			if (!exit)
			{
				//exit = true;
				// 
				// add planarity forces
				bool pExit = true;
				addPlanarityForces(compMesh, zQuadPlanar, fCenters, fNormals, tolerance, fDeviations, exit);

				// add gaussian forces
				bool gExit = true;
				//addGaussianForces();

				if (pExit && gExit) exit = true;

				// update positions
				for (int i = 0; i < compMesh.fnParticles.size(); i++)
				{
					compMesh.fnParticles[i].integrateForces(dT, intType);
					compMesh.fnParticles[i].updateParticle(true);
				}
			}
		}

		// output

		// compute planarity deviations
		computeQuadPlanarityDeviation(compMesh, fDeviations);

		//compute gaussian curvature
		computeGaussianCurvature(compMesh, vGauss);

		for (int i = 0; i < compMesh.vertexPositions.size(); i++)
		{
			outVertexPositions[i * 3 + 0] = compMesh.vertexPositions[i].x;
			outVertexPositions[i * 3 + 1] = compMesh.vertexPositions[i].y;
			outVertexPositions[i * 3 + 2] = compMesh.vertexPositions[i].z;
		}

		for (int i = 0; i < fDeviations.size(); i++)
		{
			outPlanarityDeviations[i] = fDeviations[i];
		}
				
		for (int i = 0; i < fDeviations.size(); i++)
		{
			outGaussianCurvatures[i] = vGauss(i);
		}
	}


}

