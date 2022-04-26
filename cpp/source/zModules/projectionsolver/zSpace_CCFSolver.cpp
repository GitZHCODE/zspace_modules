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

	ZSPACE_MODULES_INLINE bool ccfSolver_initialise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, double* outPlanarityDeviations, double* outGaussianCurvatures)
	{
		bool out = false;

		if (_vertexPositions && _polyCounts && _polyConnects && _triCounts && _triConnects && outPlanarityDeviations && outGaussianCurvatures)
		{
			createComputeMesh(_vertexPositions, _polyCounts, _polyConnects, numVerts, numFaces, true, ccfMesh);
			
			// set triangles and matrices for igl method call
			setTriangles(ccfMesh, numFaces, _triCounts, _triConnects);
			updateMatrixV(ccfMesh);
			updateMatrixFTris(ccfMesh);

			// set planarity type
			ccf_planarisationType = zQuadPlanar;

			// compute planarity deviations
			zDoubleArray fDeviations;
			computeQuadPlanarityDeviation(ccfMesh, fDeviations);

			//compute gaussian curvature
			VectorXd vGauss;
			computeGaussianCurvature(ccfMesh, vGauss);

			// update deviations
			for (int i = 0; i < ccfMesh.nF; i++)
			{
				outPlanarityDeviations[i] = fDeviations[i];
			}

			for (int i = 0; i < ccfMesh.nV; i++)
			{
				outGaussianCurvatures[i] = vGauss[i];
			}
			

			out = true;
		}

		return out;
	}

	//---- EXTERNAL METHODS FOR CONSTRAINTS

	ZSPACE_MODULES_INLINE void ccfSolver_compute(int numIterations, double tolerance, double* outVertexPositions, double* outPlanarityDeviations, double* outGaussianCurvatures)
	{
		zDoubleArray fDeviations;

		zPointArray fCenters;
		zVectorArray fNormals;

		bool exit = false;

		for (int i = 0; i < numIterations; i++)
		{
			// compute deviations
			computeQuadPlanarityDeviation(ccfMesh, fDeviations);

			if (!exit)
			{
				//exit = true;
				// 
				// add planarity forces
				bool pExit = true;
				addPlanarityForces(ccfMesh, zQuadPlanar, fCenters, fNormals, tolerance, fDeviations, exit);

				// add gaussian forces
				bool gExit = true;
				//addGaussianForces();

				if (pExit && gExit) exit = true;

				// update positions
				for (int i = 0; i < ccfMesh.fnParticles.size(); i++)
				{
					ccfMesh.fnParticles[i].integrateForces(0.5, zIntergrationType::zRK4);
					ccfMesh.fnParticles[i].updateParticle(true);
				}
			}
		}

		// output

		for (int i = 0; i < ccfMesh.vertexPositions.size(); i++)
		{
			outVertexPositions[i * 3 + 0] = ccfMesh.vertexPositions[i].x;
			outVertexPositions[i * 3 + 1] = ccfMesh.vertexPositions[i].y;
			outVertexPositions[i * 3 + 2] = ccfMesh.vertexPositions[i].z;
		}

		for (int i = 0; i < fDeviations.size(); i++)
		{
			outPlanarityDeviations[i] = fDeviations[i];
		}

		VectorXd vGauss;
		computeGaussianCurvature(ccfMesh, vGauss);
		for (int i = 0; i < fDeviations.size(); i++)
		{
			outGaussianCurvatures[i] = vGauss(i);
		}
	}

	ZSPACE_MODULES_INLINE void ccfSolver_setFixed(int* _fixedVertices, int numFixed)
	{
		setFixed(ccfMesh, _fixedVertices, numFixed);
	}


}

