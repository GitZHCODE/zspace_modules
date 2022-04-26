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

	ZSPACE_MODULES_INLINE bool msSolver_initialise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, bool minAreaSolver, double* outMeanCurvature)
	{
		bool out = false;

		if (_vertexPositions && _polyCounts && _polyConnects && _triCounts && _triConnects && outMeanCurvature)
		{
			createComputeMesh(_vertexPositions, _polyCounts, _polyConnects, numVerts, numFaces, true, msMesh);

			// set triangles and matrices for igl method call
			setTriangles(msMesh, numFaces, _triCounts, _triConnects);
			updateMatrixV(msMesh);
			updateMatrixFTris(msMesh);
			
			// set ms solver type
			if (minAreaSolver) 	minimiseType = zMinimiseArea;
			else minimiseType = zRestlength;
			
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

	ZSPACE_MODULES_INLINE int msSolver_compute(int numIterations, double tolerance, double* outVertexPositions, double* outMeanCurvatures)
	{
		bool exit = false;

		VectorXd vMeanCurvature;
		MatrixXd vMeanNormal;

		float springconstant = 1;		

		zFloatArray restLengths;
		restLengths.assign(msMesh.nE, 0.001);

		for (int i = 0; i < msMesh.nV; i++)
		{
			msMesh.fnParticles[i].setMass(1);
		}
		
		for (int i = 0; i < numIterations; i++)
		{
			if (!exit)
			{
				//// Minimize Area Method
				if (minimiseType == zMinimiseArea)
				{
					exit = true;
					addMinimizeAreaForces(msMesh, tolerance, vMeanCurvature, exit);
				}
								
				//Restlength Relaxation Method			
				if (minimiseType == zRestlength)
				{
					addSpringForce(msMesh, restLengths, springconstant);
				}

				// update positions
				for (int i = 0; i < msMesh.fnParticles.size(); i++)
				{					
					msMesh.fnParticles[i].integrateForces(0.1, zIntergrationType::zRK4);
					msMesh.fnParticles[i].updateParticle(true,true,true);
				}

				// update matrix positions
				updateMatrixV(msMesh);
			}
		}

		// output
		computeMeanCurvature(msMesh, vMeanCurvature, vMeanNormal);

		for (int i = 0; i < msMesh.vertexPositions.size(); i++)
		{
			outVertexPositions[i * 3 + 0] = msMesh.vertexPositions[i].x;
			outVertexPositions[i * 3 + 1] = msMesh.vertexPositions[i].y;
			outVertexPositions[i * 3 + 2] = msMesh.vertexPositions[i].z;
		}

		for (int i = 0; i < msMesh.nV; i++)
		{
			outMeanCurvatures[i] = vMeanCurvature(i);

			if (vMeanCurvature(i) < tolerance) exit = false;
		}



		return exit;
	}

	//---- EXTERNAL METHODS FOR CONSTRAINTS

	ZSPACE_MODULES_INLINE void msSolver_setFixed(int* _fixedVertices, int numFixed)
	{	

		setFixed(msMesh, _fixedVertices, numFixed);
	}

}

