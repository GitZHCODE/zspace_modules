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

	ZSPACE_MODULES_INLINE int computeMesh_minSrf(int numIterations, bool minAreaSolver, double tolerance, double* outVertexPositions, double* outMeanCurvatures)
	{
		bool exit = false;

		zMSSolverType minimiseType = (minAreaSolver) ? zMinimiseArea : minimiseType;

		VectorXd vMeanCurvature;
		MatrixXd vMeanNormal;

		float springconstant = 1;		

		zFloatArray restLengths;
		restLengths.assign(compMesh.nE, 0.01);

		
		for (int i = 0; i < numIterations; i++)
		{
			if (!exit)
			{
				//// Minimize Area Method
				if (minimiseType == zMinimiseArea)
				{
					exit = true;
					addMinimizeAreaForces(compMesh, tolerance, vMeanCurvature, exit);
				}
								
				//Restlength Relaxation Method			
				if (minimiseType == zRestlength)
				{
					addSpringForce(compMesh, restLengths, springconstant);
				}

				// update positions
				for (int i = 0; i < compMesh.fnParticles.size(); i++)
				{					
					compMesh.fnParticles[i].integrateForces(dT, intType);
					compMesh.fnParticles[i].updateParticle(true,true,true);
				}

				// update matrix positions
				updateMatrixV(compMesh);
			}
		}

		// output
		computeMeanCurvature(compMesh, vMeanCurvature, vMeanNormal);

		for (int i = 0; i < compMesh.vertexPositions.size(); i++)
		{
			outVertexPositions[i * 3 + 0] = compMesh.vertexPositions[i].x;
			outVertexPositions[i * 3 + 1] = compMesh.vertexPositions[i].y;
			outVertexPositions[i * 3 + 2] = compMesh.vertexPositions[i].z;
		}

		exit = true;
		for (int i = 0; i < compMesh.nV; i++)
		{
			outMeanCurvatures[i] = vMeanCurvature(i);

			if (vMeanCurvature(i) > tolerance) exit = false;
		}

		return exit;
	}



}

