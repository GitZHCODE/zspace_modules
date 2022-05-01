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

#include <headers/zModules/projectionsolver/zSpace_PlanariseSolver.h>

namespace  zSpace
{

	//---- EXTERNAL METHODS FOR PLANARISATION

	ZSPACE_MODULES_INLINE void computeMesh_planarise(int numIterations, bool volPlanarise, double tolerance, double* outVertexPositions, double* outDeviations)
	{
		zPlanarSolverType planarisationType = (volPlanarise)? zVolumePlanar : zQuadPlanar ;

		zDoubleArray fDeviations;

		zPointArray fCenters;
		zVectorArray fNormals;

		if (planarisationType == zVolumePlanar)
		{
			bool exit = false;

			// compute face normals			
			computeFaceNormals(compMesh, fNormals);
						
			for (int i = 0; i < numIterations; i++)
			{
				// compute volumes
				computeFaceVolumes(compMesh, fCenters, fDeviations);

				if (!exit)
				{
					exit = true;

					// add planarity forces
					addPlanarityForces(compMesh, planarisationType, fCenters, fNormals, tolerance, fDeviations, exit);
										 
					// update positions
					for (int i = 0; i < compMesh.fnParticles.size(); i++)
					{
						compMesh.fnParticles[i].integrateForces(dT, intType);
						compMesh.fnParticles[i].updateParticle(true);
					}
				}

				computeFaceNormals(compMesh, fNormals);

			}

			// compute deviations if number interation is zero
			if(numIterations == 0) computeFaceVolumes(compMesh, fCenters, fDeviations);
		
		}
		else if (planarisationType == zQuadPlanar)
		{
			bool exit = false;
		
			for (int i = 0; i < numIterations; i++)
			{
				// compute deviations
				computeQuadPlanarityDeviation(compMesh, fDeviations);

				if (!exit)
				{
					exit = true;

					// compute forces
					addPlanarityForces(compMesh, planarisationType, fCenters, fNormals, tolerance, fDeviations, exit);

					// update positions
					for (int i = 0; i < compMesh.fnParticles.size(); i++)
					{
						compMesh.fnParticles[i].integrateForces(dT, intType);
						compMesh.fnParticles[i].updateParticle(true);
					}
				}

			}

			// compute deviations if number interation is zero
			if (numIterations == 0) computeQuadPlanarityDeviation(compMesh, fDeviations);

		}
		
		// output

		for (int i = 0; i < compMesh.vertexPositions.size(); i++)
		{
			outVertexPositions[i * 3 + 0] = compMesh.vertexPositions[i].x;
			outVertexPositions[i * 3 + 1] = compMesh.vertexPositions[i].y;
			outVertexPositions[i * 3 + 2] = compMesh.vertexPositions[i].z;
		}

		for (int i = 0; i < fDeviations.size(); i++)
		{
			outDeviations[i] = fDeviations[i];
		}		

	}

}

