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

	ZSPACE_INLINE int planariseSolver_initialise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, bool volPlanarise, double* outDeviations)
	{
		bool out = false;

		if (_vertexPositions && _polyCounts && _polyConnects && outDeviations)
		{
			createComputeMesh(_vertexPositions, _polyCounts, _polyConnects, numVerts, numFaces,true, planariseMesh);

			if (_triCounts && _triConnects) setTriangles(planariseMesh, numFaces, _triCounts, _triConnects);
		
			
			// if volumetric planarise add triangles to the compute mesh
			if (volPlanarise)
			{
				// set planarity type
				planarisationType = zVolumePlanar;

				zObjMesh oMesh;
				createMeshOBJ(_vertexPositions, _polyCounts, _polyConnects, numVerts, numFaces, oMesh);

				if (!_triCounts || !_triConnects)
				{
					zFnMesh fnMesh(oMesh);
					fnMesh.getMeshTriangles(planariseMesh.triangles);
				}										

				// compute face normals
				zVectorArray fNormals;
				computeFaceNormals(planariseMesh, fNormals);

				zPointArray fCenters;
				zDoubleArray fVolumes;
				
				// compute volumes
				computeFaceVolumes(planariseMesh, fCenters, fVolumes);

				// update deviations
				for (int i = 0; i < fVolumes.size(); i++)
				{
					outDeviations[i] = fVolumes[i];
				}											

			}
			else
			{	
				// set planarity type
				planarisationType = zQuadPlanar;

				// compute planarity deviations
				zDoubleArray fDeviations;
				computeQuadPlanarityDeviation(planariseMesh, fDeviations);

				// update deviations							
				for (int i = 0; i < fDeviations.size(); i++)
				{
					outDeviations[i] = fDeviations[i];
				}
			}


			out = true;
		}

		return out;
	}

	ZSPACE_INLINE void planariseSolver_compute(int numIterations, double tolerance, double* outVertexPositions, double* outDeviations)
	{
		zDoubleArray fDeviations;

		zPointArray fCenters;
		zVectorArray fNormals;

		if (planarisationType == zVolumePlanar)
		{
			bool exit = false;

			// compute face normals			
			computeFaceNormals(planariseMesh, fNormals);
						
			for (int i = 0; i < numIterations; i++)
			{
				// compute volumes
				computeFaceVolumes(planariseMesh, fCenters, fDeviations);

				if (!exit)
				{
					exit = true;

					// add planarity forces
					addPlanarityForces(planariseMesh, planarisationType, tolerance, exit, fDeviations, fCenters, fNormals);

					// update positions
					for (int i = 0; i < planariseMesh.fnParticles.size(); i++)
					{
						planariseMesh.fnParticles[i].integrateForces(0.5, zIntergrationType::zRK4);
						planariseMesh.fnParticles[i].updateParticle(true);
					}
				}

				computeFaceNormals(planariseMesh, fNormals);

			}
		
		}
		else if (planarisationType == zQuadPlanar)
		{
			bool exit = false;
		
			for (int i = 0; i < numIterations; i++)
			{
				// compute deviations
				computeQuadPlanarityDeviation(planariseMesh, fDeviations);

				if (!exit)
				{
					exit = true;

					// compute forces
					addPlanarityForces(planariseMesh, planarisationType, tolerance, exit, fDeviations, fCenters, fNormals);

					// update positions
					for (int i = 0; i < planariseMesh.fnParticles.size(); i++)
					{
						planariseMesh.fnParticles[i].integrateForces(0.5, zIntergrationType::zRK4);
						planariseMesh.fnParticles[i].updateParticle(true);
					}
				}

			}

		}
		
		// output

		for (int i = 0; i < planariseMesh.vertexPositions.size(); i++)
		{
			outVertexPositions[i * 3 + 0] = planariseMesh.vertexPositions[i].x;
			outVertexPositions[i * 3 + 1] = planariseMesh.vertexPositions[i].y;
			outVertexPositions[i * 3 + 2] = planariseMesh.vertexPositions[i].z;
		}

		for (int i = 0; i < fDeviations.size(); i++)
		{
			outDeviations[i] = fDeviations[i];
		}		

	}

	//---- EXTERNAL METHODS FOR CONSTRAINTS

	ZSPACE_INLINE void planariseSolver_setFixed(int* _fixedVertices, int numFixed)
	{
		if (numFixed >= planariseMesh.vertexPositions.size()) throw std::invalid_argument(" error: number of fixed vertices greater than number of vertices.");

		// set fixed
		if (_fixedVertices)
		{
			for (int i = 0; i < numFixed; i++)
			{
				planariseMesh.fnParticles[_fixedVertices[i]].setFixed(true);
			}
		}
	}

}

