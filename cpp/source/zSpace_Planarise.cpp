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

#include <headers/zSpace_Planarise.h>

namespace  zSpace
{
	ZSPACE_INLINE double quadPlanarityDeviation(int index, zPointArray &_vertexPositions, zInt2DArray& _polygons, zUtilsCore &_core)
	{
		
		double uA, uB;
		zPoint pA, pB;
			
		bool check = _core.line_lineClosestPoints(_vertexPositions[_polygons[index][0]], _vertexPositions[_polygons[index][2]], _vertexPositions[_polygons[index][1]], _vertexPositions[_polygons[index][3]], uA, uB, pA, pB);
					
		return  pA.distanceTo(pB);		
		
	}

	ZSPACE_INLINE void quadPlanarise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, bool updatetopology, int numIterations, double tolerance, double* outVertexPositions, double* outDeviations)
	{
		if (_vertexPositions && _polyCounts && _polyConnects && outVertexPositions && outDeviations)
		{
			zUtilsCore core;

			zDoubleArray deviations;
			deviations.assign(numFaces, 10000);

			if (updatetopology)
			{
				vertexPositions.clear();
				vertexPositions.assign(numVerts, zPoint());

				for (int i = 0; i < numVerts; i += 1)
				{
					vertexPositions[i].x = _vertexPositions[i * 3 + 0];
					vertexPositions[i].y = _vertexPositions[i * 3 + 1];
					vertexPositions[i].z = _vertexPositions[i * 3 + 2];

				}

				polygons.clear();
				polygons.assign(numFaces, vector<int>());

				int polyconnectsCurrentIndex = 0;
				for (int i = 0; i < numFaces; i++)
				{
					int num_faceVerts = _polyCounts[i];

					for (int j = 0; j < num_faceVerts; j++)
					{
						polygons[i].push_back(_polyConnects[polyconnectsCurrentIndex + j]);
					}

					polyconnectsCurrentIndex += num_faceVerts;
				}

				// create particles
				fnParticles.clear();
				o_Particles.clear();


				for (int i = 0; i < numVerts; i++)
				{
					bool fixed = false;
										

					zObjParticle p;
					p.particle = zParticle(vertexPositions[i], fixed);
					o_Particles.push_back(p);

				}

				for (int i = 0; i < o_Particles.size(); i++)
				{
					fnParticles.push_back(zFnParticle(o_Particles[i]));
				}

			}

			bool exit = false;

			for (int i = 0; i < numIterations; i++)
			{
				
				if (!exit)
				{
					exit = true;

					// compute forces
					for (int j = 0; j < numFaces; j++)
					{
						double uA, uB;
						zPoint pA, pB;

						bool check = core.line_lineClosestPoints(vertexPositions[polygons[j][0]], vertexPositions[polygons[j][2]], vertexPositions[polygons[j][1]], vertexPositions[polygons[j][3]], uA, uB, pA, pB);

						deviations[j] = pA.distanceTo(pB);

						if (deviations[j] > tolerance)
						{
							exit = false;

							zVector dir = pB - pA;
							dir.normalize();

							zVector pForceA = dir * deviations[j] * 0.5;
							zVector pForceB = dir * deviations[j] * -0.5;

							fnParticles[polygons[j][0]].addForce(pForceA);
							fnParticles[polygons[j][2]].addForce(pForceA);

							fnParticles[polygons[j][1]].addForce(pForceB);
							fnParticles[polygons[j][3]].addForce(pForceB);
						}

					}

					// update positions
					for (int i = 0; i < fnParticles.size(); i++)
					{
						fnParticles[i].integrateForces(0.5, zIntergrationType::zRK4);
						fnParticles[i].updateParticle(true);
					}
				}			

			}
			
			// output

			for (int i = 0; i < vertexPositions.size(); i++)
			{
				outVertexPositions[i * 3 + 0] = vertexPositions[i].x;
				outVertexPositions[i * 3 + 1] = vertexPositions[i].y;
				outVertexPositions[i * 3 + 2] = vertexPositions[i].z;
			}

			for (int i = 0; i < deviations.size(); i++)
			{
				outDeviations[i] = deviations[i];
			}


			
		}

	}

	ZSPACE_INLINE void planarise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, bool updatetopology, int numIterations, double tolerance, double* outVertexPositions, double* outDeviations)
	{
		
	}

}

