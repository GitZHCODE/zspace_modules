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

	ZSPACE_INLINE void computeFaceVolumes(zPointArray& _vertexPositions, zInt2DArray& _polygons, zInt2DArray& _triangles, zUtilsCore& _core, zPointArray& fCenters, zFloatArray& fVolumes)
	{
		fVolumes.clear();
		fVolumes.assign(_polygons.size(), float());

		fCenters.clear();
		fCenters.assign(_polygons.size(), zPoint());

		for (int i = 0; i < _polygons.size(); i++)
		{
			// compute face center
			for (int j = 0; j < _polygons[i].size(); j++) fCenters[i] += _vertexPositions[_polygons[i][j]];
			fCenters[i] /= _polygons[i].size();

			int nTris = floor(triangles[i].size() / 3);

			// add volume of face tris	
			for (int j = 0; j < triangles[i].size(); j += 3)
			{
				double vol = _core.getSignedTriangleVolume(_vertexPositions[_triangles[i][j + 0]], _vertexPositions[_triangles[i][j + 1]], _vertexPositions[_triangles[i][j + 2]]);
				fVolumes[i] += vol;
			}

			// add volumes of tris formed by each pair of face edge vertices and face center
			for (int j = 0; j < _polygons[i].size(); j++)
			{
				int prevId = (j - 1 + _polygons[i].size()) % _polygons[i].size();

				double vol = _core.getSignedTriangleVolume(_vertexPositions[_polygons[i][j]], _vertexPositions[_polygons[i][prevId]], fCenters[i]);
				fVolumes[i] += vol;
			}

			//absolute value
			fVolumes[i] = abs(fVolumes[i]);
		}
	}

	ZSPACE_INLINE void computeFaceNormals(zPointArray& _vertexPositions, zInt2DArray& _polygons, zInt2DArray& _triangles, zUtilsCore& _core, zVectorArray& fNormals)
	{
		fNormals.clear();
		fNormals.assign(_polygons.size(), zVector());

		for (int i = 0; i < _polygons.size(); i++)
		{
			int nTris = floor(triangles[i].size() / 3);
						
			for (int j = 0; j < triangles[i].size(); j += 3)
			{
				zVector cross = (_vertexPositions[_triangles[i][j + 1]] - _vertexPositions[_triangles[i][j + 0]]) ^ (_vertexPositions[_triangles[i][j + 2]] - _vertexPositions[_triangles[i][j + 0]]);
				cross.normalize();

				fNormals[i] += cross;
			}

			fNormals[i] /= nTris;
			fNormals[i].normalize();
		}
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
		zUtilsCore core;

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

			triangles.clear();
			triangles.assign(numFaces, vector<int>());

			int triconnectsCurrentIndex = 0;
			for (int i = 0; i < numFaces; i++)
			{
				int num_triVerts = _triCounts[i];

				for (int j = 0; j < num_triVerts * 3; j++)
				{
					triangles[i].push_back(_triConnects[triconnectsCurrentIndex + j]);
				}

				triconnectsCurrentIndex += num_triVerts * 3;
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

		// compute face normals
		zVectorArray fNormals;
		computeFaceNormals(vertexPositions, polygons, triangles, core, fNormals);

		zPointArray fCenters;
		zFloatArray fVolumes;

		for (int i = 0; i < numIterations; i++)
		{
			// compute volumes
			computeFaceVolumes(vertexPositions, polygons, triangles, core, fCenters, fVolumes);

			if (!exit)
			{
				exit = true;

				// compute forces
				for (int j = 0; j < numFaces; j++)
				{
					

					if (fVolumes[j] > tolerance)
					{						
						exit = false;

						for (int k = 0; k < polygons[j].size(); k++)
						{
							double dist = core.minDist_Point_Plane(vertexPositions[polygons[j][k]], fCenters[j], fNormals[j]);
							zVector pForce = fNormals[j] * dist * -1.0;
							fnParticles[polygons[j][k]].addForce(pForce);

						}						
					}

				}

				// update positions
				for (int i = 0; i < fnParticles.size(); i++)
				{
					fnParticles[i].integrateForces(0.5, zIntergrationType::zRK4);
					fnParticles[i].updateParticle(true);
				}
			}

			computeFaceNormals(vertexPositions, polygons, triangles, core, fNormals);

		}

		// output

		for (int i = 0; i < vertexPositions.size(); i++)
		{
			outVertexPositions[i * 3 + 0] = vertexPositions[i].x;
			outVertexPositions[i * 3 + 1] = vertexPositions[i].y;
			outVertexPositions[i * 3 + 2] = vertexPositions[i].z;
		}

		for (int i = 0; i < fVolumes.size(); i++)
		{
			outDeviations[i] = fVolumes[i];
		}
	}
}

