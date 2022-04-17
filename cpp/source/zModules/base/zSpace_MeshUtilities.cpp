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

#include <headers/zModules/base/zSpace_MeshUtilities.h>


namespace  zSpace
{
	//----  CREATE METHODS
	
	ZSPACE_INLINE void createMeshOBJ(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, zObjMesh& out_mesh)
	{
		if (!_vertexPositions || !_polyCounts || !_polyConnects) throw std::invalid_argument(" error: mesh container is empty.");

		zPointArray vPos;
		zIntArray pConnects;
		zIntArray pCounts;

		for (int i = 0; i < numVerts; i++)
		{
			zVector v;
			v = zVector(_vertexPositions[i * 3 + 0], _vertexPositions[i * 3 + 1], _vertexPositions[i * 3 + 2]) ;
			vPos.push_back(v);
		}


		int polyconnectsCurrentIndex = 0;
		for (int i = 0; i < numFaces; i++)
		{
			int num_faceVerts = _polyCounts[i];
			pCounts.push_back(_polyCounts[i]);

			for (int j = 0; j < num_faceVerts; j++)
			{
				pConnects.push_back(_polyConnects[polyconnectsCurrentIndex + j]);
			}

			polyconnectsCurrentIndex += num_faceVerts;
		}


		zFnMesh fnMesh(out_mesh);
		fnMesh.create(vPos, pCounts, pConnects);
		
	}

	ZSPACE_INLINE void createComputeMesh(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, bool makeDynamic, zComputeMesh& outMesh)
	{
		if (!_vertexPositions || !_polyCounts || !_polyConnects) throw std::invalid_argument(" error: mesh container is empty.");

		outMesh.vertexPositions.clear();
		outMesh.vertexPositions.assign(numVerts, zPoint());

		for (int i = 0; i < numVerts; i += 1)
		{
			outMesh.vertexPositions[i].x = _vertexPositions[i * 3 + 0];
			outMesh.vertexPositions[i].y = _vertexPositions[i * 3 + 1];
			outMesh.vertexPositions[i].z = _vertexPositions[i * 3 + 2];

		}

		outMesh.polygons.clear();
		outMesh.polygons.assign(numFaces, vector<int>());

		int polyconnectsCurrentIndex = 0;
		for (int i = 0; i < numFaces; i++)
		{
			int num_faceVerts = _polyCounts[i];

			for (int j = 0; j < num_faceVerts; j++)
			{
				outMesh.polygons[i].push_back(_polyConnects[polyconnectsCurrentIndex + j]);
			}

			polyconnectsCurrentIndex += num_faceVerts;
		}

		// particles
		outMesh.fnParticles.clear();
		outMesh.o_Particles.clear();

		if (makeDynamic)
		{
			for (int i = 0; i < numVerts; i++)
			{
				bool fixed = false;


				zObjParticle p;
				p.particle = zParticle(outMesh.vertexPositions[i], fixed);
				outMesh.o_Particles.push_back(p);

			}

			for (int i = 0; i < outMesh.o_Particles.size(); i++)
			{
				outMesh.fnParticles.push_back(zFnParticle(outMesh.o_Particles[i]));
			}
		}
		
		outMesh.triangles.clear();		
	}

	//----  SET METHODS

	ZSPACE_INLINE void setTriangles(zComputeMesh& inMesh, int numFaces, int* _triCounts, int* _triConnects)
	{
		if (numFaces != inMesh.polygons.size()) throw std::invalid_argument(" error: triangles container is empty.");
		if(!_triCounts || !_triConnects) throw std::invalid_argument(" error: triangles container is empty.");

		inMesh.triangles.clear();
		inMesh.triangles.assign(numFaces, vector<int>());

		int triconnectsCurrentIndex = 0;
		for (int i = 0; i < numFaces; i++)
		{
			int num_triVerts = _triCounts[i];

			for (int j = 0; j < num_triVerts * 3; j++)
			{
				inMesh.triangles[i].push_back(_triConnects[triconnectsCurrentIndex + j]);
			}

			triconnectsCurrentIndex += num_triVerts * 3;
		}
	}

	//----  COMPUTE METHODS

	ZSPACE_INLINE void computeFaceVolumes(zComputeMesh& inMesh, zPointArray& fCenters, zDoubleArray& fVolumes)
	{
		zUtilsCore core;

		fVolumes.clear();
		fVolumes.assign(inMesh.polygons.size(), float());

		fCenters.clear();
		fCenters.assign(inMesh.polygons.size(), zPoint());

		for (int i = 0; i < inMesh.polygons.size(); i++)
		{
			// compute face center
			for (int j = 0; j < inMesh.polygons[i].size(); j++) fCenters[i] += inMesh.vertexPositions[inMesh.polygons[i][j]];
			fCenters[i] /= inMesh.polygons[i].size();

			int nTris = floor(inMesh.triangles[i].size() / 3);

			// add volume of face tris	
			for (int j = 0; j < inMesh.triangles[i].size(); j += 3)
			{
				double vol = core.getSignedTriangleVolume(inMesh.vertexPositions[inMesh.triangles[i][j + 0]], inMesh.vertexPositions[inMesh.triangles[i][j + 1]], inMesh.vertexPositions[inMesh.triangles[i][j + 2]]);
				fVolumes[i] += vol;
			}

			// add volumes of tris formed by each pair of face edge vertices and face center
			for (int j = 0; j < inMesh.polygons[i].size(); j++)
			{
				int prevId = (j - 1 + inMesh.polygons[i].size()) % inMesh.polygons[i].size();

				double vol = core.getSignedTriangleVolume(inMesh.vertexPositions[inMesh.polygons[i][j]], inMesh.vertexPositions[inMesh.polygons[i][prevId]], fCenters[i]);
				fVolumes[i] += vol;
			}

			//absolute value
			fVolumes[i] = abs(fVolumes[i]);
		}
	}

	ZSPACE_INLINE void computeFaceNormals(zComputeMesh& inMesh, zVectorArray& fNormals)
	{
		if (inMesh.triangles.size() == 0) throw std::invalid_argument(" error: triangles container is empty.");

		fNormals.clear();
		fNormals.assign(inMesh.polygons.size(), zVector());

		for (int i = 0; i < inMesh.polygons.size(); i++)
		{
			int nTris = floor(inMesh.triangles[i].size() / 3);

			for (int j = 0; j < inMesh.triangles[i].size(); j += 3)
			{
				zVector cross = (inMesh.vertexPositions[inMesh.triangles[i][j + 1]] - inMesh.vertexPositions[inMesh.triangles[i][j + 0]]) ^ (inMesh.vertexPositions[inMesh.triangles[i][j + 2]] - inMesh.vertexPositions[inMesh.triangles[i][j + 0]]);
				cross.normalize();

				fNormals[i] += cross;
			}

			fNormals[i] /= nTris;
			fNormals[i].normalize();
		}
	}

	ZSPACE_INLINE void computeQuadPlanarityDeviation(zComputeMesh& inMesh, zDoubleArray& fDeviations)
	{
		if (inMesh.vertexPositions.size() == 0) throw std::invalid_argument(" error: compute mesh containers are empty.");

		zUtilsCore core;

		fDeviations.clear();
		fDeviations.assign(inMesh.polygons.size(), double());

		for (int j = 0; j < inMesh.polygons.size(); j++)
		{
			double uA, uB;
			zPoint pA, pB;

			bool check = core.line_lineClosestPoints(inMesh.vertexPositions[inMesh.polygons[j][0]], inMesh.vertexPositions[inMesh.polygons[j][2]], inMesh.vertexPositions[inMesh.polygons[j][1]], inMesh.vertexPositions[inMesh.polygons[j][3]], uA, uB, pA, pB);

			fDeviations[j] = pA.distanceTo(pB);
		}
		
	}

	ZSPACE_INLINE void computeGaussianCurvatures(zComputeMesh& inMesh, zInt2DArray& cVertices, zBoolArray& vBoundary, zDoubleArray& vGaussianCurvatures)
	{
		if (inMesh.vertexPositions.size() == 0) throw std::invalid_argument(" error: compute mesh containers are empty.");
		if (inMesh.vertexPositions.size() != vBoundary.size()) throw std::invalid_argument(" error: input containers are not the same size.");
		if (inMesh.vertexPositions.size() != cVertices.size()) throw std::invalid_argument(" error: input containers are not the same size.");

		vGaussianCurvatures.clear();
		vGaussianCurvatures.assign(inMesh.vertexPositions.size(), double());

		for (int i = 0; i < inMesh.vertexPositions.size(); i++)
		{
			double angleSum = 0;
			double cotangentSum = 0;
			double areaSum = 0;
			double areaSumMixed = 0;
			double edgeLengthSquare = 0;
			float gaussianCurv = 0;
			float gaussianAngle = 0;

			if (!vBoundary[i])
			{
				zPoint pt = inMesh.vertexPositions[i];
				float multFactor = 0.125;

				int i = 0;
				for (auto v : cVertices[i])
				{
					int next = (i + 1) % cVertices[i].size();
					int prev = (i + cVertices[i].size() - 1) % cVertices[i].size();

					zVector pt1 = inMesh.vertexPositions[v];
					zVector pt2 = inMesh.vertexPositions[cVertices[i][next]];
					zVector pt3 = inMesh.vertexPositions[cVertices[i][prev]];

					zVector p01 = pt - pt1;
					zVector p02 = pt - pt2;
					zVector p10 = pt1 - pt;
					zVector p20 = pt2 - pt;
					zVector p12 = pt1 - pt2;
					zVector p21 = pt2 - pt1;
					zVector p31 = pt3 - pt1;

					zVector cr = (p10) ^ (p20);

					float ang = (p10).angle(p20);
					angleSum += ang;
					cotangentSum += (((p20) * (p10)) / cr.length());


					float e_Length = (pt1 - pt2).length();

					edgeLengthSquare += (e_Length * e_Length);

					zVector cr_alpha = (p01) ^ (p21);
					zVector cr_beta = (p01) ^ (p31);

					float coTan_alpha = (((p01) * (p21)) / cr_alpha.length());
					float coTan_beta = (((p01) * (p31)) / cr_beta.length());

					// check if triangle is obtuse
					if ((p10).angle(p20) <= 90 && (p01).angle(p21) <= 90 && (p12).angle(p02) <= 90)
					{
						areaSumMixed += (coTan_alpha + coTan_beta) * edgeLengthSquare * 0.125;
					}
					else
					{
						double triArea = (((p10) ^ (p20)).length()) / 2;

						if ((ang) <= 90) areaSumMixed += triArea * 0.25;
						else areaSumMixed += triArea * 0.5;

					}

					i++;
				}

				gaussianCurv = (360 - angleSum) / ((0.5 * areaSum) - (multFactor * cotangentSum * edgeLengthSquare));

			}

			vGaussianCurvatures[i] = gaussianCurv;
			
		}

	}


}