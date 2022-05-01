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

	//----  SET EXTERN VARIABLES
	
	zComputeMesh compMesh;

	//----  CREATE METHODS
	
	ZSPACE_MODULES_INLINE void createMeshOBJ(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, zObjMesh& out_mesh)
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

	ZSPACE_MODULES_INLINE void createComputeMesh(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, bool makeDynamic, zComputeMesh& outMesh)
	{
		if (!_vertexPositions || !_polyCounts || !_polyConnects) throw std::invalid_argument(" error: mesh container is empty.");

		outMesh.vertexPositions.clear();
		outMesh.vertexPositions.assign(numVerts, zPoint());

		outMesh.nV = numVerts;

		for (int i = 0; i < numVerts; i += 1)
		{
			outMesh.vertexPositions[i].x = _vertexPositions[i * 3 + 0];
			outMesh.vertexPositions[i].y = _vertexPositions[i * 3 + 1];
			outMesh.vertexPositions[i].z = _vertexPositions[i * 3 + 2];

		}

		outMesh.polygons.clear();
		outMesh.polygons.assign(numFaces, vector<int>());

		outMesh.nF = numFaces;		

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

		//edges

		outMesh.edges.clear();
		outMesh.nE = 0;

		//	halfedge connectivity map
		unordered_map <string, int> existingHalfEdges;

		for (int i = 0; i < outMesh.polygons.size(); i++)
		{
			for (int j = 0; j < outMesh.polygons[i].size(); j++)
			{
				int v1 = outMesh.polygons[i][j];
				
				int next = (j + 1) % outMesh.polygons[i].size();
				int v2 = outMesh.polygons[i][next];

				string e1 = (to_string(v1) + "," + to_string(v2));
				std::unordered_map<std::string, int>::const_iterator got = existingHalfEdges.find(e1);

				if (got != existingHalfEdges.end())
				{
					// edge exists continue
				}
				else
				{
					// add edge
					outMesh.edges.push_back(zIntArray());
					outMesh.edges[outMesh.nE].push_back(v1);
					outMesh.edges[outMesh.nE].push_back(v2);
					
					// add halfedges to map
					existingHalfEdges[e1] = outMesh.nE * 2  + 0;

					string e2 = (to_string(v2) + "," + to_string(v1));
					existingHalfEdges[e2] = outMesh.nE * 2 + 1;

					outMesh.nE++;

				}
			}
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
		outMesh.nT = 0;

		//printf("\n compute Mesh : %i %i %i", outMesh.nV, outMesh.nE, outMesh.nF);
	}

	//----  SET METHODS

	ZSPACE_MODULES_INLINE void setTriangles(zComputeMesh& inMesh, int numFaces, int* _triCounts, int* _triConnects)
	{
		if (numFaces != inMesh.polygons.size()) throw std::invalid_argument(" error: triangles container is empty.");
		if(!_triCounts || !_triConnects) throw std::invalid_argument(" error: triangles container is empty.");

		inMesh.triangles.clear();
		inMesh.triangles.assign(numFaces, vector<int>());

		inMesh.nT = 0;

		int triconnectsCurrentIndex = 0;
		for (int i = 0; i < numFaces; i++)
		{
			inMesh.nT += _triCounts[i];
			int num_triVerts = _triCounts[i];

			for (int j = 0; j < num_triVerts * 3; j++)
			{
				inMesh.triangles[i].push_back(_triConnects[triconnectsCurrentIndex + j]);
			}

			triconnectsCurrentIndex += num_triVerts * 3;
		}
	}

	//----  UPDATE METHODS

	ZSPACE_MODULES_INLINE void updateMatrixV(zComputeMesh& inMesh)
	{
		MatrixXd V(inMesh.nV, 3);
		
		// fill vertex matrix
		for (int i = 0; i < inMesh.nV; i++)
		{
			V(i, 0) = inMesh.vertexPositions[i].x;
			V(i, 1) = inMesh.vertexPositions[i].y;
			V(i, 2) = inMesh.vertexPositions[i].z;
		}

		inMesh.V = V;
	}

	ZSPACE_MODULES_INLINE void updateMatrixFTris(zComputeMesh& inMesh)
	{
		MatrixXi FTris(inMesh.nT, 3);

		int nTris = 0;
		for (int i = 0; i < inMesh.triangles.size(); i++)
		{
			for (int j = 0; j < inMesh.triangles[i].size(); j += 3)
			{
				FTris(nTris, 0) = inMesh.triangles[i][j + 0];
				FTris(nTris, 1) = inMesh.triangles[i][j + 1];
				FTris(nTris, 2) = inMesh.triangles[i][j + 2];

				nTris++;
			}
		}

		inMesh.FTris = FTris;
	}

	//----  COMPUTE METHODS

	ZSPACE_MODULES_INLINE void computeFaceVolumes(zComputeMesh& inMesh, zPointArray& fCenters, zDoubleArray& fVolumes)
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

	ZSPACE_MODULES_INLINE void computeFaceNormals(zComputeMesh& inMesh, zVectorArray& fNormals)
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

	ZSPACE_MODULES_INLINE void computeQuadPlanarityDeviation(zComputeMesh& inMesh, zDoubleArray& fDeviations)
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

	ZSPACE_MODULES_INLINE void computeBounds(zComputeMesh& inMesh, zDomainVector& bounds)
	{
		bounds.min.x = bounds.min.y = bounds.min.z = std::numeric_limits<float>::max();
		bounds.max.x = bounds.max.y = bounds.max.z = std::numeric_limits<float>::min();

		for (auto& p : inMesh.vertexPositions)
		{
			if (p.x < bounds.min.x) bounds.min.x = p.x;
			if (p.y < bounds.min.y) bounds.min.y = p.y;
			if (p.z < bounds.min.z) bounds.min.z = p.z;

			if (p.x > bounds.max.x) bounds.max.x = p.x;
			if (p.y > bounds.max.y) bounds.max.y = p.y;
			if (p.z > bounds.max.z) bounds.max.z = p.z;
		}
	}

	//----  EXTERN METHODS
	ZSPACE_MODULES_INLINE int computeMesh_initialise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces)
	{
		bool out = false;

		if (_vertexPositions && _polyCounts && _polyConnects)
		{
			createComputeMesh(_vertexPositions, _polyCounts, _polyConnects, numVerts, numFaces, true, compMesh);
			updateMatrixV(compMesh);

			// set triangles and matrices for igl method call
			if (_triCounts && _triConnects)
			{
				setTriangles(compMesh, numFaces, _triCounts, _triConnects);
				updateMatrixFTris(compMesh);				
			}

			out = true;
		}

		return out;
	}

	ZSPACE_MODULES_INLINE int heMesh_initialise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces)
	{
		bool out = false;
		if (_vertexPositions && _polyCounts && _polyConnects )
		{
			// create mesh obj
			createMeshOBJ(_vertexPositions, _polyCounts, _polyConnects, numVerts, numFaces, o_Mesh);

			zFnMesh fnMesh(o_Mesh);
			if (fnMesh.numVertices() > 0) out = true;
		}

		return out;
	}

}