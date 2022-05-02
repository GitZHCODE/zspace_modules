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

#include <headers/zModules/projectionsolver/zSpace_ProjectionForces.h>

namespace  zSpace
{

	//---- PROJECTION FORCE METHODS
	
	ZSPACE_MODULES_INLINE void addGravityForce(zComputeMesh& inMesh, zVector& gForce)
	{
		for (int j = 0; j < inMesh.nV; j++)
		{
			inMesh.fnParticles[j].addForce(gForce);
		}
	}

	ZSPACE_MODULES_INLINE void addDragForce(zComputeMesh& inMesh, float drag)
	{
		for (int j = 0; j < inMesh.nV; j++)
		{
			zVector v = inMesh.fnParticles[j].getVelocity();
			zVector pForce =  v * drag * -1;
			inMesh.fnParticles[j].addForce(pForce);
		}
	}

	ZSPACE_MODULES_INLINE void addSpringForce(zComputeMesh& inMesh, zFloatArray& restLength, double strength)
	{
		for (int j = 0; j < inMesh.nE; j++)
		{
			int v1 = inMesh.edges[j][0];
			int v2 = inMesh.edges[j][1];

			zVector e = inMesh.vertexPositions[v2] - inMesh.vertexPositions[v1];
			float eLen = e.length();
			e.normalize();
					
			float restLen = restLength[j];

			float val = strength * (eLen - restLen);
			zVector pForce_v1 = e * (val * 0.5);				

			zVector pForce_v2 = pForce_v1 * -1;			

			inMesh.fnParticles[v1].addForce(pForce_v1);
			inMesh.fnParticles[v2].addForce(pForce_v2);
		}
	}

	ZSPACE_MODULES_INLINE void addSmoothnessForce(zComputeMesh& inMesh, double strength)
	{
		for (int j = 0; j < inMesh.nV; j++)
		{
			int numCVerts = inMesh.cVertices[j].size();
			
			zPoint Avg;
			zVectorArray Vecs;
			Vecs.assign(numCVerts, zVector()); 

			for (int i = 0; i < numCVerts; i++)
			{
				Avg = Avg + inMesh.vertexPositions[inMesh.cVertices[j][i]];;
				Vecs[i] = inMesh.vertexPositions[inMesh.cVertices[j][i]] - inMesh.vertexPositions[j];
			}

			double Inv = 1.0 / (numCVerts);
			Avg = Avg * Inv;
			zVector Smooth =  (Avg - inMesh.vertexPositions[j]) * 0.5;

			zVector Normal;
			for (int i = 0; i < Vecs.size(); i++)
			{
				Normal += (Vecs[i] ^ Vecs[(i + 1) % Vecs.size()]);
			}
			Normal.normalize();
			Smooth -= Normal * (Normal * Smooth);

			zVector pForce1 = Smooth * strength;
			inMesh.fnParticles[j].addForce(pForce1);
			
			/*Smooth *= -Inv;
			zVector pForce2 = Smooth * smoothnessConstant;
			for (int i = 0; i < numCVerts; i++)
			{
				inMesh.fnParticles[inMesh.cVertices[j][i]].addForce(pForce2);
			}*/
		}

		
	}

	ZSPACE_MODULES_INLINE void addPlanarityForces(zComputeMesh& inMesh, zPlanarSolverType type, zPointArray& targetCenters, zVectorArray& targetNormals, double& tolerance, zDoubleArray& planarityDeviations, bool& exit)
	{
		zUtilsCore core;
		if (type == zQuadPlanar)
		{
			for (int j = 0; j < inMesh.nF; j++)
			{
				double uA, uB;
				zPoint pA, pB;

				bool check = core.line_lineClosestPoints(inMesh.vertexPositions[inMesh.polygons[j][0]], inMesh.vertexPositions[inMesh.polygons[j][2]], inMesh.vertexPositions[inMesh.polygons[j][1]], inMesh.vertexPositions[inMesh.polygons[j][3]], uA, uB, pA, pB);

				planarityDeviations[j] = pA.distanceTo(pB);

				if (planarityDeviations[j] > tolerance)
				{
					exit = false;

					zVector dir = pB - pA;
					dir.normalize();

					zVector pForceA = dir * planarityDeviations[j] * 0.5;
					zVector pForceB = dir * planarityDeviations[j] * -0.5;

					inMesh.fnParticles[inMesh.polygons[j][0]].addForce(pForceA);
					inMesh.fnParticles[inMesh.polygons[j][2]].addForce(pForceA);

					inMesh.fnParticles[inMesh.polygons[j][1]].addForce(pForceB);
					inMesh.fnParticles[inMesh.polygons[j][3]].addForce(pForceB);
				}

			}
		}
		
		else if (type == zVolumePlanar)
		{
			for (int j = 0; j < inMesh.nF; j++)
			{
				if (planarityDeviations[j] > tolerance)
				{
					exit = false;

					for (int k = 0; k < inMesh.polygons[j].size(); k++)
					{
						double dist = core.minDist_Point_Plane(inMesh.vertexPositions[inMesh.polygons[j][k]], targetCenters[j], targetNormals[j]);
						zVector pForce = targetNormals[j] * dist * -1.0;
						inMesh.fnParticles[inMesh.polygons[j][k]].addForce(pForce);

					}
				}
			}
		}
	}

	ZSPACE_MODULES_INLINE void addGaussianForces(zComputeMesh& inMesh, zInt2DArray& cVertices, zBoolArray& vBoundary, double& tolerance, zDoubleArray& vGaussianCurvatures, bool& exit)
	{
		int numVerts = inMesh.vertexPositions.size();

		//compute gaussian curvature
		VectorXd vGauss;
		computeGaussianCurvature(inMesh, vGauss);

		// compute gradient and force
		for (int i = 0; i < inMesh.nV; i++)
		{
			zVector gForce;
			// boundary vertices

			// interval vertices

			inMesh.fnParticles[i].addForce(gForce);
		}
		
	}

	ZSPACE_MODULES_INLINE void addMinimizeAreaForces(zComputeMesh& inMesh, double strength)
	{	
		int currentIndex = 0;
		for (int i = 0; i < inMesh.nF; i++)
		{
			int numTrisVerts = inMesh.triangles[i].size();

			for (int j = 0; j < numTrisVerts; j += 3)
			{
				zPoint PA = inMesh.vertexPositions[inMesh.triangles[i][j + 0]];
				zPoint PB = inMesh.vertexPositions[inMesh.triangles[i][j + 1]];
				zPoint PC = inMesh.vertexPositions[inMesh.triangles[i][j + 2]];

				zVector AB = PB - PA;
				zVector BC = PC - PB;
				zVector CA = PA - PC;

				zVector Normal = AB ^ BC; 
				Normal.normalize();

				zVector V0 = (BC ^ Normal) * 0.5;
				zVector V1 = (CA ^ Normal) * 0.5;  
				
				zVector pForce0 = V0 * strength;
				inMesh.fnParticles[inMesh.triangles[i][j + 0]].addForce(pForce0);

				zVector pForce1 = V1 * strength;
				inMesh.fnParticles[inMesh.triangles[i][j + 1]].addForce(pForce1);

				zVector pForce2 = ((V0 * -1) - V1) * strength;
				inMesh.fnParticles[inMesh.triangles[i][j + 2]].addForce(pForce2);				
			}
		}

		
		
	}

}

