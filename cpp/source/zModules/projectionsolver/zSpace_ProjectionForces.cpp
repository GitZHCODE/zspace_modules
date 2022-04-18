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
		int numVerts = inMesh.vertexPositions.size();

		for (int j = 0; j < numVerts; j++)
		{
			inMesh.fnParticles[j].addForce(gForce);
		}
	}

	ZSPACE_MODULES_INLINE void addPlanarityForces(zComputeMesh& inMesh, zPlanarType type, zPointArray& targetCenters, zVectorArray& targetNormals, double& tolerance, zDoubleArray& planarityDeviations, bool& exit)
	{
		zUtilsCore core;
		int numFaces = inMesh.polygons.size();

		if (type == zQuadPlanar)
		{
			for (int j = 0; j < numFaces; j++)
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
			for (int j = 0; j < numFaces; j++)
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

		}
		
	}

	//---- EXTERNAL METHODS FOR CONSTRAINTS

	ZSPACE_MODULES_INLINE void setFixed(zComputeMesh& inMesh, int* _fixedVertices, int numFixed)
	{
		if (numFixed >= inMesh.vertexPositions.size()) throw std::invalid_argument(" error: number of fixed vertices greater than number of vertices.");

		// set fixed
		if (_fixedVertices)
		{
			for (int i = 0; i < numFixed; i++)
			{
				inMesh.fnParticles[_fixedVertices[i]].setFixed(true);
			}
		}
	}

	ZSPACE_MODULES_INLINE void setMass(zComputeMesh& inMesh, int* _vMass, int numVerts)
	{
		if (numVerts != inMesh.vertexPositions.size()) throw std::invalid_argument(" error: number of fixed vertices greater than number of vertices.");

		// set mass
		if (_vMass)
		{
			for (int i = 0; i < numVerts; i++)
			{
				inMesh.fnParticles[i].setMass(_vMass[i]);
			}
		}
	}

}

