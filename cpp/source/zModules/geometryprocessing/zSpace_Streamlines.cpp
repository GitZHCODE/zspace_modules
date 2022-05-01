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

#include <headers/zModules/geometryprocessing/zSpace_Streamlines.h>

namespace  zSpace
{
	//----  METHODS 

	ZSPACE_MODULES_INLINE void createStreamGraph(zObjGraph& streamGraphObj, zVector& seedPoint)
	{
		vector<zVector> positions;
		vector<int> edgeConnects;


		// move forward
		if (streamType == zForward || streamType == zForwardBackward)
		{
			bool exit = false;

			zVector startForward = seedPoint;

			zObjParticle p;
			p.particle = zParticle(startForward);

			zFnParticle seedForward(p);

			double currentLength = 0.0;

			while (!exit)
			{
				bool firstVertex = (startForward == seedPoint) ? true : false;

				zVector curPos = seedForward.getPosition();

				if (firstVertex)
				{

					positions.push_back(curPos);
				}

				// get field focrce
				zVector fieldForce;
				bool checkBounds;// = fnField.getFieldValue(curPos, zFieldNeighbourWeighted, fieldForce);



				if (!checkBounds)
				{
					exit = true;
					continue;
				}

				// local minima or maxima point
				if (fieldForce.length() == 0)
				{
					exit = true;
					continue;
				}

				// update particle force
				fieldForce.normalize();
				fieldForce *= (dSep * 1.0);

				zVector axis(0, 0, 1);

				double rotateAngle = angle;
				fieldForce = fieldForce.rotateAboutAxis(axis, rotateAngle);

				seedForward.addForce(fieldForce);

				// update seed particle
				//seedForward.integrateForces(*dT, integrationType);
				seedForward.updateParticle(true);
				zVector newPos = seedForward.getPosition();
				//checkBounds = checkFieldBounds(newPos);

				if (!checkBounds)
				{
					//printf("\n bounds 2 working!");

					//printf("\n %1.2f %1.2f %1.2f ", newPos.x, newPos.y, newPos.z);

					exit = true;
					continue;



				}

				int index = -1;
				bool checkRepeat;// = coreUtils.checkRepeatElement(newPos, positions, index);
				if (checkRepeat)
				{
					//printf("\n repeat working!");

					exit = true;
					continue;



				}


				bool validStreamPoint;// = checkValidStreamPosition(newPos, dTest);

				if (!validStreamPoint)
				{
					exit = true;
				}

				// check length
				if (currentLength + curPos.distanceTo(newPos) > maxLength)
				{
					exit = true;
				}

				// add new stream point
				if (!exit)
				{
					if (positions.size() > 0)
					{
						edgeConnects.push_back(positions.size());
						edgeConnects.push_back(positions.size() - 1);
					}

					positions.push_back(newPos);

					currentLength += curPos.distanceTo(newPos);

				}
			}
		}
	}

	ZSPACE_MODULES_INLINE bool getFieldValue(zPoint& samplePos, zVector& fieldValue)
	{
		return ZSPACE_MODULES_INLINE bool();
	}

	//---- EXTERNAL METHODS 

	ZSPACE_MODULES_INLINE int streamlines_initialise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, double* _vField, int numVerts, int numFaces)
	{
		bool out = false;
		if (_vertexPositions && _polyCounts && _polyConnects && _vField)
		{
			// create mesh obj
			createMeshOBJ(_vertexPositions, _polyCounts, _polyConnects, numVerts, numFaces, o_Mesh);

			// set vector field per vertex
			vFieldVectors.clear();
			vFieldVectors.assign(numVerts, zVector());

			for (int i = 0; i < numVerts; i += 1)
			{
				vFieldVectors[i].x = _vField[i * 3 + 0];
				vFieldVectors[i].y = _vField[i * 3 + 1];
				vFieldVectors[i].z = _vField[i * 3 + 2];
			}

			out = true;
		}

		return out;
	}

	ZSPACE_MODULES_INLINE void computeMesh_streamlines(double* _start_seedPoints, int numStartSeeds, double _dSep, double _dTest, double _angle, double* outStreamPositions)
	{
		zUtilsCore core;

		zFnMesh fnMesh(o_Mesh);

		zObjGraphArray streams;
		streams.clear();
		streams.assign(100, zObjGraph());

		zInt2DArray childgraphs;
		zIntArray parentGraph;
		bool alternate = false;

		zPointArray seedPoints;
		zBoolArray validStream;

		int streamCounter = 0;

		// make first stream line
		if (streamCounter == 0)
		{
			bool noStartSeedPoints = false;

			if (numStartSeeds == 0)
			{
				noStartSeedPoints = true;

				zVector minBB, maxBB;
				fnMesh.getBounds(minBB, maxBB);

				zPoint seedPoint = zPoint(core.randomNumber_double(minBB.x, maxBB.x), core.randomNumber_double(minBB.y, maxBB.y), 0);
				seedPoints.push_back(seedPoint);
			}
			else
			{
				seedPoints.assign(numStartSeeds, zPoint());

				for (int i = 0; i < numStartSeeds; i += 1)
				{
					seedPoints[i].x = _start_seedPoints[i * 3 + 0];
					seedPoints[i].y = _start_seedPoints[i * 3 + 1];
					seedPoints[i].z = _start_seedPoints[i * 3 + 2];
				}
			}

			for (int i = 0; i < seedPoints.size(); i++)
			{
				//streams.push_back(zStreamLine());

				bool chk = false;// = createStreamGraph_Influence(streams[streamCounter].graphObj, start_seedPoints[i], fnInfluenceField, min_Power, max_Power);

				if (chk)
				{
					//streams.push_back(temp);

					//streams[streamCounter].isValid = true;
					//streamCounter++;

					alternate = !alternate;

				}
			}

		}

	}

}