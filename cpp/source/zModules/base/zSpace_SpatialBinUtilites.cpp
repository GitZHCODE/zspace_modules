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

#include <headers/zModules/base/zSpace_SpatialBinUtilities.h>

namespace  zSpace
{
	//----  SET EXTERN VARIABLES

	zPointArray bin_positions = zPointArray();

	zInt2DArray vin_vIds = zInt2DArray();

	int bin_RES = 100;

	zDomainVector bin_bounds = zDomainVector();

	zVector bin_unitVec = zVector();

	//----  METHODS 

	ZSPACE_MODULES_INLINE void spatialBin(zComputeMesh& inMesh, int _res)
	{
		bin_RES = _res;
		vin_vIds.clear();
		vin_vIds.assign(_res * _res * _res, zIntArray());

		bin_positions.clear();
		
		//compute bounds
		computeBounds(inMesh, bin_bounds);

		float unit_X = (bin_bounds.max.x - bin_bounds.min.x) / bin_RES;
		float unit_Y = (bin_bounds.max.y - bin_bounds.min.y) / bin_RES;
		float unit_Z = (bin_bounds.max.y - bin_bounds.min.y) / bin_RES;

		bin_unitVec = zVector(unit_X, unit_Y, unit_Z);
		zVector startPt = bin_bounds.min;

		// compute bin positions
		for (int i = 0; i < bin_RES; i++)
		{
			for (int j = 0; j < bin_RES; j++)
			{
				for (int k = 0; k < bin_RES; k++)
				{
					zVector pos;
					pos.x = startPt.x + i * bin_unitVec.x;
					pos.y = startPt.y + j * bin_unitVec.y;
					pos.z = startPt.z + k * bin_unitVec.z;

					bin_positions.push_back(pos);

				}
			}
		}

		// bin mesh positions
		for (int i =0; i < inMesh.vertexPositions.size(); i++)
		{
			int index= -1;
			bool chk = spatialBin_getID(inMesh.vertexPositions[i], index);

			if (chk)
			{
				vin_vIds[index].push_back(i);
			}
		}
	}

	ZSPACE_MODULES_INLINE bool spatialBin_getID(zPoint& inPoint, int& outID)
	{
		int out = -1;

		int _index_X = floor((inPoint.x - bin_bounds.min.x) / bin_unitVec.x);
		int _index_Y = floor((inPoint.y - bin_bounds.min.y) / bin_unitVec.y);
		int _index_Z = floor((inPoint.z - bin_bounds.min.z) / bin_unitVec.z);

		outID = _index_X * (bin_RES * bin_RES) + (_index_Y * bin_RES) + _index_Z;

		if (outID < 0 && outID >= bin_positions.size()) return false;
		else return true;		
	}

	ZSPACE_MODULES_INLINE void spatialBin_getNeighbourRing(int binID, zIntArray& ringNeighbours)
	{
		int numRings = 1;
		int numBins = bin_positions.size();

		ringNeighbours.clear();

		int idX = floor(binID / (bin_RES * bin_RES));
		int idY = floor((binID - (idX * bin_RES * bin_RES)) / (bin_RES));
		int idZ = binID % bin_RES;

		int startIdX = -numRings;
		if (idX == 0) startIdX = 0;

		int startIdY = -numRings;
		if (idY == 0) startIdY = 0;

		int startIdZ = -numRings;
		if (idZ == 0) startIdZ = 0;

		int endIdX = numRings;
		if (idX == bin_RES - 1) endIdX = 0;

		int endIdY = numRings;
		if (idY == bin_RES - 1) endIdY = 0;

		int endIdZ = numRings;
		if (idZ == bin_RES - 1) endIdZ = 0;


		for (int i = startIdX; i <= endIdX; i++)
		{
			for (int j = startIdY; j <= endIdY; j++)
			{

				for (int k = startIdZ; k <= endIdZ; k++)
				{
					int newId_X = idX + i;
					int newId_Y = idY + j;
					int newId_Z = idZ + k;

					int newId = (newId_X * (bin_RES * bin_RES)) + (newId_Y * bin_RES) + newId_Z;

					if (newId < numBins && newId >= 0) ringNeighbours.push_back(newId);
				}
			}

		}
	}


}

