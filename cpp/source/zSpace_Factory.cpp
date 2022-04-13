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

#include <headers/zSpace_Factory.h>


namespace  zSpace
{
	
	
	ZSPACE_MODULES void ConstructTopology(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, zObjMesh& out_mesh)
	{

		zPointArray vPos;
		zIntArray pConnects;
		zIntArray pCount;

		for (int i = 0; i < numVerts; i+=3)
		{
			zVector v;
			v = zVector(_vertexPositions[i], _vertexPositions[i + 1], _vertexPositions[i + 3]) ;
			vPos.push_back(v);
		}


		int polyconnectsCurrentIndex = 0;
		for (int i = 0; i < numFaces; i++)
		{
			int num_faceVerts = _polyCounts[i];
			pCount.push_back(num_faceVerts);

			for (int j = 0; j < num_faceVerts; j++)
			{
				pConnects.push_back(_polyConnects[polyconnectsCurrentIndex + j]);
			}

			polyconnectsCurrentIndex += num_faceVerts;
		}


		zFnMesh fnMesh(out_mesh);
		fnMesh.create(vPos, pConnects, pCount);
		
	}

}