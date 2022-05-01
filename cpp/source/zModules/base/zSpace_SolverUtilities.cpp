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

#include <headers/zModules/base/zSpace_SolverUtilities.h>


namespace  zSpace
{
	//----  SET EXTERN VARIABLES

	zBoolArray vFixed = vector<bool>();

	float dT = 0.5;

	zIntergrationType intType = zRK4;

	//----  EXTERN METHODS

	ZSPACE_MODULES_INLINE void computeMesh_setFixed(int* _fixedVertices, int numFixed)
	{
		if (numFixed >= compMesh.nV) throw std::invalid_argument(" error: number of fixed vertices greater than number of vertices.");

		// set fixed
		if (_fixedVertices)
		{
			for (int i = 0; i < numFixed; i++)
			{
				compMesh.fnParticles[_fixedVertices[i]].setFixed(true);
			}
		}

		vFixed.clear();
		vFixed.assign(compMesh.nV, false);		

		for (int i = 0; i < numFixed; i++)
		{
			vFixed[_fixedVertices[i]] = true;
		}
	}

	ZSPACE_MODULES_INLINE void computeMesh_setDT(double _dT)
	{
		dT = _dT;
	}

	ZSPACE_MODULES_INLINE void computeMesh_setMass(int* _vMass, int numVerts)
	{
		if (numVerts != compMesh.vertexPositions.size()) throw std::invalid_argument(" error: number of fixed vertices greater than number of vertices.");

		// set mass
		if (_vMass)
		{
			for (int i = 0; i < numVerts; i++)
			{
				compMesh.fnParticles[i].setMass(_vMass[i]);
			}
		}
	}
}