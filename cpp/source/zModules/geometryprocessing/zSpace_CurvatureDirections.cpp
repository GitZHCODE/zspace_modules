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

#include <headers/zModules/geometryprocessing/zSpace_CurvatureDirections.h>

namespace  zSpace
{

	//---- EXTERNAL METHODS 

	ZSPACE_INLINE void trimesh_curvaturedirections(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, double* outK1Directions, double* outK2Directions)
	{
		Eigen::MatrixXd V (numVerts, 3);
		Eigen::MatrixXi F (numFaces, 3);

		// fill vertex matrix
		for (int i = 0; i < numVerts; i++)
		{
			V(i, 0) = _vertexPositions[i * 3 + 0];
			V(i, 1) = _vertexPositions[i * 3 + 1];
			V(i, 2) = _vertexPositions[i * 3 + 2];
		}

		// fill face matrix
		for (int i = 0; i < numFaces; i++)
		{
			F(i, 0) = _polyConnects[i * 3 + 0];
			F(i, 1) = _polyConnects[i * 3 + 1];
			F(i, 2) = _polyConnects[i * 3 + 2];
		}

		// Compute curvature directions via quadric fitting
		MatrixXd PD1, PD2;
		VectorXd PV1, PV2;
		igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
			

		for (int i = 0; i < numVerts; i++)
		{
			outK1Directions[i * 3 + 0] = PD1(i, 0);
			outK1Directions[i * 3 + 1] = PD1(i, 1);
			outK1Directions[i * 3 + 2] = PD1(i, 2);

			outK2Directions[i * 3 + 0] = PD2(i, 0);
			outK2Directions[i * 3 + 1] = PD2(i, 1);
			outK2Directions[i * 3 + 2] = PD2(i, 2);
		}
		
	}

}

