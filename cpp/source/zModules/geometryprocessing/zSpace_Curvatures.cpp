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

#include <headers/zModules/geometryprocessing/zSpace_Curvatures.h>

namespace  zSpace
{
	//----  METHODS 

	ZSPACE_MODULES_INLINE void computeCurvatureDirections(zComputeMesh& inMesh, VectorXd& PV1, VectorXd& PV2, MatrixXd& PD1, MatrixXd& PD2)
	{
		igl::principal_curvature(inMesh.V, inMesh.FTris, PD1, PD2, PV1, PV2);
	}

	ZSPACE_MODULES_INLINE void computeGaussianCurvature(zComputeMesh& inMesh, VectorXd& K)
	{
		// Compute integral of Gaussian curvature
		igl::gaussian_curvature(inMesh.V, inMesh.FTris, K);
		// Compute mass matrix
		SparseMatrix<double> M, Minv;
		igl::massmatrix(inMesh.V, inMesh.FTris, igl::MASSMATRIX_TYPE_DEFAULT, M);
		igl::invert_diag(M, Minv);
		// Divide by area to get integral average
		K = (Minv * K).eval();
	}

	ZSPACE_MODULES_INLINE void computeMeanCurvature(zComputeMesh& inMesh, VectorXd& HV, MatrixXd& HN)
	{
		// Alternative discrete mean curvature
		SparseMatrix<double> L, M, Minv;
		igl::cotmatrix(inMesh.V, inMesh.FTris, L);
		igl::massmatrix(inMesh.V, inMesh.FTris, igl::MASSMATRIX_TYPE_VORONOI, M);
		igl::invert_diag(M, Minv);
		
		// Laplace-Beltrami of position
		HN = -Minv * (L * inMesh.V);
		
		// Extract magnitude as mean curvature
		HV = HN.rowwise().norm();
	}

	//---- EXTERNAL METHODS 

	ZSPACE_MODULES_INLINE bool curvatureDirections( double* outpV1, double* outpV2, double* outpD1, double* outpD2)
	{
		bool out = false;

		if ( outpV1 && outpV2 && outpD1 && outpD2)
		{
			// Compute curvature directions via quadric fitting
			MatrixXd PD1, PD2;
			VectorXd PV1, PV2;
			computeCurvatureDirections(compMesh, PV1, PV2, PD1, PD2);
			
			for (int i = 0; i < compMesh.nV; i++)
			{
				outpV1[i] = PV1(i);
				outpV2[i] = PV2(i);

				outpD1[i * 3 + 0] = PD1(i, 0);
				outpD1[i * 3 + 1] = PD1(i, 1);
				outpD1[i * 3 + 2] = PD1(i, 2);

				outpD2[i * 3 + 0] = PD2(i, 0);
				outpD2[i * 3 + 1] = PD2(i, 1);
				outpD2[i * 3 + 2] = PD2(i, 2);
			}

			out = true;
		}
		
		return out;
	}

	ZSPACE_MODULES_INLINE bool gaussianCurvature(double* outGV)
	{
		bool out = false;

		if ( outGV)
		{		

			VectorXd K;
			computeGaussianCurvature(compMesh, K);

			for (int i = 0; i < compMesh.nV; i++)
			{
				outGV[i] = K[i];
			}

			out = true;
		}

		return out;
	}

	ZSPACE_MODULES_INLINE bool meanCurvature( double* outHV, double* outHN)
	{
		bool out = false;

		if (outHV && outHN)
		{
			// Alternative discrete mean curvature
			MatrixXd HN;
			VectorXd HV;

			computeMeanCurvature(compMesh, HV, HN);

			for (int i = 0; i < compMesh.nV; i++)
			{
				outHV[i] = HV(i);

				outHN[i * 3 + 0] = HN(i, 0);
				outHN[i * 3 + 1] = HN(i, 1);
				outHN[i * 3 + 2] = HN(i, 2);
			}

			out = true;
		}

		return out;
	}

}

