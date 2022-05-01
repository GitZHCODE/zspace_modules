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

#include <headers/zModules/statics/zSpace_FDM.h>

namespace  zSpace
{

	//---- METHODS

	ZSPACE_MODULES_INLINE zSparseMatrix getEdgeNodeMatrix(int numRows)
	{
		int n_v = compMesh.nV;

		zSparseMatrix out(numRows, n_v);
		out.setZero();

		int coefsCounter = 0;

		vector<zTriplet> coefs; // 1 for from vertex and -1 for to vertex

		for (auto e : compMesh.edges)
		{
			int v1 = e[0];
			int v2 = e[1];


			if (!vFixed[v1] || !vFixed[v2])
			{
				coefs.push_back(zTriplet(coefsCounter, v1, 1));
				coefs.push_back(zTriplet(coefsCounter, v2, -1));

				coefsCounter++;
			}

		}

		out.setFromTriplets(coefs.begin(), coefs.end());

		return out;
	}

	ZSPACE_MODULES_INLINE zSparseMatrix subMatrix(zSparseMatrix& C, vector<int>& nodes)
	{
		zSparseMatrix C_sub(C.rows(), nodes.size());
		for (int i = 0; i < nodes.size(); i++)C_sub.col(i) = C.col(nodes[i]);
		return C_sub;
	}

	ZSPACE_MODULES_INLINE MatrixXd subMatrix(MatrixXd& X, vector<int>& nodes)
	{
		MatrixXd X_sub(nodes.size(), X.cols());
		for (int i = 0; i < nodes.size(); i++)X_sub.row(i) = X.row(nodes[i]);
		return X_sub;
	}

	//---- EXTERNAL METHODS FOR FDM

	ZSPACE_MODULES_INLINE int computeMesh_fdm(double* vForceDensities, double* vMass,  double* outVertexPositions)
	{

		bool positiveDensities = true;

		int n_v = compMesh.nV;
		int n_edges = compMesh.nE;

		for (auto e: compMesh.edges)
		{		
			int v1 = e[0];
			int v2 = e[1];

			if (vFixed[v1] && vFixed[v2]) n_edges--;
		}

		// POSITION MATRIX
		MatrixXd X = compMesh.V;

		// EDGE NODE MATRIX
		zSparseMatrix C = getEdgeNodeMatrix(n_edges);

		// FORCE DENSITY VECTOR
		VectorXd q(n_edges);
		int FD_EdgesCounter = 0;

		for (int i =0; i< compMesh.edges.size(); i++)
		{
			int v1 = compMesh.edges[i][0];
			int v2 = compMesh.edges[i][1];

			if (!vFixed[v1] || !vFixed[v2])
			{
				q[FD_EdgesCounter] = (vForceDensities[v1] + vForceDensities[v2])  * 0.5;
			
				if (q[FD_EdgesCounter] < 0) positiveDensities = false;

				FD_EdgesCounter++;
			}

		}

		zDiagonalMatrix Q = q.asDiagonal();

		// LOAD VECTOR
		VectorXd p(n_v);

		for (int i = 0; i < n_v; i++)
		{
			p[i] = vMass[i];
		}


		MatrixXd P(n_v, 3);
		P.setConstant(0.0);
		P.col(2) = p.col(0);

		// SUB MATRICES
		vector<int> freeVertices;
		vector<int> fixedVertices;

		for (int j = 0; j < vFixed.size(); j++)
		{
			if (!vFixed[j]) freeVertices.push_back(j);
			else fixedVertices.push_back(j);
		}

		zSparseMatrix Cn = subMatrix(C, freeVertices);
		zSparseMatrix Cf = subMatrix(C, fixedVertices);
		MatrixXd Xf = subMatrix(X, fixedVertices);
		MatrixXd Pn = subMatrix(P, freeVertices);

		zSparseMatrix Cn_transpose;
		Cn_transpose = Cn.transpose();

		//CHOLESKY DECOMPOSITION

		zSparseMatrix Dn = Cn_transpose * Q * Cn;
		zSparseMatrix Df = Cn_transpose * Q * Cf;

		MatrixXd  B = Pn - Df * Xf;

		// solve
		MatrixXd Xn;

		if (positiveDensities)
		{
			SimplicialLLT< zSparseMatrix > solver; // sparse cholesky solver
			solver.compute(Dn); // compute cholesky factors

			if (solver.info() != Eigen::Success)
				return false;


			Xn = solver.solve(B); // solve AX = B ;
			if (solver.info() != Eigen::Success)
				return false;

			//cout << endl << Xn;
		}
		else
		{
			MatrixXd denseDn;
			denseDn = MatrixXd(Dn);
			Xn = denseDn.ldlt().solve(B);

			// convergence error check.
			double relative_error = (denseDn * Xn - B).norm() / B.norm(); // norm() is L2 norm
			cout << endl << relative_error << " FDM - negative" << endl;

		}

		// POSITIONS OF NON FIXED VERTICES
		for (int i = 0; i < freeVertices.size(); i++)
		{
			int id = freeVertices[i];			

			compMesh.vertexPositions[id].x = Xn(i, 0);
			compMesh.vertexPositions[id].y = Xn(i, 1);
			compMesh.vertexPositions[id].z = Xn(i, 2);
		}
		
			
		// output
		for (int i = 0; i < compMesh.vertexPositions.size(); i++)
		{
			outVertexPositions[i * 3 + 0] = compMesh.vertexPositions[i].x;
			outVertexPositions[i * 3 + 1] = compMesh.vertexPositions[i].y;
			outVertexPositions[i * 3 + 2] = compMesh.vertexPositions[i].z;
		}
					
		return true;
	}



}

