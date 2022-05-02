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

#ifndef ZSPACE_MODULES_CURVATURE_H
#define ZSPACE_MODULES_CURVATURE_H

#pragma once
#include <vector>
#include <string>
using namespace std;

#include <headers/zInterface/functionsets/zFnParticle.h>

#include <headers/zModules/base/zSpace_Modules.h>
#include <headers/zModules/base/zSpace_ComputeMesh.h>
#include <headers/zModules/base/zSpace_MeshUtilities.h>

#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/gaussian_curvature.h>
#include <igl/read_triangle_mesh.h>


namespace  zSpace
{
	//--------------------------
	//---- GLOBAL VARIABLES
	//--------------------------
	
	//--------------------------
	//----  METHODS
	//--------------------------
	
	/*! \brief This method computes the principal curvature directions of the input mesh.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	*  	\param	[out]	PV1						- output container of principal curvature value 1. 
	*  	\param	[out]	PV2						- output container of principal curvature value 2. 
	*  	\param	[out]	PD1						- output matrix of principal curvature direction 1. 
	*  	\param	[out]	PD2						- output matrix of principal curvature direction 2. 
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void computeCurvatureDirections(zComputeMesh& inMesh, VectorXd& PV1, VectorXd& PV2, MatrixXd &PD1, MatrixXd &PD2);

	/*! \brief This method computes the gaussian curvature of the input mesh.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	*	\param	[out]	K						- output container of gaussian curvature per vertex.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void computeGaussianCurvature(zComputeMesh &inMesh, VectorXd& K);

	/*! \brief This method computes the pricipal curvature directions of the input triangle mesh.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	*	\param	[out]	HV						- output container of mean curvature per vertex.
	*	\param	[out]	HN						- output matrix of mean curvature normal per vertex.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void computeMeanCurvature(zComputeMesh& inMesh, VectorXd& HV, MatrixXd& HN);

	//--------------------------
	//---- EXTERNAL METHODS
	//--------------------------

	/*! \brief This method computes the principal curvature directions of the global compute mesh.
	*
	*  	\param	[out]	outpV1					- output container of principal curvature value 1. 1D array of size numVerts.
	*  	\param	[out]	outpV2					- output container of principal curvature value 2. 1D array of size numVerts.
	*  	\param	[out]	outpD1					- output container of principal curvature direction 1. Collapsed 1D array of size numVerts * 3.
	*  	\param	[out]	outpD2					- output container of principal curvature direction 2. Collapsed 1D array of size numVerts * 3.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES bool computeMesh_curvatureDirections( double* outpV1, double* outpV2, double* outpD1, double* outpD2);

	/*! \brief This method computes the gaussian curvature of the global compute mesh.
	*
	*  	\param	[out]	outGV					- output container of principal curvature value 1. 1D array of size numVerts.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES bool computeMesh_gaussianCurvature(double* outGV);

	/*! \brief This method computes the mean curvature of the global compute mesh.
	*
	*  	\param	[out]	outHV					- output container of principal curvature value 1. 1D array of size numVerts.
	*  	\param	[out]	outHV					- output container of principal curvature value 1. 1D array of size numVerts.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES bool computeMesh_meanCurvature(double* outHV, double* outHN);

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/geometryprocessing/zSpace_Curvatures.cpp>
#endif

#endif
