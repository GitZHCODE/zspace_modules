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

	/*! \brief This method computes the principal curvature directions of the input triangle mesh.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	_triCounts				- input container of number of triangles per polygon of the mesh.
	*	\param	[in]	_triConnects			- input container of triangle connectivity. Collapsed 1D array of size numFaces * (numtriangles per face * 3).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	*  	\param	[out]	outpV1					- output container of principal curvature value 1. 1D array of size numVerts.
	*  	\param	[out]	outpV2					- output container of principal curvature value 2. 1D array of size numVerts.
	*  	\param	[out]	outpD1					- output container of principal curvature direction 1. Collapsed 1D array of size numVerts * 3.
	*  	\param	[out]	outpD2					- output container of principal curvature direction 2. Collapsed 1D array of size numVerts * 3.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES int curvatureDirections(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, double* outpV1, double* outpV2, double* outpD1, double* outpD2);

	/*! \brief This method computes the gaussian curvature of the input mesh.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	_triCounts				- input container of number of triangles per polygon of the mesh.
	*	\param	[in]	_triConnects			- input container of triangle connectivity. Collapsed 1D array of size numFaces * (numtriangles per face * 3).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	*  	\param	[out]	outGV					- output container of principal curvature value 1. 1D array of size numVerts.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES int gaussianCurvature(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, double* outGV);

	/*! \brief This method computes the mean curvature of the input mesh.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	_triCounts				- input container of number of triangles per polygon of the mesh.
	*	\param	[in]	_triConnects			- input container of triangle connectivity. Collapsed 1D array of size numFaces * (numtriangles per face * 3).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	*  	\param	[out]	outHV					- output container of principal curvature value 1. 1D array of size numVerts.
	*  	\param	[out]	outHV					- output container of principal curvature value 1. 1D array of size numVerts.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES int meanCurvature(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int* _triCounts, int* _triConnects, int numVerts, int numFaces, double* outHV, double* outHN);

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/geometryprocessing/zSpace_Curvatures.cpp>
#endif

#endif
