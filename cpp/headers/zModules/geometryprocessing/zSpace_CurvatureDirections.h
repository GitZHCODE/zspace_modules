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

#ifndef ZSPACE_MODULES_CURVATUREDIRECTIONS_H
#define ZSPACE_MODULES_CURVATUREDIRECTIONS_H

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
#include <igl/read_triangle_mesh.h>


namespace  zSpace
{
	//--------------------------
	//---- GLOBAL VARIABLES
	//--------------------------
	zComputeMesh curvatureMesh;
	
	//--------------------------
	//---- EXTERNAL METHOD
	//--------------------------

	/*! \brief This method initialises the solver for planarisation.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	*  	\param	[out]	outPlanarityDeviations	- output container of planarity deviations per face/polygon.
	*  	\param	[out]	outGaussianCurvatures	- output container of gaussian curvatures per vertex.
	*	\warning		works with triangle meshes.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES void trimesh_curvaturedirections(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, double* outK1Directions, double* outK2Directions);

}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/geometryprocessing/zSpace_CurvatureDirections.cpp>
#endif

#endif
