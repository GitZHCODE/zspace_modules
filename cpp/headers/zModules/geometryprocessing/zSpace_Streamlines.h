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

#ifndef ZSPACE_MODULES_STREAMLINES_H
#define ZSPACE_MODULES_STREAMLINES_H

#pragma once
#include <vector>
#include <string>
using namespace std;

#include <headers/zModules/base/zSpace_MeshUtilities.h>

#include <headers/zModules/base/zSpace_SolverUtilities.h>


namespace  zSpace
{
	//--------------------------
	//---- GLOBAL VARIABLES
	//--------------------------
	
	/*!<\brief per vertex vectors of the field. The length should be same as the input mesh vertices. */
	zVectorArray vFieldVectors;

	/*!<\brief streamType - zForwardbackward / zForward/ zBackward.*/
	zFieldStreamType streamType;

	/*!<\brief seperation distance between stream lines.*/
	double dSep;

	/*!<\brief test seperation distance between stream lines.*/
	double dTest;

	/*!<\brief minimum length of stream.*/
	double minLength;

	/*!<\brief maximum length of stream.*/
	double maxLength;


	/*!<\brief angle of stream rotation.*/
	double angle;

	/*!<\brief boolean is true if the backward direction is flipped.*/
	bool flipBackward;

	//--------------------------
	//----  METHODS
	//--------------------------
	
	ZSPACE_MODULES_INLINE void createStreamGraph(zObjGraph& streamGraphObj, zVector& seedPoint);

	ZSPACE_MODULES_INLINE bool getFieldValue(zPoint& samplePos, zVector& fieldValue);

	//--------------------------
	//---- EXTERNAL METHODS
	//--------------------------

	/*! \brief This method initialises the solver for planarisation.
	*
	*	\param	[in]	_vertexPositions		- input container of vertex positions. Collapsed 1D array of size numVerts * 3.
	*	\param	[in]	_polyCounts				- input container of number of vertices per polygon of the mesh.
	*	\param	[in]	_polyConnects			- input container of polygon connectivity. Collapsed 1D array of size numFaces * (numVerts per face).
	* 	\param	[in]	_triCounts				- input container of number of triangles per polygon of the mesh.
	*	\param	[in]	_triConnects			- input container of triangle connectivity. Collapsed 1D array of size numFaces * (numtriangles per face * 3).
	* 	\param	[in]	numVerts				- input number of vertices in the mesh.
	*  	\param	[in]	numFaces				- input number of faces/polygons in the mesh.
	*  	\param	[out]	outPlanarityDeviations	- output container of planarity deviations per face/polygon.
	*  	\param	[out]	outGaussianCurvatures	- output container of gaussian curvatures per vertex.
	*	\return			bool					- output boolean - true if setup is successful.
	*	\since version 0.0.4
	*/
	extern "C" ZSPACE_MODULES int streamlines_initialise(double* _vertexPositions, int* _polyCounts, int* _polyConnects, double* _vField, int numVerts, int numFaces);


	/*! \brief This method creates the stream lines and stores them as a graph.
	*
	*	\param	[out]	streams							- container of streams.
	*	\param	[in]	start_seedPoints				- container of start seed positions. If empty a random position in the field is considered.
	*	\param	[in]	seedStreamsOnly					- generates streams from the seed points only if true.
	*	\since version 0.0.1
	*/
	extern "C" ZSPACE_MODULES void computeMesh_streamlines( double* start_seedPoints, int numStartSeeds, double _dSep, double _dTest, double _angle, double* outStreamPositions );


}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/geometryprocessing/zSpace_Streamlines.cpp>
#endif

#endif
