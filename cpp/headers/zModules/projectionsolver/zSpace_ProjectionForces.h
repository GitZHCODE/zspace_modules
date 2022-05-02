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

#ifndef ZSPACE_MODULES_PROJECTIONFORCES_H
#define ZSPACE_MODULES_PROJECTIONFORCES_H

#pragma once

#include <headers/zModules/base/zSpace_MeshUtilities.h>
#include <headers/zModules/geometryprocessing/zSpace_Curvatures.h>

namespace  zSpace
{

	//--------------------------
	//---- PROJECTION FORCE METHODS
	//--------------------------

	/*! \brief This method adds the gravitational force to the input mesh.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	*	\param	[in]	gForce					- gravitational force vector.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void addGravityForce(zComputeMesh& inMesh, zVector &gForce);

	/*! \brief This method adds the drag force to the input mesh.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	*	\param	[in]	drag					- drag constant.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void addDragForce(zComputeMesh& inMesh, float drag);

	/*! \brief This method adds the drag force to the input mesh.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	*	\param	[in]	restLength				- input container of restlengths per edge.
	*	\param	[in]	springConstant			- input spring constant.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void addSpringForce(zComputeMesh& inMesh, zFloatArray &restLength , double strength = 1);

	/*! \brief This method adds the drag force to the input mesh.
	*	\details based on https://github.com/Dan-Piker/K2Goals/blob/master/TangentialSmooth.cs
	*	\param	[in]	inMesh					- input compute mesh object.
	*	\param	[in]	restLength				- input container of restlengths per edge.
	*	\param	[in]	springConstant			- input spring constant.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void addSmoothnessForce(zComputeMesh& inMesh, double strength = 1);


	/*! \brief This method adds the planarisation forces to the input mesh.
	*
	*	\param	[in]	inMesh					- input compute mesh object.		
	*	\param	[in]	type					- input planarisation type - zQuadPlanar or zVolumePlanar.
	* 	\param	[in]	tolerance				- input tolerance value belwo which the force isnt applied.
	*  	\param	[out]	exit					- output boolean true if all the planarity deviations are below tolerance.
	*  	\param	[out]	planarityDeviations		- output container of planarity deviations per face.
	*  	\param	[in]	targetCenters			- container of target origin per face. Used only for zVolumePlanar planarisation type.
	*  	\param	[in]	targetNormals			- container of target normals per face. Used only for zVolumePlanar planarisation type.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void addPlanarityForces(zComputeMesh &inMesh, zPlanarSolverType type, zPointArray& targetCenters, zVectorArray& targetNormals, double &tolerance, zDoubleArray &planarityDeviations, bool& exit);
	
	/*! \brief This method adds the planarisation forces to the input mesh.
	*
	*	\param	[in]	inMesh					- input compute mesh object.
	*	\param	[in]	type					- input planarisation type - zQuadPlanar or zVolumePlanar.
	* 	\param	[in]	tolerance				- input tolerance value belwo which the force isnt applied.
	*  	\param	[out]	exit					- output boolean true if all the planarity deviations are below tolerance.
	*  	\param	[out]	planarityDeviations		- output container of planarity deviations per face.
	*  	\param	[in]	targetCenters			- container of target origin per face. Used only for zVolumePlanar planarisation type.
	*  	\param	[in]	targetNormals			- container of target normals per face. Used only for zVolumePlanar planarisation type.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void addGaussianForces(zComputeMesh& inMesh, zInt2DArray& cVertices, zBoolArray& vBoundary, double& tolerance,  zDoubleArray& vGaussianCurvatures, bool& exit);

	/*! \brief This method adds the minimize area forces (minimal surface) to the input mesh.
	*	\details based on http://courses.cms.caltech.edu/cs177/hmw/Hmw2.pdf , https://github.com/Dan-Piker/K2Goals/blob/master/SoapFilm.cs
	*	\param	[in]	inMesh					- input compute mesh object.
	*	\param	[in]	type					- input planarisation type - zQuadPlanar or zVolumePlanar.
	* 	\param	[in]	strength				- input strength of the force.
	*	\since version 0.0.4
	*/
	ZSPACE_MODULES void addMinimizeAreaForces(zComputeMesh& inMesh, double strength = 1);


}

#if defined(ZSPACE_MODULES_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zModules/projectionsolver/zSpace_ProjectionForces.cpp>
#endif

#endif
