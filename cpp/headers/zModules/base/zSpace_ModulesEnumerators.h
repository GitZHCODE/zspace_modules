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


#ifndef ZSPACE_MODULES_ENUMERATORS
#define ZSPACE_MODULES_ENUMERATORS


namespace zSpace
{
	/** \addtogroup zModulesEnumerators
	*	\brief  The enumerators of the modules library.
	*  @{
	*/

	/*! \enum	zPlanarSolverType
	*	\brief	Types of planarisation solver.
	*	\since	version 0.0.4
	*/
	enum zPlanarSolverType { zQuadPlanar = 400, zVolumePlanar };

	/*! \enum	zMSSolverType
	*	\brief	Types of planarisation.
	*	\since	version 0.0.4
	*/
	enum zMSSolverType { zMinimiseArea = 410, zRestlength };

	/** @}*/
}

#endif
