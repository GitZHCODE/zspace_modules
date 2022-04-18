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


#ifdef ZSPACE_MODULES
#undef ZSPACE_MODULES
#endif

#if defined (ZSPACE_MODULES_DYNAMIC_LIBRARY)
#  define ZSPACE_MODULES_INLINE 
#else
#  define ZSPACE_MODULES_INLINE inline
#endif

#define ZSPACE_MODULES __declspec(dllexport)