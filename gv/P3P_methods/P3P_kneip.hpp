/******************************************************************************
* This programe is implemented in Visual Studio 2013.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/

#ifndef GV_P3P_METHODS_MAIN_HPP_
#define GV_P3P_METHODS_MAIN_HPP_

#include <stdlib.h>
#include <gv/types.hpp>

namespace gv
{
	namespace P3P_methods
	{
		//Kneip's method;
		gv::transformations_t p3p_kneip(const bearingVectors_t & f, const points_t & p);

	}
}

#endif /* GV_P3P_METHODS_MAIN_HPP_ */
