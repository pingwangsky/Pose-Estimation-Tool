/******************************************************************************
* This programe is implemented in Visual Studio 2013.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/

#ifndef GV_EXPERIMENT_HELPERS_HPP_
#define GV_EXPERIMENT_HELPERS_HPP_

#include "gv\random_generators.hpp"
#include "gv\types.hpp"
#include <iostream>
#include <iomanip>
#include <stdlib.h>


namespace gv
{
//**************MY programe begin******************
void GeneratePnPdata(
	points_t & worldpoints,
	Image_points & image_points,
	translation_t & position,
	rotation_t & rotation,
	//bearingVectors_t & bearingVectors,
	size_t numberPoints,
	double depth,
	double noise);

void ShowExperimentData(
	points_t & worldpoints,
	Image_points & image_points,
	const translation_t & position,
	const rotation_t & rotation,
	/*bearingVectors_t & bearingVectors,*/
	size_t numberPoints,
	double depth,
	double noise);

void ShowSolutions(transformations_t & solutions);

//**************MY programe end*******************


}

#endif /* GV_EXPERIMENT_HELPERS_HPP_ */
