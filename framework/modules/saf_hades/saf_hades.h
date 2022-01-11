/*
 * This file is part of the saf_hades module.
 * Copyright (c) 2021 - Leo McCormack & Janani Fernandez
 *
 * The saf_hades module is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * The saf_hades module is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * See <http://www.gnu.org/licenses/> for a copy of the GNU General Public
 * License.
 */

/**
 *@addtogroup HADES
 *@{
 * @file saf_hades.h
 * @brief Main header for the HADES module (#SAF_HADES_MODULE)
 *
 * The framework for binaural rendering of Hearing-Assistive/Augmented-reality
 * Devices (HADES) is described further in [1]
 *
 * @see [1] paper submitted for review.
 *
 * @author Leo McCormack and Janani Fernandez
 * @date 01.02.2021
 * @license GNU GPLv2
 */

#ifndef __SAF_HADES_H_INCLUDED__
#define __SAF_HADES_H_INCLUDED__

/*
 * This framework for binaural rendering of Hearing-Assistive/Augmented-reality
 * Devices (HADES) is divided into two main stages: analysis and synthesis.
 *
 * In the analysis, spatial parameters are estimated across time and frequency
 * based on the input microphone array signals.
 */
# include "saf_hades_analysis.h"

/*
 * In the synthesis, the output binaural signals are synthesised based on the
 * parameter and signal containers obtained from the analysis.
 */
# include "saf_hades_synthesis.h"


#endif /* __SAF_HADES_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup HADES */
