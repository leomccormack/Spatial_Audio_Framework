/**
 * @file sldoa_database.h
 * @brief A spatially-localised active-intensity (SLAI) based direction-of-
 *        arrival estimator (SLDoA)
 * @author Leo McCormack
 * @date 12.02.2018
 * @license ISC
 */

#ifndef __SLDOA_DATABASE_INCLUDED__
#define __SLDOA_DATABASE_INCLUDED__

#define NUM_GRID_DIRS ( 2562 )
 
extern const double __grid_Y[64][2562];
extern const double __grid_dirs_deg[2562][2];

#endif /* __SLDOA_DATABASE_INCLUDED__ */
