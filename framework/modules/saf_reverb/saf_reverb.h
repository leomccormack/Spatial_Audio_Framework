/*
 * Copyright 2020 Leo McCormack
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

/**
 *@addtogroup Reverb
 *@{
 * @file saf_reverb.h
 * @brief Main header for the reverb processing module (#SAF_REVERB_MODULE)
 *
 * A collection of reverb and room simulation algorithms.
 *
 * @author Leo McCormack
 * @date 06.05.2020
 * @license ISC
 */

#ifndef __SAF_REVERB_H_INCLUDED__
#define __SAF_REVERB_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                         IMS Shoebox Room Simulator                         */
/* ========================================================================== */
/*
 * Note that this simulator is based on the shoebox-roomsim MATLAB library
 * found here: https://github.com/polarch/shoebox-roomsim
 * (Copyright (c) 2017, Archontis Politis, BSD-3-Clause License)
 */

/**
 * TODO: Arbitrary array receiver option; Directional source option; finish
 * ims_shoebox_renderRIRs() function;
 */

/** Maximum number of sources supported by an instance of the IMS simulator */
#define IMS_MAX_NUM_SOURCES 128

/** Maximum number of receivers supported by an instance of the IMS simulator */
#define IMS_MAX_NUM_RECEIVERS 16

/** Output format of the rendered room impulse responses (RIR) */
typedef struct _ims_rir{
    float* data;        
    int length, nChannels;
} ims_rir;

/**
 * Creates an instance of ims_shoebox room simulator
 *
 * Here you first set up the scene parameters, room boundaries, and the wall
 * absorption coefficients per octave band.
 *
 * Note that the room is initialised to be empty. Therefore, use the
 * ims_shoebox_addSource and ims_shoebox_addReceiverX functions to add sources
 * and recievers to the simulator. The source/receiver positions should be given
 * with respect to the bottom left corner of the room (top view), with positive
 * x+ extending to the east, and positive y+ extending to the north, while z+ is
 * extending purpendicular to them towards the viewer (right-hand rule).
 *
 * \verbatim
 *
 *   length/width
 *   |----------|
 *   ^ y           .
 *   |    ^z      /height
 *   |   /       /
 *   .__/_______.           _
 *   | /        |           |
 *   |/         |           | width/length
 *   o__________.------> x  _
 *
 *\endverbatim
 *
 * @warning There is currently no checking whether the source/receiver
 *          coordinates fall inside the boundaries of the room!
 *
 * @test test__ims_shoebox_RIR(), test__ims_shoebox_TD()
 *
 * @param[in] phIms            (&) address of the ims_shoebox handle
 * @param[in] roomDimensions   Room Length x Width x Height, in meters; 3 x 1
 * @param[in] abs_wall         Absorption coefficents per octave band and wall;
 *                             FLAT: nOctBands x 6
 * @param[in] lowestOctaveBand lowest octave band centre freq, in Hz (e.g. 125)
 * @param[in] nOctBands        Number of octave bands (i.e. doublings of
 *                             "lowestOctaveBand")
 * @param[in] c_ms             Speed of sound, meters per second
 * @param[in] fs               SampleRate, Hz
 */
void ims_shoebox_create(void** phIms,
                        float roomDimensions[3],
                        float* abs_wall,
                        float lowestOctaveBand,
                        int nOctBands,
                        float c_ms,
                        float fs);

/**
 * Destroys an instance of ims_shoebox room simulator
 *
 * @param[in] phIms (&) address of the ims_shoebox handle
 */
void ims_shoebox_destroy(void** phIms);

/**
 * Computes echograms for all active source/receiver combinations
 *
 * The sources are omnidirectional point sources, whereas the receiver will
 * have the directivity of whatever they are configured to have
 *
 * @note Set either the maximum reflection order (maxN) or the maximum IR length
 *       in seconds (maxTime_s). The option you don't want to use: set to <0.
 * @note The echograms are only updated if needed, so it is OK to call this
 *       function as many times as you wish, since there will be virtually no
 *       CPU overhead incurred if no update is required.
 *
 * @param[in] hIms      ims_shoebox handle
 * @param[in] maxN      Maximum reflection order
 * @param[in] maxTime_s Maximum length of time to compute the echograms, seconds
 */
void ims_shoebox_computeEchograms(void* hIms,
                                  int maxN,
                                  float maxTime_s);

/**
 * Renders room impulse responses for all active source/receiver combinations
 *
 * @param[in] hIms                 ims_shoebox handle
 * @param[in] fractionalDelaysFLAG 0: disabled, 1: use Lagrange interpolation
 */
void ims_shoebox_renderRIRs(void* hIms,
                            int fractionalDelaysFLAG);

/**
 * Applies the currently computed echograms in the time-domain, for all
 * sources, for one specified receiver
 *
 * Note the following:
 *  - The signal pointers for all the sources and the specified receiver must be
 *    valid, and have allocated enough memory for the number of channels and be
 *    of (at least) nSamples in length.
 *  - The given receiverID must exist in the simulation. If it does not, then an
 *    assertion error is triggered.
 *
 * @param[in] hIms                 ims_shoebox handle
 * @param[in] receiverID           ID of the receiver you wish to render
 * @param[in] nSamples             Number of samples to process
 * @param[in] fractionalDelaysFLAG 0: disabled, 1: use Lagrange interpolation
 */
void ims_shoebox_applyEchogramTD(/* Input Arguments */
                                 void* hIms,
                                 long receiverID,
                                 int nSamples,
                                 int fractionalDelaysFLAG);


/* =========================== Set/Get functions ============================ */

/** Sets new room dimensions */
void ims_shoebox_setRoomDimensions(void* hIms,
                                   float new_roomDimensions[3]);

/** Sets new wall absorption coefficients per wall and per band */
void ims_shoebox_setWallAbsCoeffs(void* hIms,
                                  float* abs_wall);


/* ================== Add/Remove/Update Objects functions ==================== */

/**
 * Adds a source object to the simulator, and returns a unique ID corresponding
 * to it
 *
 * @warning There is currently no checking whether the source/receiver
 *          coordinates fall inside the boundaries of the room!
 * @warning You are responsible for making sure that the signal pointers
 *          actually point to allocated memory of sufficient size, before
 *          calling ims_shoebox_applyEchogramTD()!
 *
 * @param[in] hIms         ims_shoebox handle
 * @param[in] position_xyz Starting source position, in metres, x,y,z
 * @param[in] pSrc_sig     (&) address of the pointer to the 1-D input buffer
 *                         for this source: nSamples x 1
 *
 * @returns A unique ID corresponding to this source object
 */
int ims_shoebox_addSource(void* hIms,
                          float position_xyz[3],
                          float** pSrc_sig);

/**
 * Adds a spherical harmonic (SH) receiver object to the simulator of a given
 * order, and returns a unique ID corresponding to it
 *
 * @warning There is currently no checking whether the source/receiver
 *          coordinates fall inside the boundaries of the room!
 * @warning You are responsible for making sure that the signal pointers
 *          actually point to allocated memory of sufficient size, before
 *          calling ims_shoebox_applyEchogramTD()!
 *
 * @param[in] hIms         ims_shoebox handle
 * @param[in] sh_order     Spherical harmonic order of the receiver
 * @param[in] position_xyz Starting receiver position, in metres, x,y,z
 * @param[in] pSH_sigs     (&) address of the pointer to the 2-D output buffer
 *                         for this receiver: (sh_order+1)^2 x nSamples
 *
 * @returns A unique ID corresponding to this receiver object
 */
int ims_shoebox_addReceiverSH(void* hIms,
                              int sh_order,
                              float position_xyz[3],
                              float*** pSH_sigs);

/** Updates the position of a specific source in the simulation */
void ims_shoebox_updateSource(void* hIms,
                              int sourceID,
                              float position_xyz[3]);

/** Updates the position of a specific receiver in the simulation */
void ims_shoebox_updateReceiver(void* hIms,
                                int receiverID,
                                float position_xyz[3]);

/**
 * Removes a specific source from the simulation
 *
 * @note This does NOT free the source signal pointer 'pSrc_sig'.
 */
void ims_shoebox_removeSource(void* hIms,
                              int sourceID);

/**
 * Removes a specific receiver from the simulation
 *
 * @note This does NOT free the receiver signals pointer 'pSH_sigs'.
 */
void ims_shoebox_removeReceiver(void* hIms,
                                int receiverID);


#ifdef __cplusplus
} /* extern "C" */
#endif  /* __cplusplus */

#endif /* __SAF_REVERB_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup Reverb */
