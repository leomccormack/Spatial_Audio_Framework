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
 * @file saf_reverb.h
 * @brief Public part of the reverb processing module (saf_reverb)
 *
 * A collection of reverb and room simulation algorithms.
 *
 * @author Leo McCormack
 * @date 06.05.2020
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
 * Note that this simulator is based on the shoebox-roomsim Matlab code by
 * Archontis Politis: https://github.com/polarch/shoebox-roomsim
 */

#define IMS_MAX_NUM_SOURCES 100
#define IMS_MAX_NUM_RECEIVERS 10

/**
 * Output format of the rendered room impulse responses (RIR)
 */
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
 * @warning There is currently no checking whether the source/receiver
 *          coordinates fall inside the boundaries of the room!
 *
 * @param[in] phIms            (&) address of the ims_shoebox handle
 * @param[in] length           Length of the room in meters
 * @param[in] width            Width of the room in meters
 * @param[in] height           Height of the room in meters
 * @param[in] abs_wall         Absorption coefficents per octave band and wall;
 *                             FLAT: nOctBands x 6
 * @param[in] lowestOctaveBand lowest octave band centre freq, in Hz (e.g. 125)
 * @param[in] nOctBands        Number of octave bands (i.e. doublings of
 *                             "lowestOctaveBand")
 * @param[in] c_ms             Speed of sound, meters per second
 * @param[in] fs               SampleRate, Hz
 */
void ims_shoebox_create(void** phIms,
                        int length,
                        int width,
                        int height,
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
 * @note The echograms are only updated if needed, so it is "OK" to call this
 *       function as many times as you wish, since there will be virtually no
 *       CPU overhead incurred if no update is required.
 *
 * @param[in] hIms      ims_shoebox handle
 * @param[in] maxTime_s Maximum length of time to compute the echograms, seconds
 */
void ims_shoebox_computeEchograms(void* hIms,
                                  float maxTime_s);

/**
 * Renders spherical harmonic room impulse responses for all active source/
 * receiver combinations
 *
 * @note This function is not intended to be used for real-time dynamic scenes,
 *       rather, it is more suited for high-quality static scene use cases. Use
 *       ims_shoebox_applyEchogramTD for real-time dynamic use cases.
 *
 * @param[in] hIms                 ims_shoebox handle
 * @param[in] fractionalDelaysFLAG 0: disabled, 1: use Lagrange interpolation
 */
void ims_shoebox_renderSHRIRs(void* hIms,
                              int fractionalDelaysFLAG);

/**
 * Applies the currently computed echograms in the time-domain, for all
 * sources and one specified receiver
 *
 * Note the following:
 *  - The signal pointers must be valid, and have allocated enough memory for
 *    the number of channels and be of (at least) nSamples in length.
 *  - The given receiverID must exist in the simulation. If it does not, then an
 *    assertion error is triggered.
 *
 * @note This function is intended for real-time dynamic scene use cases.
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


/* ====================== Add/Remove/Update functions ======================= */

/**
 * Adds a source object to the simulator, and returns a unique ID corresponding
 * to it
 *
 * @warning There is currently no checking whether the source/receiver
 *          coordinates fall inside the boundaries of the room!
 * @warning You are responsible for making sure that the signal pointers
 *          actually point to allocated memory of sufficient size, before
 *          calling ims_shoebox_applyEchogramTD!
 *
 * @param[in] hIms         ims_shoebox handle
 * @param[in] position_xyz Starting source position, in metres, x,y,z
 * @param[in] pSrc_sig     (&) address of the pointer to the 1-D input buffer
 *                         for this source: nSamples x 1
 *
 * @returns A unique ID corresponding to this source object
 */
long ims_shoebox_addSource(void* hIms,
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
 *          calling ims_shoebox_applyEchogramTD!
 *
 * @param[in] hIms         ims_shoebox handle
 * @param[in] sh_order     Spherical harmonic order of the receiver
 * @param[in] position_xyz Starting receiver position, in metres, x,y,z
 * @param[in] pSH_sigs     (&) address of the pointer to the 2-D output buffer
 *                         for this receiver: (sh_order+1)^2 x nSamples
 *
 * @returns A unique ID corresponding to this receiver object
 */
long ims_shoebox_addReceiverSH(void* hIms,
                               int sh_order,
                               float position_xyz[3],
                               float*** pSH_sigs);

/**
 * Updates the position of a specific source in the simulation
 */
void ims_shoebox_updateSource(void* hIms,
                              long sourceID,
                              float position_xyz[3]);

/**
 * Updates the position of a specific receiver in the simulation
 */
void ims_shoebox_updateReceiver(void* hIms,
                                long receiverID,
                                float position_xyz[3]);

/**
 * Removes a specific source from the simulation
 *
 * @note This does NOT free the source signal 'pSrc_sig'.
 */
void ims_shoebox_removeSource(void* hIms,
                              long sourceID);

/**
 * Removes a specific receiver from the simulation
 *
 * @note This does NOT free the receiver signals 'pSH_sigs'.
 */
void ims_shoebox_removeReceiver(void* hIms,
                                long receiverID);


#ifdef __cplusplus
} /* extern "C" */
#endif  /* __cplusplus */

#endif /* __SAF_REVERB_H_INCLUDED__ */
