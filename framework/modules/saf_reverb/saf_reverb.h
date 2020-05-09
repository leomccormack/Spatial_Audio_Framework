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
 * ...
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

#define IMS_MAX_NUM_SOURCES 100
#define IMS_MAX_NUM_RECEIVERS 100
typedef float* ims_rir;

/*
 * Here set up you scene parameters, source locations, microphone locations,
 * room boundaries and target reverberation times (or directly wall
 * absorption coefficients). If there is no room (anechoic rendering) the
 * global origin is arbitrary (e.g. can be at one of the microphones),
 * however if there is a room (reverberant rendering), all receiver and
 * source positions should be given with respect to the bottom left corner
 * of the room (top view), with positive x+ extending to the east, and
 * positive y+ extending to the north, while z+ is extending purpendicular
 * to them towards the viewer (right-hand rule)
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
 * Note that there is no checking for the source-microphone coordinates
 * falling inside the boundaries of the room.
 */

// create empty room. Hint: add a reciever and source
void ims_shoebox_create(void** phIms,
                        int length,
                        int width,
                        int height,
                        float* abs_wall,  /* Absorption coefficients for each octave band, and each wall; nOctBands x 6 */
                        float lowestOctaveBand,
                        int nOctBands,
                        float c_ms,
                        float fs);

void ims_shoebox_destroy(void** phIms);
 
void ims_shoebox_renderEchogramSH(void* hIms,
                                  float maxTime_ms,
                                  int sh_order);

void ims_shoebox_renderSHRIRs(void* hIms,
                              int fractionalDelaysFLAG);




/* add/remove/update functions: */

long ims_shoebox_addSource(void* hIms,
                           float position_xyz[3]);

long ims_shoebox_addReceiver(void* hIms,
                             float rec_xyz[3]);

void ims_shoebox_updateSource(void* hIms,
                              long sourceID,
                              float position_xyz[3]);

void ims_shoebox_updateReceiver(void* hIms,
                                long receiverID,
                                float position_xyz[3]);

void ims_shoebox_removeSource(void* hIms,
                              long sourceID);

void ims_shoebox_removeReceiver(void* hIms,
                                long receiverID);


#ifdef __cplusplus
} /* extern "C" */
#endif  /* __cplusplus */

#endif /* __SAF_REVERB_H_INCLUDED__ */
