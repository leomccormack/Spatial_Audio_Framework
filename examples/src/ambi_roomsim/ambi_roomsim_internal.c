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
 * @file ambi_roomsim_internal.c
 * @brief A simple shoebox room Ambisonic encoder.
 *
 * @author Leo McCormack
 * @date 10.08.2020
 * @license ISC
 */

#include "ambi_roomsim.h"
#include "ambi_roomsim_internal.h"
 
const float ambi_roomsim_default_abs_wall[6] = { 0.341055000f, 0.431295000f, 0.351295000f, 0.344335000f, 0.401775000f, 0.482095000f};

const float ambi_roomsim_default_room_dims[3] = { 9.1f, 8.0f, 3.0f };

const int ambi_roomsim_defaultNumSources = 1;

const float ambi_roomsim_defaultSourcePositions[ROOM_SIM_MAX_NUM_SOURCES][3] = {
    {5.2f, 1.5f, 1.4f}, {2.1f, 1.0f, 1.3f}, {3.1f, 5.0f, 2.3f}, {7.1f, 2.0f, 1.4f},
    {4.2f, 1.5f, 1.8f}, {3.1f, 1.0f, 1.2f}, {4.1f, 5.0f, 2.3f}, {6.1f, 2.0f, 1.4f},
    {3.2f, 1.5f, 2.2f}, {4.1f, 1.0f, 1.1f}, {5.1f, 5.0f, 2.3f}, {5.1f, 2.0f, 1.4f},
    {2.2f, 1.5f, 2.6f}, {5.1f, 1.0f, 1.0f}, {6.1f, 5.0f, 2.3f}, {4.1f, 2.0f, 1.4f}
};

const int ambi_roomsim_defaultNumReceivers = 1;

const float ambi_roomsim_defaultReceiverPositions[ROOM_SIM_MAX_NUM_RECEIVERS][3] = {
    {5.2f, 3.5f, 1.4f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f},
    {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f},
    {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f},
    {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}
};
