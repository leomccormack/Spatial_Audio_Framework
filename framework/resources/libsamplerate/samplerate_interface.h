
/* A custom minimal interface to libsamplerate, which condenses the library's
 * functionality to offer only what SAF requires.
 *
 * Internal changes include: raising maximum supported channels from 128 to 8192
 *
 * Leo McCormack, 2020
 */

#ifndef SAMPLERATE_INTERFACE_H
#define SAMPLERATE_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif	/* __cplusplus */

/** Mirror of the enum options sent internally to src_set_converter()
 * (DO NOT CHANGE TAG IDs!) */
typedef enum _RESAMPLE_QUALITY_OPTIONS
{
    RESAMPLE_BEST_QUALITY   = 0,
    RESAMPLE_MEDIUM_QUALITY = 1,
    RESAMPLE_FASTEST        = 2,
    RESAMPLE_ORDER_HOLD     = 3,
    RESAMPLE_LINEAR         = 4

}RESAMPLE_QUALITY_OPTIONS;

/**
 * Resamples
 *
 * @param[in]  insig  Input signal; FLAT: nChannels x length_insig
 * @param[in] <szfhdxgjhc
 */
void sampleratelib_resample(/* Input Arguments */
                            float* insig,
                            int length_insig,
                            int length_outsig,
                            int input_fs,
                            int output_fs,
                            int nChannels,
                            RESAMPLE_QUALITY_OPTIONS quality,
                            /* Output Arguments */
                            float* outsig);



#ifdef __cplusplus
}		/* extern "C" */
#endif	/* __cplusplus */

#endif	/* SAMPLERATE_INTERFACE_H */

