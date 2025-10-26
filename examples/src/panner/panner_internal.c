/*
 * Copyright 2017-2018 Leo McCormack
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
 * @file panner_internal.c
 * @brief A frequency-dependent 3D panner based on the Vector-base Amplitude
 *        Panning (VBAP) method [1], with an optional spread control [2].
 *
 * Depending on the listening room, it may be beneficial to employ amplitude-
 * normalised gains for low frequencies, and energy-normalised gains for high
 * frequencies. Therefore, this VBAP implementation also uses the method
 * described in [3], to do just that.
 *
 * @see [1] Pulkki, V. (1997). Virtual sound source positioning using vector
 *          base amplitude panning. Journal of the audio engineering society,
 *          45(6), 456-466.
 * @see [2] Pulkki, V. (1999). Uniform spreading of amplitude panned virtual
 *          sources. In Proceedings of the 1999 IEEE Workshop on Applications of
 *          Signal Processing to Audio and Acoustics. WASPAA'99 (Cat. No.
 *          99TH8452) (pp. 187-190). IEEE.
 * @see [3] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V.
 *          (2014). Gain normalisation in amplitude panning as a function of
 *          frequency and room reverberance. 55th International Conference of
 *          the AES. Helsinki, Finland.
 *
 * @author Leo McCormack
 * @date 25.09.2017
 * @license ISC
 */

#include "panner.h"
#include "panner_internal.h"

const int panner_defaultNumSources = 1;

/* First 1 source(s) are for mono [0 0], the rest are for a 128 point sphere covering */
const float panner_defaultSourceDirections[MAX_NUM_INPUTS][2] =
{ { 0.0f,    0.0f},
  { -120.415474946769f,    69.7642225220975f},
  { 169.426697678052f,    -17.3466306179260f},
  { 107.629281568279f,    -10.0793856720614f},
  { -59.5674172539630f,    -69.7763953095265f},
  { -107.629493132161f,    10.0801907586439f},
  { 59.5651193014662f,    69.7750019669784f},
  { -72.3715136529922f,    -10.0909592824967f},
  { 120.418832045969f,    -69.7630732527928f},
  { 72.3705387467529f,    10.0908849630425f},
  { 10.5699823270029f,    17.3489915451677f},
  { -169.426851150349f,    17.3476277272366f},
  { -30.7728814706980f,    68.2546206344998f},
  { 101.529181734204f,    18.5653877408229f},
  { 149.250234431564f,    68.2619792676765f},
  { 18.9183978612362f,    -10.9197242075287f},
  { 78.4729799557400f,    -18.5558892226653f},
  { -18.9161109968480f,    10.9195275229078f},
  { -78.4723066816396f,    18.5565169424192f},
  { -161.085243560336f,    -10.9220975953767f},
  { -101.530063994217f,    -18.5648176173994f},
  { 161.085855454581f,    10.9237834974412f},
  { 30.7790837424989f,    -68.2562996173874f},
  { -149.255648697558f,    -68.2594892022576f},
  { 56.4638592172755f,    41.2645742763346f},
  { 46.4626904606488f,    24.5386896277818f},
  { -123.528839291756f,    41.2565889648213f},
  { 32.1932295461248f,    38.8027347956192f},
  { 133.534545076742f,    -24.5310453346374f},
  { -32.1935458377634f,    -38.8041545028931f},
  { -133.533477606995f,    24.5319420428588f},
  { -147.797881221742f,    38.7980990097138f},
  { -46.4626630742447f,    -24.5400749997299f},
  { 147.799544204027f,    -38.7965984357137f},
  { -56.4649576930965f,    -41.2658215291547f},
  { 123.530547454259f,    -41.2555316199835f},
  { 84.7385596473317f,    27.3132538556364f},
  { 27.4072860474693f,    4.67593317316088f},
  { -95.2612516264878f,    27.3036967398252f},
  { 10.0607591758466f,    62.2287426217158f},
  { 152.592710471291f,    -4.67007594116539f},
  { -10.0604693351594f,    -62.2317767049329f},
  { -152.592615409799f,    4.67188744845125f},
  { -169.918387889403f,    62.2272712110016f},
  { -27.4047652820400f,    -4.67684895484045f},
  { 169.918493477140f,    -62.2246416140826f},
  { -84.7390238025025f,    -27.3120311727342f},
  { 95.2622793653641f,    -27.3027426350066f},
  { 136.272846725200f,    -0.726139068263599f},
  { -1.05456977390436f,    -46.2698051615206f},
  { -43.7263256019408f,    -0.737666959763443f},
  { -91.0084089447853f,    43.7170295842697f},
  { -178.936314291048f,    46.2682278552731f},
  { 91.0100541861312f,    -43.7163112498585f},
  { 178.938490558618f,    -46.2666218304286f},
  { 88.9860902890987f,    43.7285680951606f},
  { 1.05095788619312f,    46.2682407577363f},
  { -88.9851906686018f,    -43.7278600822776f},
  { -136.273314433943f,    0.726484363994244f},
  { 43.7273191822901f,    0.733967099281741f},
  { 55.2321077221357f,    10.8187456245361f},
  { 13.0867814081219f,    34.0654299770448f},
  { -124.767898552651f,    10.8099467408727f},
  { 71.4751302534919f,    53.7984109774923f},
  { 166.905549455530f,    -34.0625862919868f},
  { -71.4765158530160f,    -53.8002805142714f},
  { -166.904587570338f,    34.0633858486074f},
  { -108.515564585463f,    53.7864186881779f},
  { -13.0869589355977f,    -34.0644403811701f},
  { 108.517274377512f,    -53.7853093207457f},
  { -55.2337689585088f,    -10.8200079671294f},
  { 124.767519417537f,    -10.8096092876561f},
  { -105.493359369733f,    -68.1345228940484f},
  { -111.151973933085f,    -5.71277564515985f},
  { 74.5244447329726f,    -68.1205741555505f},
  { -173.885132403874f,    -21.0401043953034f},
  { -68.8521018079939f,    5.70313885775370f},
  { 173.885050918764f,    21.0419786777516f},
  { 68.8513569996117f,    -5.70297394874866f},
  { 6.12041448982713f,    -21.0387333986971f},
  { 111.151515131908f,    5.71357243375474f},
  { -6.11913621739966f,    21.0397756114589f},
  { 105.487055344890f,    68.1331981897040f},
  { -74.5188619841755f,    68.1214716326340f},
  { 35.2822222163987f,    -15.1801841162885f},
  { -25.1732516880600f,    51.9826514887145f},
  { -144.722913967824f,    -15.1859965048827f},
  { 108.386507250934f,    33.8824333892584f},
  { -154.839566662519f,    -51.9856067440744f},
  { -108.385758587361f,    -33.8826426573830f},
  { 154.840828582605f,    51.9884126722963f},
  { -71.6110721306039f,    33.8714084774643f},
  { 25.1699762263473f,    -51.9843011897216f},
  { 71.6126173250777f,    -33.8705490198003f},
  { -35.2810594358470f,    15.1763090469469f},
  { 144.722474074749f,    15.1863112389402f},
  { -125.277803641269f,    -28.5561977838671f},
  { -146.317489805034f,    -30.4888082061806f},
  { 54.7220842889970f,    -28.5463137051834f},
  { -133.292600372594f,    -45.8199818171681f},
  { -33.6854220654234f,    30.4794807422719f},
  { 133.295169860873f,    45.8200619710773f},
  { 33.6830420572493f,    -30.4828689919959f},
  { 46.7098267598019f,    -45.8106991831118f},
  { 146.318736830907f,    30.4889588527258f},
  { -46.7107531920589f,    45.8077510168776f},
  { 125.278543193551f,    28.5567151663926f},
  { -54.7226679373581f,    28.5424675755949f},
  { -144.402199049652f,    54.7141577397485f},
  { 112.382216591782f,    -28.0100385212576f},
  { 35.5834085793474f,    54.7193853715013f},
  { -29.9155874729000f,    -19.6471282266303f},
  { 67.6163420587218f,    28.0211878299760f},
  { 29.9151078604395f,    19.6461577548916f},
  { -67.6164804090156f,    -28.0221553885221f},
  { 150.081679085309f,    -19.6396879622134f},
  { -112.381365317985f,    28.0105570315277f},
  { -150.081127534239f,    19.6412999777121f},
  { 144.402936805264f,    -54.7127597303097f},
  { -35.5837994140428f,    -54.7203576209200f},
  { 68.5348504928164f,    -52.8527273534571f},
  { -54.8218918027732f,    12.7570858971645f},
  { -111.464947340552f,    -52.8652964226689f},
  { 164.510735700738f,    34.1892226850382f},
  { -125.181243085238f,    -12.7710900404128f},
  { -164.510172672569f,    -34.1865270720665f},
  { 125.180742749445f,    12.7714932994183f},
  { -15.4967238792034f,    34.1855682977571f}};

const int panner_defaultNumLoudspeakers = 2;

/* First 2 loudspeakers are for stereo, the rest are for a 128 point sphere covering */
const float panner_defaultLoudspeakerDirections[MAX_NUM_OUTPUTS][2] =
{ { 30.0f,    0.0f},
  { -30.0f,    0.0f},
  { 169.426697678052f,    -17.3466306179260f},
  { 107.629281568279f,    -10.0793856720614f},
  { -59.5674172539630f,    -69.7763953095265f},
  { -107.629493132161f,    10.0801907586439f},
  { 59.5651193014662f,    69.7750019669784f},
  { -72.3715136529922f,    -10.0909592824967f},
  { 120.418832045969f,    -69.7630732527928f},
  { 72.3705387467529f,    10.0908849630425f},
  { 10.5699823270029f,    17.3489915451677f},
  { -169.426851150349f,    17.3476277272366f},
  { -30.7728814706980f,    68.2546206344998f},
  { 101.529181734204f,    18.5653877408229f},
  { 149.250234431564f,    68.2619792676765f},
  { 18.9183978612362f,    -10.9197242075287f},
  { 78.4729799557400f,    -18.5558892226653f},
  { -18.9161109968480f,    10.9195275229078f},
  { -78.4723066816396f,    18.5565169424192f},
  { -161.085243560336f,    -10.9220975953767f},
  { -101.530063994217f,    -18.5648176173994f},
  { 161.085855454581f,    10.9237834974412f},
  { 30.7790837424989f,    -68.2562996173874f},
  { -149.255648697558f,    -68.2594892022576f},
  { 56.4638592172755f,    41.2645742763346f},
  { 46.4626904606488f,    24.5386896277818f},
  { -123.528839291756f,    41.2565889648213f},
  { 32.1932295461248f,    38.8027347956192f},
  { 133.534545076742f,    -24.5310453346374f},
  { -32.1935458377634f,    -38.8041545028931f},
  { -133.533477606995f,    24.5319420428588f},
  { -147.797881221742f,    38.7980990097138f},
  { -46.4626630742447f,    -24.5400749997299f},
  { 147.799544204027f,    -38.7965984357137f},
  { -56.4649576930965f,    -41.2658215291547f},
  { 123.530547454259f,    -41.2555316199835f},
  { 84.7385596473317f,    27.3132538556364f},
  { 27.4072860474693f,    4.67593317316088f},
  { -95.2612516264878f,    27.3036967398252f},
  { 10.0607591758466f,    62.2287426217158f},
  { 152.592710471291f,    -4.67007594116539f},
  { -10.0604693351594f,    -62.2317767049329f},
  { -152.592615409799f,    4.67188744845125f},
  { -169.918387889403f,    62.2272712110016f},
  { -27.4047652820400f,    -4.67684895484045f},
  { 169.918493477140f,    -62.2246416140826f},
  { -84.7390238025025f,    -27.3120311727342f},
  { 95.2622793653641f,    -27.3027426350066f},
  { 136.272846725200f,    -0.726139068263599f},
  { -1.05456977390436f,    -46.2698051615206f},
  { -43.7263256019408f,    -0.737666959763443f},
  { -91.0084089447853f,    43.7170295842697f},
  { -178.936314291048f,    46.2682278552731f},
  { 91.0100541861312f,    -43.7163112498585f},
  { 178.938490558618f,    -46.2666218304286f},
  { 88.9860902890987f,    43.7285680951606f},
  { 1.05095788619312f,    46.2682407577363f},
  { -88.9851906686018f,    -43.7278600822776f},
  { -136.273314433943f,    0.726484363994244f},
  { 43.7273191822901f,    0.733967099281741f},
  { 55.2321077221357f,    10.8187456245361f},
  { 13.0867814081219f,    34.0654299770448f},
  { -124.767898552651f,    10.8099467408727f},
  { 71.4751302534919f,    53.7984109774923f},
  { 166.905549455530f,    -34.0625862919868f},
  { -71.4765158530160f,    -53.8002805142714f},
  { -166.904587570338f,    34.0633858486074f},
  { -108.515564585463f,    53.7864186881779f},
  { -13.0869589355977f,    -34.0644403811701f},
  { 108.517274377512f,    -53.7853093207457f},
  { -55.2337689585088f,    -10.8200079671294f},
  { 124.767519417537f,    -10.8096092876561f},
  { -105.493359369733f,    -68.1345228940484f},
  { -111.151973933085f,    -5.71277564515985f},
  { 74.5244447329726f,    -68.1205741555505f},
  { -173.885132403874f,    -21.0401043953034f},
  { -68.8521018079939f,    5.70313885775370f},
  { 173.885050918764f,    21.0419786777516f},
  { 68.8513569996117f,    -5.70297394874866f},
  { 6.12041448982713f,    -21.0387333986971f},
  { 111.151515131908f,    5.71357243375474f},
  { -6.11913621739966f,    21.0397756114589f},
  { 105.487055344890f,    68.1331981897040f},
  { -74.5188619841755f,    68.1214716326340f},
  { 35.2822222163987f,    -15.1801841162885f},
  { -25.1732516880600f,    51.9826514887145f},
  { -144.722913967824f,    -15.1859965048827f},
  { 108.386507250934f,    33.8824333892584f},
  { -154.839566662519f,    -51.9856067440744f},
  { -108.385758587361f,    -33.8826426573830f},
  { 154.840828582605f,    51.9884126722963f},
  { -71.6110721306039f,    33.8714084774643f},
  { 25.1699762263473f,    -51.9843011897216f},
  { 71.6126173250777f,    -33.8705490198003f},
  { -35.2810594358470f,    15.1763090469469f},
  { 144.722474074749f,    15.1863112389402f},
  { -125.277803641269f,    -28.5561977838671f},
  { -146.317489805034f,    -30.4888082061806f},
  { 54.7220842889970f,    -28.5463137051834f},
  { -133.292600372594f,    -45.8199818171681f},
  { -33.6854220654234f,    30.4794807422719f},
  { 133.295169860873f,    45.8200619710773f},
  { 33.6830420572493f,    -30.4828689919959f},
  { 46.7098267598019f,    -45.8106991831118f},
  { 146.318736830907f,    30.4889588527258f},
  { -46.7107531920589f,    45.8077510168776f},
  { 125.278543193551f,    28.5567151663926f},
  { -54.7226679373581f,    28.5424675755949f},
  { -144.402199049652f,    54.7141577397485f},
  { 112.382216591782f,    -28.0100385212576f},
  { 35.5834085793474f,    54.7193853715013f},
  { -29.9155874729000f,    -19.6471282266303f},
  { 67.6163420587218f,    28.0211878299760f},
  { 29.9151078604395f,    19.6461577548916f},
  { -67.6164804090156f,    -28.0221553885221f},
  { 150.081679085309f,    -19.6396879622134f},
  { -112.381365317985f,    28.0105570315277f},
  { -150.081127534239f,    19.6412999777121f},
  { 144.402936805264f,    -54.7127597303097f},
  { -35.5837994140428f,    -54.7203576209200f},
  { 68.5348504928164f,    -52.8527273534571f},
  { -54.8218918027732f,    12.7570858971645f},
  { -111.464947340552f,    -52.8652964226689f},
  { 164.510735700738f,    34.1892226850382f},
  { -125.181243085238f,    -12.7710900404128f},
  { -164.510172672569f,    -34.1865270720665f},
  { 125.180742749445f,    12.7714932994183f},
  { -15.4967238792034f,    34.1855682977571f}};

void panner_setCodecStatus(void* const hPan, CODEC_STATUS newStatus)
{
    panner_data *pData = (panner_data*)(hPan);
    if(newStatus==CODEC_STATUS_NOT_INITIALISED){
        /* Pause until current initialisation is complete */
        while(pData->codecStatus == CODEC_STATUS_INITIALISING)
            SAF_SLEEP(10);
    }
    pData->codecStatus = newStatus;
}

void panner_initGainTables(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
#ifndef FORCE_3D_LAYOUT
    int i;
    float sum_elev;
    
    /* determine dimensionality */
    sum_elev = 0.0f;
    for(i=0; i<pData->nLoudpkrs; i++)
        sum_elev += fabsf(pData->loudpkrs_dirs_deg[i][1]); 
    if(sum_elev < 0.01f)
        pData->output_nDims = 2;
    else
        pData->output_nDims = 3;
#endif
    
    strcpy(pData->progressBarText,"Generating VBAP Gain Table");
    pData->progressBar0_1 = 0.7f;
    
    /* generate VBAP gain table */
    free(pData->vbap_gtable);
    pData->vbapTableRes[0] = 1;
    pData->vbapTableRes[1] = 1;
#ifdef FORCE_3D_LAYOUT
    pData->output_nDims = 3;
    generateVBAPgainTable3D((float*)pData->loudpkrs_dirs_deg, pData->nLoudpkrs, pData->vbapTableRes[0], pData->vbapTableRes[1], 1, 1, pData->spread_deg,
                            &(pData->vbap_gtable), &(pData->N_vbap_gtable), &(pData->nTriangles));
#else
    if(pData->output_nDims==2)
        generateVBAPgainTable2D((float*)pData->loudpkrs_dirs_deg, pData->nLoudpkrs, pData->vbapTableRes[0],
                                &(pData->vbap_gtable), &(pData->N_vbap_gtable), &(pData->nTriangles));
    else{
        generateVBAPgainTable3D((float*)pData->loudpkrs_dirs_deg, pData->nLoudpkrs, pData->vbapTableRes[0], pData->vbapTableRes[1], 1, 1, pData->spread_deg,
                                &(pData->vbap_gtable), &(pData->N_vbap_gtable), &(pData->nTriangles));
        if(pData->vbap_gtable==NULL){
            /* if generating vbap gain tabled failed, re-calculate with 2D VBAP */
            pData->output_nDims = 2;
            panner_initGainTables(hPan);
        }
    }
#endif
}

void panner_initTFT
(
    void* const hPan
)
{
    panner_data *pData = (panner_data*)(hPan);
    
    if(pData->hSTFT==NULL)
        afSTFT_create(&(pData->hSTFT), pData->new_nSources, pData->new_nLoudpkrs, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    else if (pData->new_nSources!=pData->nSources || pData->new_nLoudpkrs!=pData->nLoudpkrs){
        afSTFT_channelChange(pData->hSTFT, pData->new_nSources, pData->new_nLoudpkrs);
        afSTFT_clearBuffers(pData->hSTFT); 
    }
    pData->nSources = pData->new_nSources;
    pData->nLoudpkrs = pData->new_nLoudpkrs;
}

void panner_loadSourcePreset
(
    SOURCE_CONFIG_PRESETS preset,
    _Atomic_FLOAT32 dirs_deg[MAX_NUM_INPUTS][2],
    _Atomic_INT32* newNCH,
    int* nDims
)
{
    float sum_elev;
    int ch, i, nCH;
    
    switch(preset){
        default:
        case SOURCE_CONFIG_PRESET_DEFAULT:
            nCH = 1;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = 0.0f;
            break;
        case SOURCE_CONFIG_PRESET_MONO:
            nCH = 1;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __mono_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_STEREO:
            nCH = 2;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __stereo_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_5PX:
            nCH = 5;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __5pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_7PX:
            nCH = 7;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __7pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_8PX:
            nCH = 8;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __8pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_9PX:
            nCH = 9;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __9pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_10PX:
            nCH = 10;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __10pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_11PX:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_11PX_7_4:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_7_4_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_13PX:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __13pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_22PX:
            nCH = 22;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __22pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_22P2_9_10_3:
            nCH = 24;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __9_10_3p2_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_AALTO_MCC:
            nCH = 45;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_MCC_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_AALTO_MCC_SUBSET:
            nCH = 37;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_MCCsubset_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_AALTO_APAJA:
            nCH = 29;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_Apaja_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_AALTO_LR:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_LR_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_DTU_AVIL:
            nCH = 64;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __DTU_AVIL_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_4:
            nCH = 4;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_2_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_12:
            nCH = 12;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_4_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_24:
            nCH = 24;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_6_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_36:
            nCH = 36;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_8_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_48:
            nCH = 48;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_9_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_60:
            nCH = 60;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_10_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_9:
            nCH = 9;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_9_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_16:
            nCH = 16;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_16_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_25:
            nCH = 25;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_25_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_49:
            nCH = 49;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_49_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_64:
            nCH = 64;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_64_dirs_deg[ch][i];
            break;
    }
    
    /* Fill remaining slots with default coords */
    for(; ch<MAX_NUM_INPUTS; ch++){
        for(i=0; i<2; i++){
            dirs_deg[ch][i] = __default_LScoords128_deg[ch][i];
        }
    }
    
    /* For dynamically changing the number of TFT channels */
    (*newNCH) = nCH;
    
    /* estimate number of dimensions. (Obviously fails if using 2D setups thare are on an angle.
       However, in these cases, triangulation should fail and revert to 2D anyway) */
    sum_elev = 0.0f;
    for(i=0; i<nCH; i++)
        sum_elev += fabsf((float)dirs_deg[i][1]);
    if(sum_elev < 0.01f)
        (*nDims) = 2;
    else
        (*nDims) = 3;
}

void panner_loadLoudspeakerPreset
(
    LOUDSPEAKER_ARRAY_PRESETS preset,
    _Atomic_FLOAT32 dirs_deg[MAX_NUM_INPUTS][2],
    _Atomic_INT32* newNCH,
    int* nDims
)
{
    float sum_elev;
    int ch, i, nCH;

    switch(preset){
        default:
        case LOUDSPEAKER_ARRAY_PRESET_DEFAULT:
        case LOUDSPEAKER_ARRAY_PRESET_STEREO:
            nCH = 2;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __stereo_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_5PX:
            nCH = 5;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __5pX_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_7PX:
            nCH = 7;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __7pX_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_8PX:
            nCH = 8;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __8pX_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_9PX:
            nCH = 9;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __9pX_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_10PX:
            nCH = 10;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __10pX_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_11PX:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_11PX_7_4:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_7_4_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_13PX:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __13pX_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_22PX:
            nCH = 22;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __22pX_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_22P2_9_10_3:
            nCH = 24;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __9_10_3p2_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_AALTO_MCC:
            nCH = 45;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_MCC_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_AALTO_MCC_SUBSET:
            nCH = 37;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_MCCsubset_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_AALTO_APAJA:
            nCH = 29;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_Apaja_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_AALTO_LR:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_LR_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_DTU_AVIL:
            nCH = 64;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __DTU_AVIL_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_4:
            nCH = 4;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_2_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_12:
            nCH = 12;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_4_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_24:
            nCH = 24;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_6_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_36:
            nCH = 36;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_8_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_48:
            nCH = 48;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_9_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_60:
            nCH = 60;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_10_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_SPH_COV_9:
            nCH = 9;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_9_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_SPH_COV_16:
            nCH = 16;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_16_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_SPH_COV_25:
            nCH = 25;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_25_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_SPH_COV_49:
            nCH = 49;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_49_dirs_deg[ch][i];
            break;
        case LOUDSPEAKER_ARRAY_PRESET_SPH_COV_64:
            nCH = 64;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_64_dirs_deg[ch][i];
            break;
    }

    /* Fill remaining slots with default coords */
    for(; ch<MAX_NUM_INPUTS; ch++){
        for(i=0; i<2; i++){
            dirs_deg[ch][i] = __default_LScoords128_deg[ch][i];
        }
    }

    /* For dynamically changing the number of TFT channels */
    (*newNCH) = nCH;

    /* estimate number of dimensions. (Obviously fails if using 2D setups thare are on an angle.
       However, in these cases, triangulation should fail and revert to 2D anyway) */
    sum_elev = 0.0f;
    for(i=0; i<nCH; i++)
        sum_elev += fabsf((float)dirs_deg[i][1]);
    if(sum_elev < 0.01f)
        (*nDims) = 2;
    else
        (*nDims) = 3;
}
 
