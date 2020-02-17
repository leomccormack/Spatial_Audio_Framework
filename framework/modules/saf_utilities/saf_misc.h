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

/*
 * Filename: saf_misc.h
 * --------------------
 * Miscellaneous functions
 *
 * Dependencies:
 *     none
 * Author, date created:
 *     Leo McCormack, 29.01.2020
 */

#ifndef SAF_MISC_H_INCLUDED
#define SAF_MISC_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Cross-platform sleep macro from (slightly modified):
 https://cboard.cprogramming.com/c-programming/170381-cross-platform-wait-sleep.html */
#ifdef _WIN32
    /* For Windows (32- and 64-bit) */
# include <windows.h>
# define SAF_SLEEP(msecs) Sleep(msecs)
#elif defined(__unix) || defined(__APPLE__)
    /* For linux, OSX, and other unixes */
# ifndef _POSIX_C_SOURCE
#  define _POSIX_C_SOURCE 199309L /* or greater */
# endif
# include <time.h>
# define SAF_SLEEP(msecs) do {            \
   struct timespec ts;             \
   ts.tv_sec = msecs/1000;         \
   ts.tv_nsec = msecs%1000*1000;   \
   nanosleep(&ts, NULL);           \
} while (0)
#else
# error "Unknown system"
#endif


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_MISC_H_INCLUDED */
