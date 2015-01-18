/* Mersenne Twister PRNG */
/* Copyright (c) 2015 Sugizaki Yukimasa; derived from libbrahe 1.3.2 */

/*
    Brahe is a heterogenous collection of mathematical tools,  written in Standard C.

    Copyright 2011 Scott Robert Ladd. All rights reserved.

    Brahe is user-supported open source software. Its continued development is dependent
    on financial support from the community. You can provide funding by visiting the Brahe
    website at:

        http://www.coyotegulch.com

    You may license Brahe in one of two fashions:

    1) Simplified BSD License (FreeBSD License)

    Redistribution and use in source and binary forms, with or without modification, are
    permitted provided that the following conditions are met:

    1.  Redistributions of source code must retain the above copyright notice, this list of
        conditions and the following disclaimer.

    2.  Redistributions in binary form must reproduce the above copyright notice, this list
        of conditions and the following disclaimer in the documentation and/or other materials
        provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY SCOTT ROBERT LADD ``AS IS'' AND ANY EXPRESS OR IMPLIED
    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
    FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SCOTT ROBERT LADD OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
    ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
    ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    The views and conclusions contained in the software and documentation are those of the
    authors and should not be interpreted as representing official policies, either expressed
    or implied, of Scott Robert Ladd.

    2) Closed-Source Proprietary License

    If your project is a closed-source or proprietary project, the Simplified BSD License may
    not be appropriate or desirable. In such cases, contact the Brahe copyright holder to
    arrange your purchase of an appropriate license.

    The author can be contacted at:

          scott.ladd@coyotegulch.com
          scott.ladd@gmail.com
          http:www.coyotegulch.com
*/

/*
    ORIGINAL ALGORITHM COPYRIGHT

    Copyright (C) 1997, 2002 Makoto Matsumoto and Takuji Nishimura.
    Any feedback is very welcome. For any question, comments, see
    http://www.math.keio.ac.jp/matumoto/emt.html or email
    matumoto@math.keio.ac.jp

    Much as this algorithm is popular, I find it slower than Marsaglia's
    multiply-with-carry generators (implemented below). Still, this is
    an excellent generator that provide a good alternative when testing
    with multiple PRNGs.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <inttypes.h>
#include <errno.h>

struct mtwister_state_t {
	uint32_t *m;
	uint32_t m_seed;
	int i;
};

void mtwister_init(struct mtwister_state_t *stp);
void mtwister_finalize(struct mtwister_state_t *stp);
uint32_t mtwister_next(struct mtwister_state_t *stp);

#define MTWISTER_N 624
#define MTWISTER_M 397
#define MTWISTER_MAG 0x9908b0dfUL
#define MTWISTER_UPPER_MASK 0x80000000UL
#define MTWISTER_LOWER_MASK 0x7fffffffUL

void mtwister_init(struct mtwister_state_t *stp)
{
	if (stp == NULL) {
		fprintf(stderr, "%s:%d: error: invalid input st\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	stp->m = malloc(MTWISTER_N * sizeof(uint32_t));
	if (stp->m == NULL) {
		fprintf(stderr, "%s:%d: MTWISTER_N=%d: error: failed to malloc stp->m\n", __FILE__, __LINE__, MTWISTER_N);
		exit(EXIT_FAILURE);
	}

	{
		int fd;
		unsigned int seed;
		ssize_t rc;

		if ((fd = open("/dev/urandom", O_RDONLY)) == -1) {
			fprintf(stderr, "%s:%d: error: open: /dev/urandom: %s\n", __FILE__, __LINE__, strerror(errno));
			exit(EXIT_FAILURE);
		}

		if ((rc = read(fd, &seed, sizeof(seed))) == -1) {
			fprintf(stderr, "%s:%d: error: read: %s\n", __FILE__, __LINE__, strerror(errno));
			exit(EXIT_FAILURE);
		}
		srandom(seed);

		if ((rc = close(fd)) == -1) {
			fprintf(stderr, "%s:%d: error: close: %s\n", __FILE__, __LINE__, strerror(errno));
			exit(EXIT_FAILURE);
		}
	}

	stp->m[0] = random();

	for (stp->i = 1; stp->i < MTWISTER_N; stp->i++)
		stp->m[stp->i] = 1812433253UL * (stp->m[stp->i - 1] ^ (stp->m[stp->i - 1] >> 30)) + stp->i;

	{
		int i;

		for (i = 0; i < 100; i++)
			(void) mtwister_next(stp);
	}
	
	return;
}

void mtwister_finalize(struct mtwister_state_t *stp)
{
	if (stp == NULL) {
		fprintf(stderr, "%s:%d: error: invalid input st\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	free(stp->m);

	return;
}

uint32_t mtwister_next(struct mtwister_state_t *stp)
{
	size_t kk;
	uint32_t v;

	if ((stp == NULL) || (stp->m == NULL)) {
		fprintf(stderr, "%s:%d: error: invalid input st\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	if (stp->i >= MTWISTER_N) {
		uint32_t tmp;

		for (kk = 0; kk < MTWISTER_N - MTWISTER_M; kk++) {
			tmp = (stp->m[kk] & MTWISTER_UPPER_MASK) | (stp->m[kk + 1] & MTWISTER_LOWER_MASK);
			stp->m[kk] = stp->m[kk + MTWISTER_M] ^ (tmp >> 1) ^ (MTWISTER_MAG * (tmp & 0x1));
		}
		for (; kk < MTWISTER_N - 1; kk++) {
			tmp = (stp->m[kk] & MTWISTER_UPPER_MASK) | (stp->m[kk + 1] & MTWISTER_LOWER_MASK);
			stp->m[kk] = stp->m[kk + (MTWISTER_M - MTWISTER_N)] ^ (tmp >> 1) ^ (MTWISTER_MAG * (tmp & 0x1));
		}

		tmp = (stp->m[MTWISTER_N - 1] & MTWISTER_UPPER_MASK) | (stp->m[0] & MTWISTER_LOWER_MASK);
		stp->m[MTWISTER_N - 1] = stp->m[MTWISTER_M - 1] ^ (tmp >> 1) ^ (MTWISTER_MAG * (tmp & 0x1));

		stp->i = 0;
	}

	v = stp->m[stp->i++];
	v ^= (v >> 11);
	v ^= (v << 7) & 0x9d2c5680UL;
	v ^= (v << 15) & 0xefc60000UL;
	v ^= (v >> 18);

	return v;
}

#ifdef TEST
int main()
{
	struct mtwister_state_t st;
	
	mtwister_init(&st);

	printf("0x%08"PRIx32"\n", mtwister_next(&st));

	mtwister_finalize(&st);

	return 0;
}
#endif /* TEST */
