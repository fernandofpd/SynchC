// This is the ``Mersenne Twister'' random number generator MT19937, which
// generates pseudorandom integers uniformly distributed in 0..(2^32 - 1)
// starting from any odd seed in 0..(2^32 - 1).  This version is a recode
// by Shawn Cokus (Cokus@math.washington.edu) on March 8, 1998 of a version by
// Takuji Nishimura (who had suggestions from Topher Cooper and Marc Rieffel in
// July-August 1997).
// 

typedef unsigned long uint32;

#define NORM 4294967295
#define N (624)
#define M (397)
#define K (0x9908B0DFU)
#define hiBit(u) ((u) & 0x80000000U)
#define loBit(u) ((u) & 0x00000001U)
#define loBits(u) ((u) & 0x7FFFFFFFU)
#define mixBits(u,v) (hiBit(u)|loBits(v))

static uint32 state[N+1];
static uint32 *next;
static int left = -1;

void seedMT(uint32 seed);
uint32 reloadMT(void);
uint32 randomMT(void);
