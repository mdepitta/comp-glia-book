#ifndef PNTSTSH
#define PNTSTSH 1

#ifdef __cplusplus
	extern "C" {
#endif

#define MAXBLOCKS 10000

typedef point site;
typedef Coord* normalp;
point	site_blocks[MAXBLOCKS];
int	num_blocks;

#ifdef __cplusplus
}
#endif
