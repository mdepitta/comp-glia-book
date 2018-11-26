#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

/* from rand.c */

double double_rand(void);
void init_rand(long seed);


int main(int argc, char *argv[])

{	int nsites,d,i;
	if (!argv[1] || !argv[2]) exit(3);
	sscanf(argv[1], "%d", &nsites);
	sscanf(argv[2], "%d", &d);
	init_rand(0);
	while(nsites>0){
		for (i=0;i<d;i++) printf("%6.0f ", floor(1e5*double_rand()));
		printf("\n");
		nsites--;
	}
}
