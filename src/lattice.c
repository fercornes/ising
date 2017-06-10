#include "lattice.h"
#include <stdio.h>
#include <stdlib.h>

int fill_lattice(int *lattice, int n, float p) {
	int il;
	for (il = 0; il < n*n; il++){
		double e;
		e = (double)rand()/(double)RAND_MAX;
		if (e <= p){
			lattice[il]=1;
		}
		else {
			lattice[il]=-1;
		}
	}
 return 0;
}

int print_lattice(int *lattice, int n) {
	int ic,jc;
	printf("\n\n");
	for (ic = 0; ic < n; ic++)
	{
		for (jc = 0; jc < n; jc++)
		{
			if (lattice[jc+ic*n]<0) printf("%i ",lattice[jc+ic*n]);
			else printf(" %i ",lattice[jc+ic*n]);
		}
	printf("\n");
	}
	printf("\n\n");

  return 0;
}
