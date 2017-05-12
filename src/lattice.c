#include "lattice.h"
#include <stdio.h>
#include <stdlib.h>

int fill_lattice(int *lattice, int n, float p) {
	int i;
	for (i = 0; i < n*n; i++){
		double e;
		e = (double)rand()/(double)RAND_MAX;
		if (e <= p){
			lattice[i]=1;
		}
		else {
			lattice[i]=-1;
		}
	}
 return 0;
}

int print_lattice(int *lattice, int n) {
	int i,j;
	printf("\n\n");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (lattice[j+i*n]<0) printf("%i ",lattice[j+i*n]);
			else printf(" %i ",lattice[j+i*n]);
		}
	printf("\n");
	}
	printf("\n\n");

  return 0;
}
