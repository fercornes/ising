#include "metropolis.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int metropolis(int *lattice, int n, float T) {
   
//calculo la energía inicial (utilizo bloques tipo P)	
	int fila,columna,e_inicial,actual,abajo,derecha,i,indice;
	e_inicial=0;	
	for (fila=0;fila<n-1;fila++)
	{
		for(columna=0;columna<n-1;columna++)
		{	
		actual=fila*n+columna;
		abajo=actual+n;
		derecha=actual+1;	
		e_inicial=e_inicial+lattice[actual]*lattice[abajo]+lattice[actual]*lattice[derecha];
		}
	}
	
	//calculo la interacción del spin ubicado en el extremo inferior derecho de la red
	e_inicial=e_inicial+lattice[actual+n+1]*lattice[abajo]+lattice[actual+n+1]*lattice[derecha];

	//Condiciones periódicas
	for (i=0;i<n;i++) 
	{
		e_inicial=e_inicial+lattice[i]*lattice[(n-1)*n+i]+lattice[n*i]*lattice[n*i+n-1];
	}

	indice=pick_site(lattice,n);
	flip(lattice,n,T,indice);
	//printf("energia=%i\n",e_inicial);

   return 0;
}

int pick_site(int *lattice, int n) {
	double e;
	e = n*n*(double)rand()/(double)RAND_MAX;
	//printf("%g\n",e);
  return floor(e);
}

int flip(int *lattice, int n, float T, int indice) {
	//calculo la energía con sus vecinos sin "flipear" ("sumo" las orientaciones de los vecinos)
	int energia,ei;
	float prob;
	energia = 0;
	ei=0;	
	prob=0.0;

	//la sumatoria la organizo así: izquierda derecha abajo arriba
	if (indice<n) //primera fila
	{
		if (indice==0) ei=lattice[indice+n-1]+lattice[indice+1]+lattice[indice+n]+lattice[indice+n*(n-1)];
		else if (indice==n-1) ei=lattice[indice-1]+lattice[indice-n+1]+lattice[indice+n]+lattice[indice+n*(n-1)];
		else ei=lattice[indice-1]+lattice[indice+1]+lattice[indice+n]+lattice[indice+n*(n-1)];
	}
	else if (indice>=n*(n-1)) //última fila
	{
		if (indice==n*(n-1)) ei=lattice[indice+n-1]+lattice[indice+1]+lattice[indice-n*(n-1)]+lattice[indice-n];
		else if (indice==n*n-1) ei=lattice[indice-1]+lattice[indice-n+1]+lattice[indice-n*(n-1)]+lattice[indice-n];
		else ei=lattice[indice-1]+lattice[indice+1]+lattice[indice-n*(n-1)]+lattice[indice-n];
	}
	else  //filas interiores
	{
		if (indice % n == 0) ei=lattice[indice+n-1]+lattice[indice+1]+lattice[indice+n]+lattice[indice-n];
		else if ((indice+1) % n == 0) ei=lattice[indice-1]+lattice[indice-n+1]+lattice[indice+n]+lattice[indice-n];
		else ei=lattice[indice-1]+lattice[indice+1]+lattice[indice+n]+lattice[indice-n];
	}
	
	energia=ei*lattice[indice]; //energía sin flipear
	//printf("indice=%i\tei=%i\tenergia=%i\n",indice,ei,energia);
	
	//la energía si "flipeo" es -energia. Así que me fijo si "energía" (sin flipear) es positiva o negativa. Si es positiva entonces cuando flipeo llego a un estado de menor energía. Y viceversa. 

	if (energia>0) lattice[indice]=-lattice[indice]; //si flipeo llego a un estado de menor energia, entonces lo "acepto"
	else 
	{
		prob=exp(2*energia/T); //en este caso Ef>=Ei ya que energia<=0
		double e;
		e = (double)rand()/(double)RAND_MAX;
		//printf("%g\n",e);
		if (e<=prob) lattice[indice]=-lattice[indice];
		else lattice[indice]=lattice[indice];
	}

  return 0;
}
