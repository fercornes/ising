#include "metropolis.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int metropolis(int *lattice, int n, int j, int b, int seg, double *vec_exp, double *energ, double *m) {
   	
	//int im;
	//int indice;
	//indice = 0;
	//for (im = 0; im < 2; im++) indices[im] = 0;	
	//pick_site(lattice,n,&indice);

	int indices[2];
	pick_site(lattice,n,indices);
	flip(lattice,n,j,b,seg,indices,vec_exp,energ,m);

   return 0;
}

int pick_site(int *lattice, int n, int *indices) {
	
	for (int ia = 0; ia < 2; ia++) indices[ia] = rand()%n;

   return 0;
}

int flip(int *lattice, int n, int j, int b, int seg, int *indices, double *vec_exp, double *energ, double *m) {
	
	//calculo la energía con sus vecinos sin "flipear" ("sumo" las orientaciones de los vecinos)
	int energia_inic,dE,dM,fila,columna,actual,abajo,arriba,derecha,izquierda,diag_ab,diag_ar,diag_iz,diag_de;
	double p,es;	
	
	//la sumatoria la organizo así: izquierda derecha abajo arriba
	
	energia_inic = 0;
	fila=indices[0];
	columna=indices[1];
	actual=fila*n+columna;
	abajo=(fila-1+n)%n;
	arriba=(fila+1+n)%n;
	izquierda=(columna-1+n)%n;
	derecha=(columna+1+n)%n;

	if (seg != 0) //ingresa solamente si considero la interacción a segundos vecinos
	{
		diag_ab=(fila-1+n)%n;
		diag_ar=(fila+1+n)%n;
		diag_iz=(columna-1+n)%n;
		diag_de=(columna+1+n)%n;
		energia_inic=j*lattice[actual]*(lattice[(fila*n+diag_iz+(n-1)*n)%(n*n)]+lattice[(fila*n+diag_de+(n-1)*n)%(n*n)]+lattice[(fila*n+diag_iz+n)%(n*n)]+lattice[(fila*n+diag_de+n)%(n*n)]);
	}

	//energía inicial del sistemita de 5 espines sin flipear
	energia_inic=energia_inic-j*lattice[actual]*(lattice[fila*n+izquierda]+lattice[fila*n+derecha]+lattice[arriba*n+columna]+lattice[abajo*n+columna]);	

	//si "flipeo" el spin central la energía final es -energía
	dE=-2*energia_inic-2*b*lattice[actual]; 
	//printf("dE=%i\n",dE);	
	
	dM=lattice[actual];
	
	if (dE < 0) 
	{
		lattice[actual]=-lattice[actual];	
		dM=lattice[actual]-dM;
	}
	else //en ising.c dE=j*(-8,-4,0,4,8)
	{
		if (b == 0 && seg == 0 && j>0) p=vec_exp[dE/4+2];
		else if (b == 0 && seg == 0 && j<0) p=vec_exp[2-dE/4];
		//else if (b == 0 && seg != 0) p=vec_exp[dE/4+4];
		//else if (j == 0) p=vec_exp[0]; //dE>0

		es = (double)rand()/(double)RAND_MAX;

		//printf("p=%g\tes=%g\n",p,es);
		if (es <= p) 
		{
			lattice[actual]=-lattice[actual];
			dM=lattice[actual]-dM;
		}
		else 
		{
			dE=0;
			dM=0;
		}
	}
	//printf("dE=%i\tdM=%i\n",dE,dM);

	*energ=*energ+(double)dE/(double)(n*n);
	*m=*m+(double)dM/(double)(n*n);	
	//printf("energ=%g\tmag=%g\n",*energ,*m);

  return 0;
}
