#include "stdlib.h"
#include <stdio.h>
#include "time.h"
#include <math.h>

#include "metropolis.h"
#include "lattice.h"
#include "termo.h"

void escribir(double *termo_e_m);

int main(int argc, char *argv[]) {

  int n,b,j,nsamp,ns,nterm,nt,seg;
  float T,prob;
  double energ,m,jprima;
 
  n = 16;
  b = 0;
  j = -1;
  seg = 0;  
  nsamp = 2*n*n;
  ns = 10000;
  nterm=10*n*n;
  //nterm=6;
  nt=1000;
  prob = 0.75;	

  int *lattice = malloc(n * n * sizeof(int));
  double *vec_termo = malloc(2 * ns * sizeof(double));
  
  double vec_exp[5];
  //double vec_exp[9];
  double termo_e_m[4];
  //double vec_termo[2*ns];

  
  srand(time(NULL));
  fill_lattice(lattice, n, prob);

  /* if (argc==3)
  {
  	sscanf(argv[1],"%d",&b); //primer argumento: campo magnetico
	sscanf(argv[2],"%d",&j); //segundo argumento: acomplamiento
  }
  */

  for (int t = 0; t < nt; t++) //hago un barrido en T
  {
	T=5.0*(1.0-t/(float)nt);
	jprima=(double)j/(double)T;

  	for (int i = 0; i < 5; i++) vec_exp[i]=exp((float)(-j*(4*i-8)-2*b)/T); //dE=Ef-Ei=j*(4*i-8)+2*b
	//for (int i = 0; i < 9; i++) vec_exp[i]=exp((float)(-j*(4*i-8)*(2*s)-2*b)/T); //dE=Ef-Ei=j*(4*i-8)+2*b
	
	//Calculo la energia y la magnetizacion inicial (por sitio) de la red
 	energ=energia(lattice,n,j,b,seg);
	m=magnetizacion(lattice,n);
	//printf("ener=%g\tmagnet=%g\n",energ,m);

	if (t>0) nterm = n*n;

  	//primero dejo que "termalice el sistema" 
  	for (int ij = 0; ij < nterm; ij++) {
		//print_lattice(lattice, n);
		metropolis(lattice, n, j, b, seg, vec_exp, &energ, &m);
  	}

  	//ahora sampleo despues de nsamp pasos un total de ns veces
  	for (int s = 0; s < ns; s++)
  	{
		//for (int iii = 0; iii < 4; iii++) termo_e_m[iii]=0.0;		
		
		for (int is = 0; is < nsamp; is++) 
		{
			metropolis(lattice, n, j, b, seg, vec_exp, &energ, &m);
		}
		vec_termo[s]=energ;
		vec_termo[s+ns]=m;
  	}
	func_resp(vec_termo,ns,n,termo_e_m);
	escribir(termo_e_m);
	if (t==0 || t==200 || t==400 || t==600 || t==800) print_lattice(lattice, n);

  }
  //print_lattice(lattice, n);

  free (lattice);
  free (vec_termo);
  return 0;
}

void escribir(double *termo_e_m){
	int ip;	
	FILE *fp;
	fp=fopen("termo_n_32_j_1_b_0_prueba.txt","a");
	//fp=fopen("termo_n_16_j_1_b_0_s_1.txt","a");
	for(ip = 0; ip < 4; ip++) fprintf(fp,"%.6f\t",termo_e_m[ip]);
	fprintf(fp,"\n");
	fclose(fp);
}
