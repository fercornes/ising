#include "termo.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

double energia(int *lattice, int n, int j, int b, int seg) {
		
	int filas,columnas,e_inicial,act,aba,der,ie,ib,actual,diag_ab,diag_ar,diag_iz,diag_de;
	e_inicial=0;	
	
	if (j != 0)
	{
		//calculo la contribucion a la energía inicial mediante la interaccion spin-spin (utilizo bloques tipo P)
		for (filas = 0; filas < n-1; filas++)
		{
			for(columnas = 0; columnas < n-1; columnas++)
			{	
				act=filas*n+columnas;
				aba=act+n;
				der=act+1;	
				e_inicial=e_inicial-j*lattice[act]*lattice[aba]-j*lattice[act]*lattice[der];
			}
		}

		//calculo la interacción del spin ubicado en el extremo inferior derecho de la red
		e_inicial=e_inicial-j*lattice[act+n+1]*lattice[aba]-j*lattice[act+n+1]*lattice[der];

		//calculo la interaccion a segundos vecinos
		if (seg != 0)
		{
			for (filas = 0; filas < n; filas++)
			{
				for(columnas = 0; columnas < n; columnas++)
				{
					actual=filas*n+columnas;
					diag_ab=(filas-1+n)%n;
					diag_ar=(filas+1+n)%n;
					diag_iz=(columnas-1+n)%n;
					diag_de=(columnas+1+n)%n;
					e_inicial=e_inicial+j*lattice[actual]*(lattice[(filas*n+diag_iz+(n-1)*n)%(n*n)]+lattice[(filas*n+diag_de+(n-1)*n)%(n*n)]+lattice[(filas*n+diag_iz+n)%(n*n)]+lattice[(filas*n+diag_de+n)%(n*n)]);
	
					//printf("%i\t%i\t%i\t%i\n",(filas*n+diag_iz+(n-1)*n)%(n*n),(filas*n+diag_de+(n-1)*n)%(n*n),(filas*n+diag_iz+n)%(n*n),(filas*n+diag_de+n)%(n*n));
				}
			}
		}

		//Condiciones periódicas
		for (ie = 0; ie < n; ie++) 
		{
			e_inicial=e_inicial-j*lattice[ie]*lattice[(n-1)*n+ie]-j*lattice[n*ie]*lattice[n*ie+n-1];
		}
		//printf("e_inicial=%i\n",e_inicial);
	}

	//calculo la contribución a la energía inicial si hay campo externo
	if (b != 0) 
	{
		for (ib = 0; ib < n*n; ib++) e_inicial=e_inicial-b*lattice[ib];
	}

   return (double)e_inicial/(double)(n*n);
}

double magnetizacion(int *lattice, int n) {
	int iim,magnet;
	magnet=0;	
	
	for (iim = 0; iim < n*n; iim++) magnet=magnet+lattice[iim];

   return (double)magnet/(double)(n*n);
}

double func_resp(double *vec_termo, int ns, int n, double *termo_e_m) {
	int ii;
	double energia_media,energia_media_cuadr,magnet_media,magnet_media_cuadr;
	energia_media = 0.0;
	energia_media_cuadr = 0.0;
	magnet_media = 0.0;
	magnet_media_cuadr = 0.0;

	for (ii = 0; ii < ns; ii++) 
	{
		energia_media=energia_media+vec_termo[ii];
		energia_media_cuadr=energia_media_cuadr+vec_termo[ii]*vec_termo[ii];
		magnet_media=magnet_media+vec_termo[ii+ns];
		magnet_media_cuadr=magnet_media_cuadr+vec_termo[ii+ns]*vec_termo[ii+ns];
	}
	termo_e_m[0]=energia_media/(double)ns;
	termo_e_m[1]=energia_media_cuadr/(double)ns;
	termo_e_m[2]=magnet_media/(double)ns;
	termo_e_m[3]=magnet_media_cuadr/(double)ns;

	//if (termo_e_m[0]>5 || termo_e_m[1]>5 ||termo_e_m[2]>5 || termo_e_m[3]>5) printf("energia_media=%.6f\tmagnet_media=%.6f\n",energia_media,magnet_media);

	//printf("%g\t%g\t%g\t%g\n",termo_e_m[0],termo_e_m[1],termo_e_m[2],termo_e_m[3]);
	
   return 0;
}


double correlacion(double *vec_termo, int nsamp, int n_vec_corr, double *vec_corr, int ns) {
	int ic,k,ki,iiii;
	double e_media,m_media,num_pk_e,num_pk_m,den_pk_e,den_pk_m;
	e_media = 0.0;
	m_media = 0.0;
	den_pk_e = 0.0;
	den_pk_m = 0.0;
		
	for (ic = 0; ic < nsamp; ic++) 
	{	
		e_media=e_media+vec_termo[ic];
		m_media=m_media+vec_termo[ic+nsamp];
	}
	
	e_media=e_media/(double)nsamp;
	m_media=m_media/(double)nsamp;

	for (iiii = 0; iiii < nsamp; iiii++) 
	{
		den_pk_e=den_pk_e+(vec_termo[iiii]-e_media)*(vec_termo[iiii]-e_media);	
		den_pk_m=den_pk_m+(vec_termo[iiii+nsamp]-m_media)*(vec_termo[iiii+nsamp]-m_media);		
	}

	for (k = 0; k < n_vec_corr; k++)
	{
		num_pk_e = 0.0;
		num_pk_m = 0.0;
		for (ki = 0; ki < nsamp-k; ki++)
		{
			num_pk_e=num_pk_e+(vec_termo[ki]-e_media)*(vec_termo[ki+k]-e_media);
			num_pk_m=num_pk_m+(vec_termo[ki+nsamp]-m_media)*(vec_termo[ki+nsamp+k]-m_media);
		}
		vec_corr[k]=vec_corr[k]+num_pk_e/(double)(den_pk_e*ns);
		vec_corr[k+n_vec_corr]=vec_corr[k+n_vec_corr]+num_pk_m/(double)(den_pk_m*ns);
	}
	//printf("vec_corr[0]=%g\n",vec_corr[0]);

   return 0;
}

