
/*
 * Include Files
 *
 */
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#define SIMPLIFIED_RTWTYPES_COMPATIBILITY
#include "rtwtypes.h"
#undef SIMPLIFIED_RTWTYPES_COMPATIBILITY
#endif



/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
#include <math.h>
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 4
#define u_1_width 4
#define u_2_width 4
#define u_3_width 16
#define y_width 4

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* extern double func(double a); */
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output function
 *
 */
void MPC_Params_Outputs_wrapper(const real_T *x_,
			const real_T *p_max,
			const real_T *v_,
			const real_T *gbus,
			real_T *u,
			const real_T *h_, const int_T p_width0,
			const real_T *n_, const int_T p_width1,
			const real_T *ni_, const int_T p_width2,
			const real_T *paso_, const int_T p_width3,
			const real_T *pv_, const int_T p_width4,
			const real_T *Vbase_, const int_T p_width5,
			const real_T *Pbase_, const int_T p_width6,
			const real_T *T1, const int_T p_width7,
			const real_T *T2, const int_T p_width8,
			const real_T *T3, const int_T p_width9,
			const real_T *T4, const int_T p_width10,
			const real_T *z1, const int_T p_width11,
			const real_T *z2, const int_T p_width12,
			const real_T *z3, const int_T p_width13,
			const real_T *z4, const int_T p_width14,
			const real_T *px_, const int_T p_width15)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
int a = 0, k = 0, kk = 0;
	//
    double h = *(h_); // tiempo de discretizacion
	int n = *n_;
    // Parametros del control
	int ni = *ni_;
	double paso = *paso_;
	double pv = *pv_; //peso de funcion de voltaje
    double px = *px_;
    double Vbase = *Vbase_; // Voltaje Base
	double Pbase = *Pbase_; // Potencia Base
	double *tau; tau = calloc(n, sizeof(double)); // parametro tau
	double *droop; droop = calloc(n, sizeof(double)); // parametro droop
    *(tau+0) = *T1;
    *(tau+1) = *T2;
    *(tau+2) = *T3;
    *(tau+3) = *T4;
    *(droop+0) = *z1*(Pbase/Vbase);
    *(droop+1) = *z2*(Pbase/Vbase);
    *(droop+2) = *z3*(Pbase/Vbase);
    *(droop+3) = *z4*(Pbase/Vbase);
	double *M; M = calloc(n*n, sizeof(double)); // declarando M
	double *Q; Q = calloc(n*n, sizeof(double)); // declarando Q
	double *x; x = calloc(n, sizeof(double));
	double *pmax; pmax = calloc(n, sizeof(double));
    double *v; v = calloc(n, sizeof(double));
    double *u_; u_ = calloc(n, sizeof(double)); // variable de control en pu
    double *y; y = calloc(n*n, sizeof(double));	
    double *aux0; aux0 = calloc(n, sizeof(double)); // vectores
    double *aux1; aux1 = calloc(n, sizeof(double)); 
    double *aux2; aux2 = calloc(n, sizeof(double));
    double *aux5; aux5 = calloc(n, sizeof(double));  
	double *aux6; aux6 = calloc(n, sizeof(double));  
    double *aux3; aux3 = calloc(n*n, sizeof(double)); // matriz
    double *aux4; aux4 = calloc(n*n, sizeof(double));
    double *v_nodal; v_nodal = calloc(n, sizeof(double));
    double *dx; dx = calloc(n, sizeof(double));
    double *dv; dv = calloc(n, sizeof(double));
    double *c; c = calloc(n, sizeof(double));
    double *b; b = calloc(n, sizeof(double));
    double ones[4] = {1,1,1,1};
    // Toma de medidas
	for(a = 0; a < n; a++)
	{
		// Convirtiendo en pu
		*(pmax+a) = *(p_max+a)/Pbase;
		*(x+a) = *(x_+a)/Pbase;
		*(v+a) = *(v_+a)/Vbase;
	}
    //Creando M
	for(a = 0; a<n; a++)
	{
		*(aux1+a) = 1/ *(pmax+a);
		*(aux2+a) = *(aux1+a)/n; 
	}
    // diag(1./pmax)
    int i, r, rr = 0;
    for(i = 0; i<n; i++)
	{
		for(r = 0; r<n; r++)
		{
			if(i == r)
			{
				*(aux3+rr) = *(aux1+i);
			}
			else
			{
				*(aux3+rr) = 0;
			}
			rr++;
		}
	}
    rr = 0;
    // kron(1./pmax'/n,ones(n,1))
    // aux4 = kron(aux2,ones(n),n);
	for(i = 0; i<n; i++)
	{
		for(r = 0; r<n; r++)
		{
			*(aux4+rr) = *(aux2+i)**(ones+r);
			rr++;
		}
	}
	rr = 0;
    // M = diag(1./pmax)-kron(1./pmax'/n,ones(n,1));
    for(i = 0; i < n*n; i++)
    {
        *(M+i) = *(aux3+i)- *(aux4+i);
    }
    double *Ma; Ma = calloc(n*n, sizeof(double)); // declarando M
	i=0, r=0, a = 0;
	for(i = 0; i<n; i++)
	{
		for(r = 0; r<n; r++)
		{
			*(Ma+a) = *(M+i+n*r);
			a++;
		}
	}    
    // Q = M'*M
    int j=0;
    a = 0, i = 0, k = 0;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			*(y+i+j) = 0; //inicializar
			for(k = 0; k < n; k++)
			{
				*(y+i+j) = *(y+i+j) + *(M+i*n+k)**(Ma+j+k*n);
			}
			*(Q+a) = *(y+i+j);
			a = a+1;
		}
	}
    // u_ = ones(n,1);
    for(i=0; i<n; i++)
    {
        *(u_+i) = *(ones+i);
    }
	// Modelo de la microrred
    for(a = 0; a < n; a++)
	{
		*(v_nodal+a) = *(u_+a) - *(droop+a)*(*(x+a)- *(pmax+a)); // v_nodal = u_ - (droop).*(x-pmax)
	}
    // gbus*v_nodal
    j = 0, a = 0, i = 0, k = 0;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < 1; j++)
		{
			*(y+i+j) = 0; //inicializar
			for(k = 0; k < n; k++)
			{
				*(y+i+j) = *(y+i+j) + *(gbus+i*n+k)**(v_nodal+j+k*1);
			}
			*(aux1+a) = *(y+i+j);
			a = a+1;
		}
	}
	for(a = 0; a < n; a++)
	{
		*(dx+a) = (1/ *(tau+a))*((*(v_nodal+a)**(aux1+a))- *(x+a)); // dx = 1./tau.*(v_nodal.*gbus*v_nodal-x);
		*(x+a)  = *(x+a) + *(dx+a)*h; // x = x + dx*h
	// Modelo de Optimizacion
		*(c+a)  = (1-(h/ *(tau+a)))**(x+a); // c = (1-h./tau).*x
		*(b+a)  = h/ *(tau+a); // b = h./tau
	}
    for(kk = 0; kk < ni; kk++)
    {
        // t0 = gbus*v
        j = 0, a = 0, i = 0, k = 0;
	    for(i = 0; i < n; i++)
	    {
		    for(j = 0; j < 1; j++)
		    {
			    *(y+i+j) = 0; //inicializar
			    for(k = 0; k < n; k++)
			    {
				    *(y+i+j) = *(y+i+j) + *(gbus+i*n+k)**(v+j+k*1);
			    }
			    *(aux0+a) = *(y+i+j);
			    a = a+1;
		    }
	    }
        for(a = 0; a < n; a++)
		{	
			*(aux1+a) = *(b+a)**(v+a); // t1 = b.*v
			*(aux6+a) = *(c+a) + *(aux1+a)**(aux0+a); // (c+t1.*t0)
		}
        // t2 = Q*(c+t1.*t0)
        j = 0, a = 0, i = 0, k = 0;
	    for(i = 0; i < n; i++)
	    {
		    for(j = 0; j < 1; j++)
		    {
			    *(y+i+j) = 0; //inicializar
			    for(k = 0; k < n; k++)
			    {
				    *(y+i+j) = *(y+i+j) + *(Q+i*n+k)**(aux6+j+k*1);
			    }
			    *(aux2+a) = *(y+i+j);
			    a = a+1;
		    }
	    }
		for(a = 0; a < n; a++)
		{
			*(aux6+a) = *(b+a)*(*(aux2+a)**(aux0+a)); // b.*(t2.*t0)
			*(aux3+a) = pv*(*(v+a)-1); // pv*(v-1)
			*(aux4+a) = *(aux2+a)**(aux1+a); // (t2.*t1)
		} 
        // gbus*(t2.*t1)
        j = 0, a = 0, i = 0, k = 0;
	    for(i = 0; i < n; i++)
	    {
		    for(j = 0; j < 1; j++)
		    {
			    *(y+i+j) = 0; //inicializar
			    for(k = 0; k < n; k++)
			    {
				    *(y+i+j) = *(y+i+j) + *(gbus+i*n+k)**(aux4+j+k*1);
			    }
			    *(aux5+a) = *(y+i+j);
			    a = a+1;
		    }
	    }
        for(a=0; a<n; a++)
		{
			*(dv+a) = px*(*(aux6+a)+ *(aux5+a))+ *(aux3+a); // dv = b.*(t2.*t0)+gbus*(t2.*t1)+pv*(v-1)
			*(v+a) = *(v+a) - paso**(dv+a);
		}
    }
    for(a = 0; a < n; a++)
    {
        *(u_+a) = *(v+a) + *(droop+a)*(*(x+a)- *(pmax+a));
		*(u+a) = *(u_+a)*Vbase;
        //*(u+a) = *(v_nodal+a);
    }
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


