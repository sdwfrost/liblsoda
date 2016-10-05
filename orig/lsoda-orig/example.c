#include <stdio.h>
#include <math.h>

static void *fex (neq, t, y, ydot) double t, *y, *ydot; int neq;
{
ydot[1] = 1.0E4*y[2]*y[3] - .04E0*y[1] ;
ydot[3] = 3.0E7*y[2]*y[2];
ydot[2] = -1.0 * (ydot[1] + ydot[3]);
}

extern void lsoda();

main()
{
double rwork1, rwork5, rwork6, rwork7;
double atol[4], rtol[4], t, tout, y[4];
int iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9;
int neq = 3;
int itol, itask, istate, iopt, jt, iout, j, jdum;

iwork1= iwork2= iwork5= iwork6= iwork7= iwork8= iwork9= 0;
rwork1= rwork5= rwork6= rwork7= 0.0;
y[1] = 1.0E0;
y[2] = 0.0E0;
y[3] = 0.0E0;
t = 0.0E0;
tout = 0.4E0;
itol = 2;
rtol[0] = 0.0; atol[0] = 0.0;
rtol[1] = rtol[3] = 1.0E-4;  rtol[2]= 1.0E-8;
atol[1] = 1.0E-6;
atol[2] = 1.0E-10;
atol[3] = 1.0E-6;
itask = 1;
istate = 1;
iopt = 0;
jt = 2;

for (iout=1; iout <= 12; iout++) {
  lsoda(fex,neq,y,&t,tout,itol,rtol,atol,itask,&istate,iopt,jt,
	iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9,
	rwork1,rwork5,rwork6,rwork7);

	printf (" at t= %12.4e y= %14.6e %14.6e %14.6e\n", t,y[1],y[2],y[3]);
  if (istate <= 0) {
		printf("error istate = %d\n",istate);
		exit(0);
		}
	tout = tout*10.0E0;
	}
exit(1);
}
