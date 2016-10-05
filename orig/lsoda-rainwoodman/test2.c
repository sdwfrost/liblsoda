#include <stdio.h>
#include "common.h"
#include "lsoda.h"

int function(double t, double * y, double * ydot, void * data) {
	double gamma_HI = 6067.4902268982732;
	double alpha_HII_A = 4966736.4918453256;
	double ion = 0.13580749322316402;
	double alpha_HII_AB = 2441214.829750936;
	double xHII = 1.0 - y[0];
	double xHI = y[0];
	double ye = xHII;
	ydot[0] = -ion + (- gamma_HI * ye * xHI + alpha_HII_A * ye * xHII);
	return 1;
}

void main(void) {
	int i;
		test();
}
int test(void) {
	struct lsoda_context_t ctx = {
		.function = function,
		.neq = 1,
		.data = NULL,
		.state = 1,
	};
	static double rtol[7] = {1e-7, 1e-7, 1e-7, 1e-7, 1e-7, 1e-7, 1e-7};
	static double atol[7] = {1e-7, 1e-7, 1e-7, 1e-7, 1e-7, 1e-7, 1.};

	static struct lsoda_opt_t opt = {
		.ixpr = 1,
		.rtol = rtol,
		.atol = atol,
		.itask = 1,
		.ixpr = 1,
	};
	double x[1];
	x[0] = 1.0;
	double step = 7.3633639519195151e-06;
	int i;
	int k;
	for(k = 0; k < 1000; k++) {
		lsoda_prepare(&ctx, &opt);

		for(i = 0; i < 1001; i ++) {
			double t = 0;
			lsoda_reset(&ctx);
			ctx.state = 1;
			lsoda(&ctx, x, &t, step);
		}
		//printf ("%g %g\n", i * step, x[0]);
		lsoda_free(&ctx);
	}
}
