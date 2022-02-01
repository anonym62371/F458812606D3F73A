#ifndef __donut_h
#define __donut_h

struct Circle {
	double o[3];
	double r, sq_r;
	double n[3], n_unit[3];
	double u1[3], u2[3];
	double ca[3], cb[3];
	double ca_pt[3][3]; // ca_pt[0] is the point when circle reaches smallest x
	double cb_pt[3][3];
};

struct Box {
	double a[3], b[3];
};

void pre_compute(struct Circle *cir);

// return the min dist from a 3D box (defined by a[3], b[3]) to a 3D circle (defined by o[3], n[3], r)
// a, b, o, n, r must be filled before calling this function
double min_dist_donut(struct Box *box, struct Circle *cir, double *time1, double *time2, double *time3);

double donut_dist_sq(double *x, struct Circle *cir);

#endif