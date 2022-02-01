#include <math.h>
#include <stdio.h>
#include <float.h>
#include "donut.h"
#include "poly34.h"

#ifndef WIN32
#include <sys/resource.h>
#include <sys/times.h>
#endif

double calculateExecutionTime(struct rusage *myTimeStart, struct rusage *myTimeEnd)
{
    return (
        ((double) (myTimeEnd->ru_utime.tv_sec  - myTimeStart->ru_utime.tv_sec)) +
        ((double) (myTimeEnd->ru_utime.tv_usec - myTimeStart->ru_utime.tv_usec)) * 1e-6
    );
}

#define PI      3.1415926535897932384626;
#define TWO_PI  6.2831853071795864769253;
#define HALF_PI 1.5707963267948966192313;

#define SQ(x) ((x)*(x))

// return vector ab
void vec(double *a, double *b, double *ret) {
	int i;
	for (i = 0; i < 3; i++) {
		ret[i] = b[i] - a[i];
	}
}

void sub(double *a, double *b, double *ret) {
	int i;
	for (i = 0; i < 3; i++) {
		ret[i] = a[i] - b[i];
	}
}

void add(double *a, double *b, double *ret) {
	int i;
	for (i = 0; i < 3; i++) {
		ret[i] = b[i] + a[i];
	}
}

void sca(double *p, double v, double *ret) {
	int i;
	for (i = 0; i < 3; i++) {
		ret[i] = p[i] * v;
	}
}

void opp_it(double *p) {
	int i;
	for (i = 0; i < 3; i++) {
		p[i] = -p[i];
	}
}

void cop(double *p, double *ret) {
	int i;
	for (i = 0; i < 3; i++) {
		ret[i] = p[i];
	}
}

double sq_mode(double *p) {
	return (SQ(p[0]) + SQ(p[1]) + SQ(p[2]));
}

double mode(double *p) {
	return sqrt(sq_mode(p));
}

double sq_dist(double *a, double *b) {
	double ret = 0.0;
	int i;
	for (i = 0; i < 3; i++) {
		ret += SQ(a[i] - b[i]);
	}
	return ret;
}

void unitify_it(double *p) {
	double norm = mode(p);
	int i;
	for (i = 0; i < 3; i++) {
		p[i] /= norm;
	}
}

void unitify(double *p, double *ret) {
	double norm = mode(p);
	int i;
	for (i = 0; i < 3; i++) {
		ret[i] = p[i] / norm;
	}
}

void cross_prd(double *a, double *b, double *ret) {
	ret[0] = a[1] * b[2] - a[2] * b[1];
	ret[1] = a[2] * b[0] - a[0] * b[2];
	ret[2] = a[0] * b[1] - a[1] * b[0];
}

double dot_prd(double *a, double *b) {
	double ret = 0.0;
	int i;
	for (i = 0; i < 3; i++) {
		ret += a[i] * b[i];
	}
	return ret;
}

double cos_theta(double *a, double *b) {
	return (dot_prd(a, b) / (mode(a) * mode(b)));
}

void get_perp_unit(double *p, double *ret) {
	if (p[0] != 0.0) {
		ret[0] = - (p[1] + p[2]) / p[0];
		ret[1] = 1.0;
		ret[2] = 1.0;
	} else if (p[1] != 0.0) {
		ret[0] = 1.0;
		ret[1] = - (p[0] + p[2]) / p[1];
		ret[2] = 1.0;
	} else if (p[2] != 0.0) {
		ret[0] = 1.0;
		ret[1] = 1.0;
		ret[2] = - (p[0] + p[1]) / p[2];
	}
	unitify_it(ret);
}

double donut_dist_sq(double *x, struct Circle *cir) {
    double sq_a = sq_dist(x, cir->o);
    double sq_n = sq_mode(cir->n);
    double lambda = (dot_prd(cir->o, cir->n) - dot_prd(x, cir->n)) / sq_n;
    double sq_h = SQ(lambda) * sq_n;
    return (sq_a + cir->sq_r - 2.0 * cir->r * sqrt(sq_a - sq_h));
}

void get_circle_pt(struct Circle *cir, double t, double *u1, double *u2, double *ret) {
    double cos_t = cos(t);
    double sin_t = sin(t);
    int i;
    for (i = 0; i < 3; i++) {
    	ret[i] = cir->o[i] + cir->r * u1[i] * cos_t + cir->r * u2[i] * sin_t;
    }
}

// e.g. when dim = 0, calculate the x bound [ca, cb] of circle, and store the points when x reaches ca and cb into ca_pt and cb_pt
void get_cir_bound(int dim, struct Circle *cir, double *ca, double *cb, double *ca_pt, double *cb_pt) {
	int i;
	if (cir->u1[dim] == 0.0) { // circle is parallel to an axis-aligned plane, thus the bound in this dim is 0-length
		get_circle_pt(cir, 0, cir->u1, cir->u2, ca_pt);
		*ca = ca_pt[dim];
		*cb = ca_pt[dim];
		for (i = 0; i < 3; i++) {
			cb_pt[i] = ca_pt[i];
		}
	} else {
		double opt_t1 = atan(cir->u2[dim] / cir->u1[dim]);
		double opt_t2 = opt_t1 + 3.1415926535897932384626;
		double opt_pt1[3], opt_pt2[3];
		get_circle_pt(cir, opt_t1, cir->u1, cir->u2, opt_pt1);
		get_circle_pt(cir, opt_t2, cir->u1, cir->u2, opt_pt2);
		if (opt_pt1[dim] > opt_pt2[dim]) {
			*ca = opt_pt2[dim];
			*cb = opt_pt1[dim];
			for (i = 0; i < 3; i++) {
				ca_pt[i] = opt_pt2[i];
				cb_pt[i] = opt_pt1[i];
			}
		} else {
			*ca = opt_pt1[dim];
			*cb = opt_pt2[dim];
			for (i = 0; i < 3; i++) {
				ca_pt[i] = opt_pt1[i];
				cb_pt[i] = opt_pt2[i];
			}
		}
	}
}

int contains(double *a, double *b, double *x, double *y) {
	int i;
	for (i = 0; i < 3; i++) {
		if (!(a[i] <= x[i] && b[i] >= y[i])) {
			return 0;
		}
	}
	return 1;
}

// min dist between the face (defined by e.g. when dim = 0 (x), x = val (0/20), y in [0, 50], z in [0, 30]) and the circle
double min_dist_face_int(struct Box *box, int dim, double val, struct Circle *cir) {
	int i, j;
	double delta, phi;
	double t_plus_phi[2];
	double t_inters;
	double pt_inters[3];
	int flag_inters_in_face, flag_proj_in_face;

	if (cir->ca[dim] <= val && val <= cir->cb[dim] && cir->ca[dim] != cir->cb[dim]) {
		// printf("    the bound of dim %d of circle covers val %f, then it is possible that the face and circle intersect\n", dim, val);

		/* calculate the intersect of the circle and plane
		 * let o[dim] + r cos(t) u1[dim] + r sin(t) u2[dim] = val
		 * then sin (t + phi) = delta defined as follows: */
	    delta = (val - cir->o[dim]) / (cir->r * sqrt(SQ(cir->u1[dim]) + SQ(cir->u2[dim])));
	    phi = atan2(cir->u1[dim], cir->u2[dim]);
	    // printf("    delta: %f, phi: %f\n", delta, phi);

	    if (-1 <= delta && delta <= 1) { // delta in [-1, 1], thus the intersect must exist, this must also be fulfilled in this case, since the dim-bound of circle covers val

	        t_plus_phi[0] = asin(delta);
	        t_plus_phi[1] = 3.1415926535897932384626 - t_plus_phi[0];

	        for (i = 0; i < 2; i++) { // the intersecting t might have two solutions, check them one-by-one

	        	t_inters = t_plus_phi[i] - phi;
	        	// printf("        t inters: %f\n", t_inters);
	        	get_circle_pt(cir, t_inters, cir->u1, cir->u2, pt_inters); // get the intersecting point coordinate
	        	// printf("    intersection at (%f, %f, %f)\n", pt_inters[0], pt_inters[1], pt_inters[2]);

	        	flag_inters_in_face = 1; // default true
	        	for (j = 0; j < 3; j++) {
	        		if (j == dim) { // no need to check dim x
	        			continue;
	        		}
	        		if (!(box->a[j] <= pt_inters[j] && pt_inters[j] <= box->b[j])) { // if e.g. intersecting point y = -10 not in [0, 50]
	        			flag_inters_in_face = 0;
	        			break;
	        		}
	        	}
		        if (flag_inters_in_face) {
		        	// printf("    intersection in face!\n");
		        	return 0;
		        }
	        }
	    }
	    // if not return 0.0, then it must be the case that the min dist appears on edges
	    return -1;
	} else {
		// printf("    the bound of dim %d of circle does not cover val, then we only consider the case when the ca_pts or cb_pts projects to the face\n", dim);
		if (val < cir->ca[dim]) { // face is on the left of circle, we calculate ca.x - face.x(val)
			flag_proj_in_face = 1; // default true
			for (i = 0; i < 3; i++) {
				if (i == dim) {
					continue;
				}
				if (!(box->a[i] <= cir->ca_pt[dim][i] && cir->ca_pt[dim][i] <= box->b[i])) { // if e.g. projection point y = -10 not in [0. 50]
					flag_proj_in_face = 0;
					break;
				}
			}
			if (flag_proj_in_face) {
		        // printf("    projection in face!\n");
				return (cir->ca[dim] - val);
			}
		} else { // face is on the right of circle, we calculate face.x(val) - cb.x
			flag_proj_in_face = 1; // default true
			for (i = 0; i < 3; i++) {
				if (i == dim) {
					continue;
				}
				if (!(box->a[i] <= cir->cb_pt[dim][i] && cir->cb_pt[dim][i] <= box->b[i])) { // if e.g. projection point y = -10 not in [0. 50]
					flag_proj_in_face = 0;
					break;
				}
			}
			if (flag_proj_in_face) {
		        // printf("    projection in face!\n");
				return (val - cir->cb[dim]);
			}
		}
		// if not return the projection dist, then it must be the case that the min dist appears on edges
		return -1;
	}
}

double cal_e_for(int dim, int d1, double k1, int d2, double k2, struct Circle *cir) {
	double o1_1, o2_2, o3, u1_1, u1_2, u1_3, u2_1, u2_2, u2_3;
	o1_1 = k1 - cir->o[d1]; o2_2 = k2 - cir->o[d2]; o3 = cir->o[dim];
	u1_1 = cir->u1[d1]; u1_2 = cir->u1[d2]; u1_3 = cir->u1[dim];
	u2_1 = cir->u2[d1]; u2_2 = cir->u2[d2]; u2_3 = cir->u2[dim];

	double tan_t = (o1_1 * u1_2 - o2_2 * u1_1) / (o2_2 * u2_1 - o1_1 * u2_2);
	double t = atan(tan_t);
	double cos_t = cos(t);
	double sin_t = sin(t);
	double lambda = o1_1 / (cir->r * (cos_t * u1_1 + sin_t * u2_1));
	return (o3 + lambda * cir->r * (cos_t * u1_3 + sin_t * u2_3));
}

int solve_A_by_cost_p_B_by_sint_eq_C(double A, double B, double C, double ret[4]) {
	int i;
	double t;

	double a = -2.0 * (A + C) / B;
	double b = 0.0;
	double c = -2.0 * (A - C) / B;
	double d = -1.0;
	double x[4];
	int num_root = SolveP4(x, a, b, c, d);
	for (i = 0; i < num_root; i++) {
		ret[i] = 2.0 * atan(x[i]);
	}
	return num_root;
}

double donut_dist_sq_by_t(double t, double S, double T, double tan_tau, double sq_tan_tau, double *g_pr, double *e, double *u2, struct Circle *cir, double *f) {
    double o_f_pr = S / cos(t);
    double R = o_f_pr * sin(t);
    double sq_e_f_pr = SQ(T + R);
    double sq_dist = sq_e_f_pr * sq_tan_tau + SQ(o_f_pr) - 2.0 * cir->r * o_f_pr + cir->sq_r;

    double u2_by_R[3], f_pr[3], e_f_pr[3], vert[3];
    sca(u2, R, u2_by_R);
    add(g_pr, u2_by_R, f_pr);
    vec(e, f_pr, e_f_pr);
    sca(cir->n_unit, mode(e_f_pr) * tan_tau, vert);
    if (dot_prd(e_f_pr, u2) > 0) {
    	add(f_pr, vert, f);
    } else {
    	sub(f_pr, vert, f);
    }

    return sq_dist;
}

// min dist between the edge (defined by e.g. when dim = 0 (x), x in [0, 20], y = 50, z = 0) and the circle
// in this case, d1 = 1, k1 = 50, d2 = 2, k2 = 0
double min_dist_edge(struct Box *box, int dim, int d1, double k1, int d2, double k2, struct Circle *cir) {
	double e[3], eo[3], d[3], d_pr[3];
	int d_up;
	double cos_tau_pr, tau_pr, tau;
	double h[3], h_pr[3], n_unit_sin_tau[3];
	double eo_cross_d_pr[3];
	double g_pr[3], e_g_pr[3], o_g_pr[3];
	double S, T;
	double new_u1[3], new_u2[3];
	double tan_tau, sq_tan_tau, A, B;
	double t_ret[4];
	int num_sol, i, j;
	double f[3];

	double sq_dist;
	double curr_min_sq = DBL_MAX;
	int opt_found = 0;

	e[dim] = cal_e_for(dim, d1, k1, d2, k2, cir);
	e[d1] = k1;
	e[d2] = k2;
    vec(e, cir->o, eo);

    d[dim] = 1;
    d[d1] = 0;
    d[d2] = 0;

    d_up = 1;
    if (dot_prd(d, cir->n) < 0) { // TODO: = 0? other special cases
    	d_up = 0;
    }

    cos_tau_pr = cos_theta(d, cir->n);
    tau_pr = acos(cos_tau_pr);
    tau = (tau_pr <= 1.5707963267948966192313) ? (1.5707963267948966192313 - tau_pr) : (tau_pr - 1.5707963267948966192313);

    add(e, d, h);
    sca(cir->n_unit, sin(tau), n_unit_sin_tau);
    if (d_up) {
    	sub(h, n_unit_sin_tau, h_pr);
    } else {
    	add(h, n_unit_sin_tau, h_pr);
    }
    vec(e, h_pr, d_pr);
    unitify_it(d_pr);

    // normalize such that d and d_pr point to right when o is below the line
    cross_prd(eo, d_pr, eo_cross_d_pr);
    if (dot_prd(eo_cross_d_pr, cir->n_unit) < 0) {
        opp_it(d);
        opp_it(d_pr);
        d_up = !d_up;
    }
    if (!d_up) {
        tau = -tau;
    }

    sca(d_pr, dot_prd(eo, d_pr), e_g_pr);
    add(e, e_g_pr, g_pr);
    vec(cir->o, g_pr, o_g_pr);
    S = mode(o_g_pr);

    cop(o_g_pr, new_u1);
    unitify_it(new_u1);
    cross_prd(new_u1, cir->n, new_u2);
    unitify_it(new_u2);

    T = mode(e_g_pr);
    if (dot_prd(e_g_pr, new_u2) < 0) {
        T = -T;
    }

    tan_tau = tan(tau);
    sq_tan_tau = SQ(tan_tau);
    A = S * (1.0 + sq_tan_tau);
    B = T * sq_tan_tau;

    // j = 1: positive connected; j = -1: negative connected
    for (j = 1; j >= -1; j -= 2) {
    	num_sol = solve_A_by_cost_p_B_by_sint_eq_C(A, B * j, cir->r, t_ret);
    	for (i = 0; i < num_sol; i++) {
	    	if (-1.5707963267948966192313 < t_ret[i] && t_ret[i] < 1.5707963267948966192313) {
	            sq_dist = donut_dist_sq_by_t(t_ret[i], S, T, tan_tau, sq_tan_tau, g_pr, e, new_u2, cir, f);

	            if (f[dim] < box->a[dim] || f[dim] > box->b[dim]) { // skip if the f point is not on the edge
	            	continue;
	            }

		    	if (sq_dist == 0.0) {
		    		return 0.0;
		    	} else if (sq_dist < curr_min_sq) {
		    		curr_min_sq = sq_dist;
		    	}

	            // printf("One optimal found %f\n", sqrt(sq_dist));
	            opt_found = 1;
	    	}
    	}
    }

    if (opt_found) {
    	return sqrt(curr_min_sq);
	} else {
		return -1;
	}

}

void pre_compute(struct Circle *cir) {
	cir->sq_r = SQ(cir->r);
	// C(t) = o + r cos(t) u1 + r sin(t) u2
    get_perp_unit(cir->n, cir->u1);
    cross_prd(cir->u1, cir->n, cir->u2);
    unitify_it(cir->u2);
    // unit vector of n
    unitify(cir->n, cir->n_unit);

    int i;
    // calculate the bound of the circle
    for (i = 0; i < 3; i++) {
    	// printf("bound for dim %d:\n", i);
    	get_cir_bound(i, cir, &cir->ca[i], &cir->cb[i], cir->ca_pt[i], cir->cb_pt[i]);
		// printf("bounded %f at (%f, %f, %f)\n", cir->ca[i], cir->ca_pt[i][0], cir->ca_pt[i][1], cir->ca_pt[i][2]);
		// printf("bounded %f at (%f, %f, %f)\n", cir->cb[i], cir->cb_pt[i][0], cir->cb_pt[i][1], cir->cb_pt[i][2]);
    }
}

double min_dist_donut(struct Box *box, struct Circle *cir, double *time1, double *time2, double *time3) {

	double curr_min = DBL_MAX; //, curr_max = 0.0;

    int i, j;

    // if fully contained return 0.0 as smallest
    if (contains(box->a, box->b, cir->ca, cir->cb)) {
    	// printf("fully contained\n");
    	return 0.0;
    }

	// struct rusage ru_start, ru_end;

	// getrusage(RUSAGE_SELF, &ru_start);

    // check each of the 6 faces
    for (i = 0; i < 3; i++) {
    	// printf("cal face a for dim %d, val = %f\n", i, box->a[i]);
    	double min_face_a = min_dist_face_int(box, i, box->a[i], cir);
    	// printf("face dist: %f\n", min_face_a);
    	if (min_face_a == 0.0) {
			// getrusage(RUSAGE_SELF, &ru_end);
			// *time1 += calculateExecutionTime(&ru_start, &ru_end) * 1000.0;
    		return 0.0;
    	} else if (0.0 < min_face_a && min_face_a < curr_min) {
    		curr_min = min_face_a;
    	}

    	// printf("cal face b for dim %d, val = %f\n", i, box->b[i]);
    	double min_face_b = min_dist_face_int(box, i, box->b[i], cir);
    	// printf("face dist: %f\n", min_face_b);
    	if (min_face_b == 0.0) {
			// getrusage(RUSAGE_SELF, &ru_end);
			// *time1 += calculateExecutionTime(&ru_start, &ru_end) * 1000.0;
    		return 0.0;
    	} else if (0.0 < min_face_b && min_face_b < curr_min) {
    		curr_min = min_face_b;
    	}
    }

	// getrusage(RUSAGE_SELF, &ru_end);
	// *time1 += calculateExecutionTime(&ru_start, &ru_end) * 1000.0;

	// getrusage(RUSAGE_SELF, &ru_start);

    // check each of the 12 edges
    for (i = 0; i < 3; i++) {
    	int d1, d2;
    	if (i == 0) {
    		d1 = 1; d2 = 2;
    	} else if (i == 1) {
    		d1 = 0; d2 = 2;
    	} else if (i == 2) {
    		d1 = 0; d2 = 1;
    	}

    	for (j = 0; j < 4; j++) {
    		double k1 = (j & 1) ? box->b[d1] : box->a[d1];
    		double k2 = (j & 2) ? box->b[d2] : box->a[d2];

	    	// printf("cal edge for dim %d, k[%d] = %f, k[%d] = %f\n", i, d1, k1, d2, k2);
	    	double min_edge = min_dist_edge(box, i, d1, k1, d2, k2, cir);
	    	// printf("edge dist: %f\n", min_edge);
	    	if (min_edge == 0.0) {
				// getrusage(RUSAGE_SELF, &ru_end);
				// *time2 += calculateExecutionTime(&ru_start, &ru_end) * 1000.0;
	    		return 0.0;
	    	} else if (0.0 < min_edge && min_edge < curr_min) {
	    		curr_min = min_edge;
	    	}
    	}
    }

	// getrusage(RUSAGE_SELF, &ru_end);
	// *time2 += calculateExecutionTime(&ru_start, &ru_end) * 1000.0;

	// getrusage(RUSAGE_SELF, &ru_start);

    // check each of the 8 vertices
    double ver[3];
    for (i = 0; i < 8; i++) {
    	ver[0] = (i & 1) ? box->b[0] : box->a[0];
    	ver[1] = (i & 2) ? box->b[1] : box->a[1];
    	ver[2] = (i & 4) ? box->b[2] : box->a[2];

    	// printf("cal vertex (%f, %f, %f)\n", ver[0], ver[1], ver[2]);
    	double min_vert = sqrt(donut_dist_sq(ver, cir));
    	// printf("vertex dist: %f\n", min_vert);
    	if (min_vert == 0.0) {
			// getrusage(RUSAGE_SELF, &ru_end);
			// *time3 += calculateExecutionTime(&ru_start, &ru_end) * 1000.0;
    		return 0.0;
    	} else if (min_vert < curr_min) {
    		curr_min = min_vert;
    	}
    }
    
	// getrusage(RUSAGE_SELF, &ru_end);
	// *time3 += calculateExecutionTime(&ru_start, &ru_end) * 1000.0;

    // printf("final min dist: %f\n", curr_min);

    return curr_min;

}