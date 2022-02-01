#ifndef __UTILPOLYGON_H_
#define __UTILPOLYGON_H_

#include "TriMesh.h"
#include "Vec.h"

#include "UtilVec.hpp"
#include "UtilMatrix.hpp"

// template<int S>
// using Polygon = trimesh::point[S];
// template<int S>
// using PolygonPtr = trimesh::point(*)[S];

inline void normal(const trimesh::point (&t)[3], trimesh::vec & ret)
{
    ret = ((t[1] - t[0]) CROSS (t[2] - t[0]));
    normalize(ret);
}

inline trimesh::vec normal(const trimesh::point (&t)[3])
{
    trimesh::vec ret;
    normal(t, ret);
    return ret;
}

class TriSide
{
private:
    const trimesh::point *m_v[2]; // s, t
    trimesh::vec m_e; // st = t - s
    float m_l2; // dist2(s, t)
    void loadEdge()
    {
        this->m_e = (*this->m_v[1]) - (*this->m_v[0]);
        this->m_l2 = trimesh::len2(m_e);
    }
public:
    TriSide& set(const trimesh::point *s, const trimesh::point *t)
    {
        this->m_v[0] = s;
        this->m_v[1] = t;
        this->loadEdge();
        return (*this);
    }
    TriSide()
    {
        this->m_v[0] = nullptr;
        this->m_v[1] = nullptr;
    }
    TriSide(int ep1, int ep2, trimesh::TriMesh* mesh)
    {
        this->set(&(mesh->vertices[ep1]), &(mesh->vertices[ep2]));
    }
    TriSide(const trimesh::point *s, const trimesh::point *t)
    {
        this->set(s, t);
    }
    TriSide(const trimesh::point (&p)[2])
    {
        this->set(&(p[0]), &(p[1]));
    }
    const trimesh::point& s() const
    {
        return (*this->m_v[0]);
    }
    const trimesh::point& t() const
    {
        return (*this->m_v[1]);
    }
    const trimesh::vec& st() const
    {
        return this->m_e;
    }
    float l2() const
    {
        return this->m_l2;
    }
};

class Triangle
{
private:
    const trimesh::point *m_v[3]; // a, b, c
    trimesh::vec m_n;
    TriSide m_s[3]; // ab, ac, bc
    void loadSides()
    {
        const int comb[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
        for (int i = 0; i < 3; i++) {
            const int i1 = comb[i][0];
            const int i2 = comb[i][1];
            this->m_s[i].set(this->m_v[i1], this->m_v[i2]); // ab = b - a, ac = c - a, bc = c - b
        }
    }
public:
    Triangle()
    {
        this->m_v[0] = nullptr;
        this->m_v[1] = nullptr;
        this->m_v[2] = nullptr;
    }
    Triangle(int f, trimesh::TriMesh* mesh)
    {
        this->m_v[0] = &(mesh->vertices[mesh->faces[f][0]]);
        this->m_v[1] = &(mesh->vertices[mesh->faces[f][1]]);
        this->m_v[2] = &(mesh->vertices[mesh->faces[f][2]]);
        this->m_n = mesh->trinorm(f);
        this->loadSides();
    }
    Triangle(const trimesh::point (&tv)[3])
    {
        this->m_v[0] = &(tv[0]);
        this->m_v[1] = &(tv[1]);
        this->m_v[2] = &(tv[2]);
        normal(tv, this->m_n);
        this->loadSides();
    }
    Triangle(const trimesh::point *_v[3])
    {
        this->m_v[0] = _v[0];
        this->m_v[1] = _v[1];
        this->m_v[2] = _v[2];
        normal({*(_v[0]), *(_v[1]), *(_v[2])}, this->m_n);
        this->loadSides();
    }
    const trimesh::point& v(int i) const
    {
        return (*this->m_v[i]);
    }
    const trimesh::vec& n() const
    {
        return this->m_n;
    }
    const TriSide& s(int i) const
    {
        return this->m_s[i];
    }
    friend std::ostream & operator<<(std::ostream & os, const Triangle & t)
    {
        os << *(t.m_v[0]) << ", " << *(t.m_v[1]) << ", " << *(t.m_v[2]);
        return os;
    }
};

void centroidTri(const Triangle & tri, trimesh::point & ctd)
{
    ctd = tri.v(0) + tri.v(1) + tri.v(2);
    ctd /= 3.0;
}

float areaTri(const Triangle & tri)
{
    return 0.5 * len(cross(tri.s(0).st(), tri.s(1).st()));
}

// check whether p and q are in the same side of st
// side.e = st = t - s, side.v[0] -> s
bool sameSide(const trimesh::point & p, const trimesh::point & q, const TriSide & side)
{
	auto cp = (side.st() CROSS (p - side.s())); // cp = st CROSS sp
	auto cq = (side.st() CROSS (q - side.t())); // cq = st CROSS sq
	if ((cp DOT cq) >= 0) {
		return true;
	} else {
		return false;
	}
}

// tv = { a, b, c }, ts = { ab, ac, bc }
// require: q and tv are co-planar
bool ptInsideTri2D(const trimesh::point & q, const Triangle & tri)
{
    auto aq = q - tri.v(0);
    Matrix2 U_1;
    trimesh::vec2 aq2, ret_params;
    const int comb[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
    for (int i = 0; i < 3; i++) {
        auto i1 = comb[i][0];
        auto i2 = comb[i][1];
        Matrix2 U = {
            { tri.s(0).st()[i1], tri.s(1).st()[i1] },
            { tri.s(0).st()[i2], tri.s(1).st()[i2] }
        };
        aq2.set(aq[i1], aq[i2]);
        if (!approx_eqt(trimesh::vec2(0, 0), aq2, 1e-9f) && inverse(U, U_1)) {
        // if (!approx_eqt(trimesh::vec2(0, 0), aq2, 1e-12f) && inverse(U, U_1, true)) {
            matrixPrd(U_1, aq2, ret_params);
            // std::cout << ret_params[0] << " " << ret_params[1] << std::endl;
            // TODO: a bug here: why previously it works??? (should be 1e-6, but previously it was wrongly written as 1e6)?????? Amazing
            if (-1e-6 <= ret_params[0] && ret_params[0] <= 1 + 1e-6 && -1e-6 <= ret_params[1] && ret_params[1] <= 1 + 1e-6 && ret_params[0] + ret_params[1] <= 1 + 1e-6) {
                return true;
            } else {
                return false;
            }
        }
    }
    return false; // TODO: added to kill warning, not verified

    // // Below: a faster way (actually slower? not verified)
    // if (sameSide(q, (*tv)[0], (*tv)[1], (*tv)[2]) && sameSide(q, (*tv)[1], (*tv)[0], (*tv)[2]) && sameSide(q, (*tv)[2], (*tv)[0], (*tv)[1])) {
    //  return true;
    // } else {
    //  return false;
    // }
}

// Return the projection point from q to the plane of triangle tv
// let abc be the triangle tv, and let r be the return
// then ar = alpha * ab + beta * ac
float projDist2(const trimesh::point & q, const Triangle & tri, trimesh::point & proj)
{
    // fast method: the projection is only affected by the one vertex of the triangle when normal is supplied
    float gamma_p = (tri.n() DOT (tri.v(0) - q));
    float sq_norm_tn = len2(tri.n()); // TODO: ||n|| should always be 1
    float gamma = gamma_p / sq_norm_tn;
    proj = q + gamma * tri.n();
    return (gamma_p * gamma);
}


// tv = { a, b, c }, ts = { ab, ac, bc }
// Return whether q's projection to the plane of triangle tv is inside tv
// proj is set to the projection point
bool projInsideTri(const trimesh::point & q, const Triangle & tri, trimesh::point & proj, float & retProjDist2, float currentMinDist2 = INFTY)
{
    // std::cout << "Cal whether proj of " << q << " inside " << tri << std::endl;
    retProjDist2 = projDist2(q, tri, proj);
    // std::cout << "Calculated proj dist2: " << retProjDist2 << std::endl;
    if (retProjDist2 > currentMinDist2) {
        return false;
    } else {
        bool inside2D = ptInsideTri2D(proj, tri);
        // std::cout << "Result: " << inside2D << std::endl;
        return inside2D;
    }
}

// side.e = st = t - s, side.v[0] -> s
// Return the SQ dist between q to line segment st (t - s)
// let x be the nn
// then lambda is set such that sx = lambda * st
// lambda will be in [0, 1]
// lambda = 0 means x = s
// lambda = 1 means x = t
float p2lDist2(const trimesh::point & q, const TriSide & side, float & lambda)
{
    auto sq = q - side.s();
    lambda = (side.st() DOT sq) / side.l2();
    if (0 <= lambda && lambda <= 1) {
        return trimesh::len2(lambda * side.st() - sq);
    } else if (lambda < 0) { // nn is s
        lambda = 0;
        return trimesh::len2(sq);
    } else { // nn is t
        lambda = 1;
        return trimesh::len2(q - side.t());
    }
}

float p2lDist2(const trimesh::point & q, const trimesh::point & s, const trimesh::point & t, float & lambda)
{
    const trimesh::point p[2] = { s, t };
    TriSide side(p);
    return p2lDist2(q, side, lambda);
}

float p2lDist2_bf(const trimesh::point & q, const trimesh::point & s, const trimesh::point & t, float & lambda)
{
    float min_sq_dist = INFTY;
    const int num_steps = 10000;
    const float step = 1.0 / (float) num_steps;
    for (int i = 0; i <= num_steps; i++) {
        const float tmp_l = (float) i * step;
        const float tmp_sq_dist = dist2(q, s + tmp_l * (t - s));
        if (min_sq_dist > tmp_sq_dist) {
            min_sq_dist = tmp_sq_dist;
            lambda = tmp_l;
        }
    }
    return min_sq_dist;
}

// Return the SQ dist between trimesh::point q to triangle tv[3]
// nn is set to the NN point on the triangle to q
// float p2fDist2(const trimesh::point &q, const TrianglePtr &tv, const TriangleInSides &ts, const trimesh::vec &tn, trimesh::point &nn)
float p2fDist2(const trimesh::point & q, const Triangle & tri, trimesh::point & nn)
{
    // get projection of q to plane defined by tv[3]
    // if projection inside triangle, then the projection distance must be the shortest
    trimesh::point proj;
    float projDist2;
    if (projInsideTri(q, tri, proj, projDist2)) { // proj inside triangle
        // std::cout << "q " << q << " proj inside tri with proj dist2: " << projDist2 << std::endl;
        nn = proj;
        return projDist2;
    }
    // std::cout << "q proj not inside tri, continue searching for edges" << std::endl;

    // then, the shortest distance must appear on edge
    float minDist2 = INFTY;
    float tmpLambda, tmpDist2;
    for (int i = 0; i < 3; i++) {
        tmpDist2 = p2lDist2(q, tri.s(i), tmpLambda);
        if (minDist2 > tmpDist2) {
            minDist2 = tmpDist2;
            if (approx_eqt(tmpLambda, 0.0f)) {
                nn = tri.s(i).s();
            } else if (approx_eqt(tmpLambda, 1.0f)) {
                nn = tri.s(i).t();
            } else {
                nn = tri.s(i).s() + tmpLambda * tri.s(i).st();
            }
        }
    }

    return minDist2;
}

float p2fDist2(const trimesh::point & q, const trimesh::point (&tv)[3], trimesh::point & nn)
{
    Triangle tri(tv);
    return p2fDist2(q, tri, nn);
}

float p2fDist2_bf(const trimesh::point & q, const trimesh::point (&tv)[3], trimesh::point & nn)
{
    float min_sq_dist = INFTY;
    const int num_steps = 1000;
    const float step = 1.0 / (float) num_steps;
    for (int i = 0; i <= num_steps; i++) {
        for (int j = 0; j <= num_steps - i; j++) {
            const float tmp_al = (float) i * step;
            const float tmp_be = (float) j * step;
            const trimesh::point tmp_res = tv[0] + tmp_al * (tv[1] - tv[0]) + tmp_be * (tv[2] - tv[0]);
            const float tmp_sq_dist = dist2(q, tmp_res);
            if (min_sq_dist > tmp_sq_dist) {
                min_sq_dist = tmp_sq_dist;
                nn = tmp_res;
            }
        }
    }
    return min_sq_dist;
}

float maxP2fDist2(const trimesh::point & q, const Triangle & tri)
{
    return std::max({ dist2(q, tri.v(0)), dist2(q, tri.v(1)), dist2(q, tri.v(2)) });
}

// TODO: a more efficient way
float p2mDist2(const trimesh::point & q, const std::vector<Triangle> & triList, trimesh::point & nn, int & nnTriID)
{
    float minDist2 = INFTY;
    // std::cout << "p2fDist2:";
    for (int i = 0; i < triList.size(); i++) {
        trimesh::point currNN;
        float currDist2 = p2fDist2(q, triList[i], currNN);
        // std::cout << " " << currDist2;
        if (minDist2 > currDist2) {
            minDist2 = currDist2;
            nnTriID = i;
            nn = currNN;
        }
    }
    // std::cout << std::endl;
    return minDist2;
}

// side1: ab, side2: cd
// Return the SQ dist between line segment ab to cd
float e2eDist2(const TriSide & side1, const TriSide & side2, trimesh::point & nn1, trimesh::point & nn2)
{
    const trimesh::vec &ab = side1.st();
    const trimesh::vec &cd = side2.st();
    const trimesh::vec ca = side1.s() - side2.s();
    const float ab_dot_cd = (ab DOT cd);
    const float sq_norm_ab = side1.l2();
    const float sq_norm_cd = side2.l2();
    const float ab_dot_ca = (ab DOT ca);

    float alpha, beta;

    if (sq_norm_ab * sq_norm_cd != sq(ab_dot_cd)) { // ab not parallel to cd
        beta = (sq_norm_ab * (cd DOT ca) - (ab_dot_ca) * (ab_dot_cd)) / (sq_norm_ab * sq_norm_cd - sq(ab_dot_cd));
        alpha = (beta * ab_dot_cd - ab_dot_ca) / sq_norm_ab;

        if (0 <= alpha && alpha <= 1 && 0 <= beta && beta <= 1) {
            nn1 = side1.s() + ab * alpha;
            nn2 = side2.s() + cd * beta;
            return dist2(nn1, nn2);
        }
    }

    float lambda;
    float tmpDist2;
    float minDist2 = INFTY;

    tmpDist2 = p2lDist2(side1.s(), side2, lambda);
    if (tmpDist2 < minDist2) {
        minDist2 = tmpDist2;
        alpha = 0;
        beta = lambda;
    }

    tmpDist2 = p2lDist2(side1.t(), side2, lambda);
    if (tmpDist2 < minDist2) {
        minDist2 = tmpDist2;
        alpha = 1;
        beta = lambda;
    }
    
    tmpDist2 = p2lDist2(side2.s(), side1, lambda);
    if (tmpDist2 < minDist2) {
        minDist2 = tmpDist2;
        alpha = lambda;
        beta = 0;
    }
    
    tmpDist2 = p2lDist2(side2.t(), side1, lambda);
    if (tmpDist2 < minDist2) {
        minDist2 = tmpDist2;
        alpha = lambda;
        beta = 1;
    }

    nn1 = side1.s() + ab * alpha;
    nn2 = side2.s() + cd * beta;
    return minDist2;
}

float e2eDist2(const trimesh::point &a, const trimesh::point &b, const trimesh::point &c, const trimesh::point &d, trimesh::point &nn1, trimesh::point &nn2)
{
    trimesh::point p1[2] = { a, b };
    trimesh::point p2[2] = { c, d };
    TriSide side1(p1);
    TriSide side2(p2);
    return e2eDist2(side1, side2, nn1, nn2);
}

float e2eDist2_bf(const trimesh::point &a, const trimesh::point &b, const trimesh::point &c, const trimesh::point &d, trimesh::point &nn1, trimesh::point &nn2)
{
    float min_sq_dist = INFTY;
    float tmp_sq_dist;
    float tmp_al, tmp_be, alpha, beta;
    const int num_steps = 1000;
    const float step = 1.0 / float(num_steps);
    for (int i = 0; i <= num_steps; i++) {
        for (int j = 0; j <= num_steps; j++) {
            tmp_al = float(i * step);
            tmp_be = float(j * step);
            tmp_sq_dist = dist2(a + tmp_al * (b - a), c + tmp_be * (d - c));
            if (min_sq_dist > tmp_sq_dist) {
                min_sq_dist = tmp_sq_dist;
                alpha = tmp_al;
                beta = tmp_be;
            }
        }
    }

    nn1 = a + alpha * (b - a);
    nn2 = c + beta * (d - c);
    return min_sq_dist;
}

// side1: ab, side2: cd
bool lineIntersLine(const TriSide &side1, const TriSide &side2, trimesh::point &inters)
{
    const trimesh::vec &ab = side1.st();
    const trimesh::vec dc = -side2.st();
    const trimesh::vec ac = side2.s() - side1.s();

    if (sq((ab DOT dc)) == (side1.l2() * side2.l2())) { // ab // cd: ab DOT dc = ||ab|| ||dc||
        if (len2((ab CROSS ac)) > 1e-9) { // but not co-linear
            return false;
        } else { // co-linear
            for (int i = 0; i < 3; i++) {
                if (dc[i] == 0) {
                    continue;
                }
                const float al = ac[i] / dc[i];
                const float be = (side2.s()[i] - side1.t()[i]) / dc[i]; // c - b
                if (al < 0) {
                    if (be < 0) {
                        return false;
                    } else {
                        inters = side2.s();
                        return true;
                    }
                } else if (0 <= al && al <= 1) {
                    inters = side1.s();
                    return true;
                } else {
                    if (be <= 1) {
                        inters = side2.t();
                        return true;
                    } else {
                        return false;
                    }
                }
            }
        }
    }

    Matrix2 U_1;
    trimesh::vec2 ac2, ret_params;
    const int comb[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
    for (int i = 0; i < 3; i++) {
        auto s1 = comb[i][0];
        auto s2 = comb[i][1];

        Matrix2 U = { { ab[s1], dc[s1] }, { ab[s2], dc[s2] } };
        ac2.set(ac[s1], ac[s2]);
        
        if (inverse(U, U_1) && !approx_eqt(trimesh::vec2(0, 0), ac2, 1e-9f)) {
            matrixPrd(U_1, ac2, ret_params);
            if (0 <= ret_params[0] && ret_params[0] <= 1 && 0 <= ret_params[1] && ret_params[1] <= 1) {
                inters = side1.s() + ret_params[0] * ab;
                return true;
            } else {
                return false;
            }
        }
    }

    return false; // TODO: just to kill warning, not verified
}

// require: all the 5 points are coplanar
bool lineIntersTri2D(const TriSide &side, const Triangle &tri, trimesh::point &inters)
{
    if (ptInsideTri2D(side.s(), tri)) {
        inters = side.s();
        return true;
    }
    if (ptInsideTri2D(side.t(), tri)) {
        inters = side.t();
        return true;
    }

    for (int i = 0; i < 3; i++) {
        if (lineIntersLine(tri.s(i), side, inters)) {
            return true;
        }
    }

    return false;
}

// return whether side st intersects triangle t
bool lineIntersTri(const TriSide &side, const Triangle &tri, trimesh::point &inters)
{
    auto sgn1 = approx_sgn((tri.n() DOT (side.s() - tri.v(0))), 1e-9);
    auto sgn2 = approx_sgn((tri.n() DOT (side.t() - tri.v(0))), 1e-9);
	if (sgn1 == sgn2) {
		if (sgn1 != 0) { // the case when s and t in different sides of face (coplanar excluded)
			return false;
		} else { // the case when the 5 points are coplanar, then we perform the 2d check
			return lineIntersTri2D(side, tri, inters);
		}
	}

    // which is defined to be inters - lv1 = gamma * (lv2 - lv1)
    float gamma = (tri.n() DOT (tri.v(0) - side.s())) / (tri.n() DOT side.st());

    inters = side.s() + side.st() * gamma;
    return ptInsideTri2D(inters, tri);
}

// Return the shortest dist2 between triangle t1 and t2
// nn1 and nn2 will be set to the nn points on t1 and t2 respectively
// float nn_triangle(const trimesh::point (&tv1)[3], const trimesh::point (&tv2)[3], const trimesh::vec &tn1, const trimesh::vec &tn2, trimesh::point &nn1, trimesh::point &nn2, float time_consumption[3])
// float nn_triangle(const trimesh::point (&tv1)[3], const trimesh::point (&tv2)[3], const trimesh::vec &tn1, const trimesh::vec &tn2, trimesh::point &nn1, trimesh::point &nn2)
float nnTriDist2(const Triangle &tri1, const Triangle &tri2, trimesh::point &nn1, trimesh::point &nn2)
{
    // timer_start();

    // First check if tv1 and tv2 intersects <==> at least one adge of one triangle intersects another triangle\
    // Very first, check if any point-pair coincident
    for (int i = 0; i < 3; i++) {
    	for (int j = 0; j < 3; j++) {
    		if (approx_eqt(tri1.v(i), tri2.v(j), 1e-9f)) {
    			nn1 = nn2 = tri1.v(i);
    			return 0.0;
    		}
    	}
    }

    bool bb_inters = true;
    // // First of first, if the bounding box not intersecting, bypass
    // const float tv1_bb_min[3] = {
    //     min({tri1.v(0)[0], tri1.v(1)[0], tri1.v(2)[0]}),
    //     min({tri1.v(0)[1], tri1.v(1)[1], tri1.v(2)[1]}),
    //     min({tri1.v(0)[2], tri1.v(1)[2], tri1.v(2)[2]})
    // };
    // const float tv1_bb_max[3] = {
    //     max({tri1.v(0)[0], tri1.v(1)[0], tri1.v(2)[0]}),
    //     max({tri1.v(0)[1], tri1.v(1)[1], tri1.v(2)[1]}),
    //     max({tri1.v(0)[2], tri1.v(1)[2], tri1.v(2)[2]})
    // };
    // const float tv2_bb_min[3] = {
    //     min({tri2.v(0)[0], tri2.v(1)[0], tri2.v(2)[0]}),
    //     min({tri2.v(0)[1], tri2.v(1)[1], tri2.v(2)[1]}),
    //     min({tri2.v(0)[2], tri2.v(1)[2], tri2.v(2)[2]})
    // };
    // const float tv2_bb_max[3] = {
    //     max({tri2.v(0)[0], tri2.v(1)[0], tri2.v(2)[0]}),
    //     max({tri2.v(0)[1], tri2.v(1)[1], tri2.v(2)[1]}),
    //     max({tri2.v(0)[2], tri2.v(1)[2], tri2.v(2)[2]})
    // };
    // for (int i = 0; i < 3; i++) {
    // 	if (tv1_bb_max[i] < tv2_bb_min[i] || tv1_bb_min[i] > tv2_bb_max[i]) {
    // 		bb_inters = false; // if any dimension non-intersecting, then it must be non-intersecting
    // 		break;
    // 	}
    // }

    // if (bb_inters) { // if only the bounding box not intersecting, then we perform this intersection check
	    trimesh::point inters;
	    for (int i = 0; i < 3; i++) {
	        if (lineIntersTri(tri1.s(i), tri2, inters) || lineIntersTri(tri2.s(i), tri1, inters)) {
	            nn1 = nn2 = inters;
				// time_consumption[0] += timer_end(SECOND);
	            return 0.0;
	        }
	    }
	// }

	// time_consumption[0] += timer_end(SECOND);

    float minDist2 = INFTY;

    // timer_start();

    // then check each vertex projection to another face, then check each e2e
    trimesh::point proj;
    float projDist2;
    for (int i = 0; i < 3; i++) {
        if (projInsideTri(tri1.v(i), tri2, proj, projDist2, minDist2)) { // return means projDist2 must be smaller than minDist2
            minDist2 = projDist2;
            nn1 = tri1.v(i);
            nn2 = proj;
            // std::cout << "When proj " << tri1.v(i) << " to [" << tri2 << "]:" << std::endl;
            // std::cout << "proj hit: " << minDist2 << std::endl;
        }
        if (projInsideTri(tri2.v(i), tri1, proj, projDist2, minDist2)) {
            minDist2 = projDist2;
            nn1 = proj;
            nn2 = tri2.v(i);
            // std::cout << "When proj " << tri2.v(i) << " to [" << tri1 << "]:" << std::endl;
            // std::cout << "proj hit: " << minDist2 << std::endl;
        }
    }

	// time_consumption[1] += timer_end(SECOND);

    // timer_start();

    // then check each e2e
    float tmpE2eDist2;
    trimesh::point tmpE2eNN1, tmpE2eNN2;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            tmpE2eDist2 = e2eDist2(tri1.s(i), tri2.s(j), tmpE2eNN1, tmpE2eNN2);
            if (minDist2 > tmpE2eDist2) {
                minDist2 = tmpE2eDist2;
                nn1 = tmpE2eNN1;
                nn2 = tmpE2eNN2;
                // std::cout << "e2e hit: " << minDist2 << std::endl;
            }
        }
    }

	// time_consumption[2] += timer_end(SECOND);

    return minDist2;
}

float nnTriDist2(const Triangle &tri1, const Triangle &tri2)
{
    trimesh::point nn1, nn2;
    return nnTriDist2(tri1, tri2, nn1, nn2);
}

float nnTriDist2(const trimesh::point (&tv1)[3], const trimesh::point (&tv2)[3], trimesh::point &nn1, trimesh::point &nn2)
{
    Triangle tri1(tv1);
    Triangle tri2(tv2);
    return nnTriDist2(tri1, tri2, nn1, nn2);
}

float nnTriDist2_bf(const trimesh::point (&tv1)[3], const trimesh::point (&tv2)[3], trimesh::point &nn1, trimesh::point &nn2)
{
    float min_sq_dist = INFTY;
    // int min_i, min_j, min_k, min_l;
    const int num_steps = 100;
    const float step = 1.0 / (float) num_steps;
    for (int i = 0; i <= num_steps; i++) {
        for (int j = 0; j <= num_steps - i; j++) {
            for (int k = 0; k <= num_steps; k++) {
                for (int l = 0; l <= num_steps - k; l++) {
                    const float tmp_al = (float) i * step;
                    const float tmp_be = (float) j * step;
                    const float tmp_la = (float) k * step;
                    const float tmp_mu = (float) l * step;
                    const trimesh::point tmp_nn1 = tv1[0] + tmp_al * (tv1[1] - tv1[0]) + tmp_be * (tv1[2] - tv1[0]);
                    const trimesh::point tmp_nn2 = tv2[0] + tmp_la * (tv2[1] - tv2[0]) + tmp_mu * (tv2[2] - tv2[0]);
                    const float tmp_sq_dist = dist2(tmp_nn1, tmp_nn2);
                    if (min_sq_dist > tmp_sq_dist) {
                        min_sq_dist = tmp_sq_dist;
                        nn1 = tmp_nn1;
                        nn2 = tmp_nn2;
                        // min_i = i;
                        // min_j = j;
                        // min_k = k;
                        // min_l = l;
                    }
                }
            }
        }
    }

    // std::cout << min_i << " " << min_j << " " << min_k << " " << min_l << std::endl;

    return min_sq_dist;
}

float maxTriDist2(const Triangle &tri1, const Triangle &tri2, trimesh::point &nn1, trimesh::point &nn2)
{
    float maxDist2 = -1;
    int iMax = -1, jMax = -1;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            float currDist2 = dist2(tri1.v(i), tri2.v(j));
            if (maxDist2 < currDist2) {
                maxDist2 = currDist2;
                iMax = i;
                jMax = j;
            }
        }
    }

    nn1 = tri1.v(iMax);
    nn2 = tri1.v(jMax);
    return maxDist2;
}

float maxTriDist2(const Triangle &tri1, const Triangle &tri2)
{
    trimesh::point nn1, nn2;
    return maxTriDist2(tri1, tri2, nn1, nn2);
}


// Return the Hausdorff dist from triangle tv1 to tv2
// ref1 and ref2 will be set to the ref points on tv1 and tv2 respectively
float diHausDist2(const Triangle &tri1, const Triangle &tri2, trimesh::point &nn1, trimesh::point &nn2)
{
    float maxDist2 = -1;
    int iMax = -1;
    for (int i = 0; i < 3; i++) {
        trimesh::point currNN2;
        float currDist2 = p2fDist2(tri1.v(i), tri2, currNN2);
        if (maxDist2 < currDist2) {
            maxDist2 = currDist2;
            iMax = i;
            nn2 = currNN2;
        }
    }

    nn1 = tri1.v(iMax);
    return maxDist2;
}

float diHausDist2_bf(const Triangle &tri1, const std::vector<Triangle>& tri2List, trimesh::point &nn1, trimesh::point &nn2, int &nn2TriID)
{
    float maxP2mDist2 = -1;
    const int numSteps = 100;
    const float step = 1.0 / (float) numSteps;
    for (int i = 0; i <= numSteps; i++) {
        for (int j = 0; j <= numSteps - i; j++) {
            auto tmpQ1 = tri1.v(0) + (float(i) * step) * (tri1.v(1) - tri1.v(0)) + (float(j) * step) * (tri1.v(2) - tri1.v(0));
            trimesh::point tmpNN2;
            int tmpNN2TriID;
            float tmpP2mDist2 = p2mDist2(tmpQ1, tri2List, tmpNN2, tmpNN2TriID);
            // std::cout << "tmpP2mDist2: " << tmpP2mDist2 << std::endl;
            // if (tmpP2mDist2 > 1.0) {
            //     std::cout << "i = " << i << ", j = " << j << ", Calc tmpP2mDist2 from " << tmpQ1 << ": " << tmpP2mDist2 << std::endl;
            // }
            if (maxP2mDist2 < tmpP2mDist2) {
                maxP2mDist2 = tmpP2mDist2;
                nn1 = tmpQ1;
                nn2 = tmpNN2;
                nn2TriID = tmpNN2TriID;
            }
        }
    }
    return maxP2mDist2;
}

float diHausDist2_bf(const std::vector<Triangle>& tri1List, const std::vector<Triangle>& tri2List, trimesh::point &nn1, int &nn1TriID, trimesh::point &nn2, int &nn2TriID)
{
    float maxF2mDist2 = -1;
    for (int i = 0; i < tri1List.size(); i++) {
        trimesh::point currNN1, currNN2;
        int currNN2TriID;
        float currF2mDist2 = diHausDist2_bf(tri1List[i], tri2List, currNN1, currNN2, currNN2TriID);
        // std::cout << "currF2mDist2: " << currF2mDist2 << std::endl;
        if (maxF2mDist2 < currF2mDist2) {
            maxF2mDist2 = currF2mDist2;
            nn1 = currNN1;
            nn2 = currNN2;
            nn1TriID = i;
            nn2TriID = currNN2TriID;
        }
    }
    return maxF2mDist2;
}


float itgTriDist_bf(const Triangle &tri1, const Triangle &tri2)
{
    float sum_dist = 0.0;
    int count = 0;
    const int num_steps = 100;
    const float step = 1.0 / float(num_steps);
    trimesh::point tmp_p;
    for (int i = 0; i <= num_steps; i++) {
        for (int j = 0; j <= num_steps - i; j++) {
            const float tmp_al = float(i) * step;
            const float tmp_be = float(j) * step;
            const trimesh::point tmp_q = tri1.v(0) + tmp_al * (tri1.v(1) - tri1.v(0)) + tmp_be * (tri1.v(2) - tri1.v(0));
            sum_dist += sqrt(p2fDist2(tmp_q, tri2, tmp_p));
            count++;
        }
    }

    return (sum_dist / float(count));
}

float itgTriDist(const Triangle &tri1, const Triangle &tri2)
{
    // TODO
    return itgTriDist_bf(tri1, tri2);
}


#endif