#ifndef __UTILVEC_H_
#define __UTILVEC_H_

#include <cmath>
#include <random>
#include <cstdlib>
#include <vector>
#include <iterator>
#include <algorithm>

#include "Vec.h"
#include "TriMesh.h"

#define SQ_PI 9.869604401f
#define SQRT_3 1.73205080757f

#define APPROX_ZERO 1e-6f
#define INFTY 3.3e33f

// Their order matters
uint64_t getNumHash(uint32_t a, uint32_t b)
{
	uint64_t hash;
	uint32_t sum1 = a + b;
	hash = sum1 * (sum1 + 1) / 2 + b;
	return hash;
}

// Their order matters
uint64_t getNumHash(uint32_t a, uint32_t b, uint32_t c)
{
	uint64_t hash;
	uint32_t sum1 = a + b;
	hash = sum1 * (sum1 + 1) / 2 + b;
	uint64_t sum2 = hash + c;
	hash = sum2 * (sum2 + 1) / 2 + c;
	return hash;
}

inline int approx_sgn(float a, float zero = APPROX_ZERO)
{
	if (a > zero) {
		return 1;
	} else if (a < -zero) {
		return -1;
	} else {
		return 0;
	}
}

inline bool approx_eqt(float a, float b, float zero = APPROX_ZERO)
{
    return (abs(a - b) <= zero);
}

inline bool approx_eqt(const trimesh::vec2 & a, const trimesh::vec2 & b, float zero = APPROX_ZERO)
{
    return (approx_eqt(a[0], b[0], zero) && approx_eqt(a[1], b[1], zero));
}

inline bool approx_eqt(const trimesh::vec & a, const trimesh::vec & b, float zero = APPROX_ZERO)
{
    return (approx_eqt(a[0], b[0], zero) && approx_eqt(a[1], b[1], zero) && approx_eqt(a[2], b[2], zero));
}

inline float ran_f_0_1()
{
	return ((float) rand() / RAND_MAX);
}

// return a random float in (-scaling, scaling)
inline float ran_f(float scaling = 1.0)
{
    return (2.0 * scaling * (float) rand() / RAND_MAX - scaling);
}

// return a random trimesh::point in (-1, 1)
inline trimesh::point ran_p(float scaling = 1.0)
{
    return trimesh::point(ran_f(scaling), ran_f(scaling), ran_f(scaling));
}

inline void ran_p(trimesh::point &ret, float scaling = 1.0)
{
	ret[0] = ran_f(scaling);
	ret[1] = ran_f(scaling);
	ret[2] = ran_f(scaling);
}

inline void n_rand_k(int n, int k, std::vector<int> & ret)
{
	std::vector<int> nList;
	for (int i = 0; i < n; i++) {
		nList.push_back(i);
	}

	std::sample(nList.begin(), nList.end(), std::back_inserter(ret), k, std::mt19937{std::random_device{}()});
	sort(ret.begin(), ret.end());
}

inline float sq(float v)
{
	return v*v;
}

inline trimesh::point middle_pt(const trimesh::point & a, const trimesh::point & b)
{
    return trimesh::point(0.5*(a[0]+b[0]), 0.5*(a[1]+b[1]), 0.5*(a[2]+b[2]));
}

inline void middle_pt(const trimesh::point & a, const trimesh::point & b, trimesh::point & ret)
{
	ret[0] = 0.5*(a[0]+b[0]);
	ret[1] = 0.5*(a[1]+b[1]);
	ret[2] = 0.5*(a[2]+b[2]);
}

inline void normalize(trimesh::point & p, float to = 1.0)
{
    float length = len(p);
    float ratio = to / length;
    p *= ratio;
}

inline trimesh::point normalize(const trimesh::point & p, float to = 1.0)
{
	trimesh::point ret;
	normalize(ret, to);
    return ret;
}

inline float box_length(const trimesh::TriMesh::BBox & box)
{
	return std::max({box.max[0] - box.min[0], box.max[1] - box.min[1], box.max[2] - box.min[2]});
}

#endif