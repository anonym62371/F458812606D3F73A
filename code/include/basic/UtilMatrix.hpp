#ifndef __UTILMATRIX_H_
#define __UTILMATRIX_H_

#include <iostream>

#include <Eigen/Dense>

#include "Vec.h"

using Matrix4 = float[4][4];
using Matrix3 = float[3][3];
using Matrix2 = float[2][2];

inline void printMatrix(const Matrix3 &A)
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << " " << A[i][j];
        }
        std::cout << std::endl;
    }
}

inline void printMatrix(const Matrix2 &A)
{
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            std::cout << " " << A[i][j];
        }
        std::cout << std::endl;
    }
}

inline float determ(const Matrix2 &A)
{
    return (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
}

inline float determ(const Matrix3 &A)
{
	return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) +
	       A[0][1] * (A[1][2] * A[2][0] - A[1][0] * A[2][2]) +
	       A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}

inline float determ(const Matrix4 &A)
{
    Eigen::Matrix4d M;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            M(i, j) = A[i][j];
        }
    }
    return M.determinant();
}

// Function to get adjoint of A[N][N] in adj[N][N].
inline void adjoint(const Matrix3 &A, Matrix3 &adj)
{
	adj[0][0] = A[1][1] * A[2][2] - A[1][2] * A[2][1];
	adj[0][1] = A[0][2] * A[2][1] - A[0][1] * A[2][2];
	adj[0][2] = A[0][1] * A[1][2] - A[0][2] * A[1][1];

	adj[1][0] = A[2][0] * A[1][2] - A[1][0] * A[2][2];
	adj[1][1] = A[0][0] * A[2][2] - A[0][2] * A[2][0];
	adj[1][2] = A[0][2] * A[1][0] - A[0][0] * A[1][2];

	adj[2][0] = A[1][0] * A[2][1] - A[1][1] * A[2][0];
	adj[2][1] = A[0][1] * A[2][0] - A[0][0] * A[2][1];
	adj[2][2] = A[0][0] * A[1][1] - A[0][1] * A[1][0];
}

inline bool inverse(const Matrix2 &A, Matrix2 &inverse, bool verbose = false)
{
    float det_A = determ(A);
    if (det_A == 0.0) {
        if (verbose) {
            std::cout << "Singular 2x2 matrix, can't find its inverse" << std::endl;
        }
        return false;
    }

    inverse[0][0] =   A[1][1] / det_A;
    inverse[0][1] = - A[0][1] / det_A;
    inverse[1][0] = - A[1][0] / det_A;
    inverse[1][1] =   A[0][0] / det_A;

    return true;
}

inline bool inverse(const Matrix3 &A, Matrix3 &inverse, bool verbose = false)
{
	float det_A = determ(A);
	if (det_A == 0.0) {
        if (verbose) {
            std::cout << "Singular 3x3 matrix, can't find its inverse" << std::endl;
        }
        return false;
	}

    // Find adjoint
    Matrix3 adj;
    adjoint(A, adj);

    // Get inverse
    for (int i = 0; i < 3; i++) {
    	for (int j = 0; j < 3; j++) {
    		inverse[i][j] = adj[i][j] / det_A;
    	}
    }

    return true;
}


inline void makeMatrix(const trimesh::vec &a, const trimesh::vec &b, const trimesh::vec &c, Matrix3 &ret)
{
    for (int i = 0; i < 3; i++) {
        ret[i][0] = a[i];
        ret[i][1] = b[i];
        ret[i][2] = c[i];
    }
}

inline void matrixPrd(const Matrix2 &A, const trimesh::vec2 &p, trimesh::vec2 &ret)
{
    ret.set(
        A[0][0] * p[0] + A[0][1] * p[1],
        A[1][0] * p[0] + A[1][1] * p[1]
    );
}

inline trimesh::vec2 matrixPrd(const Matrix2 &A, const trimesh::vec2 &p)
{
    trimesh::vec2 ret;
    matrixPrd(A, p, ret);
    return ret;
}

inline void matrixPrd(const Matrix3 &A, const trimesh::vec &p, trimesh::vec &ret)
{
    ret.set(
        A[0][0] * p[0] + A[0][1] * p[1] + A[0][2] * p[2],
        A[1][0] * p[0] + A[1][1] * p[1] + A[1][2] * p[2],
        A[2][0] * p[0] + A[2][1] * p[1] + A[2][2] * p[2]
    );
}

inline trimesh::vec matrixPrd(const Matrix3 &A, const trimesh::vec &p)
{
    trimesh::vec ret;
    matrixPrd(A, p, ret);
    return ret;
}

#endif