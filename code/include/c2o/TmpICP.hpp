#ifndef __TMPICP_H_
#define __TMPICP_H_

#include "TriMesh.h"
#include "KDtree.h"
#include "XForm.h"

#include "UtilMatrix.hpp"

#include <Eigen/Dense>

#include <iostream>

#define NUM_ITR 5

struct Corr_Element
{
	trimesh::point v;
	trimesh::point n;
	Corr_Element(trimesh::point _v, trimesh::point _n) {
		v = _v;
		n = _n;
	}
};

bool align_v2v(const std::vector<Corr_Element>& elements_data, const std::vector<Corr_Element>& elements_model, trimesh::xform& xf_incremental)
{
	int n = elements_data.size();

	// calculate centroids
	trimesh::point centroid_data, centroid_model;
	for (int i = 0; i < n; i++) {
		centroid_data += elements_data[i].v;
		centroid_model += elements_model[i].v;
	}
	centroid_data /= (float) n;
	centroid_model /= (float) n;

	float U[3][3] = { { 0 } };

	for (int i = 0; i < n; i++) {
		trimesh::dvec3 v_model = trimesh::dvec3(elements_model[i].v - centroid_model);
		trimesh::dvec3 v_data = trimesh::dvec3(elements_data[i].v - centroid_data);
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				U[j][k] += v_model[j] * v_data[k];
	}
	// print_matrix(U);

	float s[3], V[3][3];
	trimesh::svd<float,3,3>(U, s, V);
	// cout << "s: " << s[0] << " " << s[1] << " " << s[2] << endl;
	for (int i = 0; i < 3; i++) {
		if (s[i] < 1e-6)
			return false;
	}

	if ((determ(U) < 0.0) ^ (determ(V) < 0.0)) {
		V[2][0] = -V[2][0];
		V[2][1] = -V[2][1];
		V[2][2] = -V[2][2];
	}
	xf_incremental = trimesh::xform::trans(centroid_model) *
	          		 trimesh::xform::fromarray(U) * trimesh::transp(trimesh::xform::fromarray(V)) *
	          		 trimesh::xform::trans(-centroid_data);

	return true;
}

// Symmetric point-to-plane alignment.
bool align_symm(const std::vector<Corr_Element>& elements_data, const std::vector<Corr_Element>& elements_model, trimesh::xform& xf_incremental)
{
	int n = elements_data.size();

	// calculate centroids
	trimesh::point centroid_data, centroid_model;
	for (int i = 0; i < n; i++) {
		centroid_data += elements_data[i].v;
		centroid_model += elements_model[i].v;
	}
	centroid_data /= (float) n;
	centroid_model /= (float) n;

	float scale = 0.0f;
	for (int i = 0; i < n; i++) {
		scale += dist2(elements_data[i].v, centroid_data);
		scale += dist2(elements_model[i].v, centroid_model);
	}
	scale = sqrt(scale / (2 * n));
	scale = 1.0f / scale;

	double A[6][6] = { { 0 } }, b[6] = { 0 };

	for (int i = 0; i < n; i++) {
		trimesh::dvec3 p_data = trimesh::dvec3(scale * (elements_data[i].v - centroid_data));
		trimesh::dvec3 p_model = trimesh::dvec3(scale * (elements_model[i].v - centroid_model));
		trimesh::dvec3 n = trimesh::dvec3(elements_data[i].n + elements_model[i].n);
		trimesh::dvec3 p = p_data + p_model;
		trimesh::dvec3 c = p CROSS n;
		trimesh::dvec3 d = p_model - p_data;

		double x[6] = { c[0], c[1], c[2], n[0], n[1], n[2] };
		double dn = d DOT n;

		for (int j = 0; j < 6; j++) {
			b[j] += dn * x[j];
			for (int k = j; k < 6; k++)
				A[j][k] += x[j] * x[k];
		}
	}

	// Make matrix symmetric
	for (int j = 1; j < 6; j++)
		for (int k = 0; k < j; k++)
			A[j][k] = A[k][j];

	// Eigen-decomposition and inverse
	double eval[6], einv[6];
	trimesh::eigdc<double,6>(A, eval);
	for (int i = 0; i < 6; i++) {
		if (eval[i] < 1e-6)
			return false;
		else
			einv[i] = 1.0 / eval[i];
	}

	// Solve system
	trimesh::eigmult<double,6>(A, einv, b);

	// Extract rotation and translation
	trimesh::dvec3 rot(b[0], b[1], b[2]), trans(b[3], b[4], b[5]);
	double rotangle = atan(len(rot));
	trans *= 1.0 / scale;
	trans *= cos(rotangle);

	trimesh::xform R = trimesh::xform::rot(rotangle, rot);
	xf_incremental = trimesh::xform::trans(centroid_model) *
	        	 R * trimesh::xform::trans(trans) * R *
	        		 trimesh::xform::trans(-centroid_data);

	return true;
}

void tmp_icp_v2v(const trimesh::TriMesh* qMesh, trimesh::TriMesh* dbMesh, trimesh::xform& xf, const trimesh::KDtree* dbKd, int num_itr = NUM_ITR, bool verbose = false)
{
	orthogonalize(xf);

	for (int i = 0; i < num_itr; i++) {
		// if (verbose) {
		// 	std::cout << "Itr(" << i << "): ";
		// }

		std::vector<Corr_Element> elements_data;
		std::vector<Corr_Element> elements_model;

		// TODO: random select subset
		for (int j = 0; j < qMesh->vertices.size(); j++) {
			trimesh::point qV = xf * qMesh->vertices[j];
			auto nn = dbKd->closest_to_pt((float *) &qV);
			if (!nn) {
				continue;
			}

			int nnIdx = (const trimesh::point *) nn - &(dbMesh->vertices[0]);
			// std::cout << "nnIdx: " << nnIdx << std::endl;
			// std::cout << "dbMesh->vertices.size():" << dbMesh->vertices.size() << std::endl;
			// std::cout << "dbMesh->normals.size():" << dbMesh->normals.size() << std::endl;

			elements_data.emplace_back(qV, trimesh::point());
			elements_model.emplace_back(dbMesh->vertices[nnIdx], trimesh::point());
		}

		trimesh::xform xf_incremental;
		auto aligned = align_v2v(elements_data, elements_model, xf_incremental);
		if (aligned) {
			xf = xf_incremental * xf;
		}
	}
}

void tmp_icp_symm(const trimesh::TriMesh* qMesh, trimesh::TriMesh* dbMesh, trimesh::xform& xf, const trimesh::KDtree* dbKd, int num_itr = NUM_ITR, bool verbose = false)
{
	orthogonalize(xf);

	for (int i = 0; i < num_itr; i++) {
		// if (verbose) {
		// 	std::cout << "Itr(" << i << "): ";
		// }

		std::vector<Corr_Element> elements_data;
		std::vector<Corr_Element> elements_model;

		auto nxf = norm_xf(xf);
		// TODO: random select subset
		for (int j = 0; j < qMesh->vertices.size(); j++) {
			trimesh::point qV = xf * qMesh->vertices[j];
			auto nn = dbKd->closest_to_pt((float *) &qV);
			if (!nn) {
				continue;
			}

			int nnIdx = (const trimesh::point *) nn - &(dbMesh->vertices[0]);
			// std::cout << "nnIdx: " << nnIdx << std::endl;
			// std::cout << "dbMesh->vertices.size():" << dbMesh->vertices.size() << std::endl;
			// std::cout << "dbMesh->normals.size():" << dbMesh->normals.size() << std::endl;

			elements_data.emplace_back(qV, nxf * qMesh->normals[j]);
			elements_model.emplace_back(dbMesh->vertices[nnIdx], dbMesh->normals[nnIdx]);
		}

		trimesh::xform xf_incremental;
		auto aligned = align_symm(elements_data, elements_model, xf_incremental);
		if (aligned) {
			xf = xf_incremental * xf;
		}
	}
}

#endif