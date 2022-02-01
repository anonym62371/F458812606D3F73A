#ifndef __UTILTRANS_H_
#define __UTILTRANS_H_

#include <iostream>

#include <Eigen/Dense>

#include "Vec.h"
#include "XForm.h"

#include "UtilPolygon.hpp"


void print_trans(const trimesh::xform & t, std::ostream & os)
{
    for (int i = 0; i < 4; i++) {
        os << t[i*4] << " " << t[i*4+1] << " " << t[i*4+2] << " " << t[i*4+3] << std::endl;
    } // a normal way
}

void trans_pt(const trimesh::xform & xf, const trimesh::point & p, trimesh::point & xf_p)
{
    xf_p = xf * p;
}

void trans_tri(const trimesh::xform & xf, const Triangle & tri, Triangle & xf_tri)
{
    const trimesh::point* tri_v[3];
    for (int i = 0; i < 3; i++) {
        tri_v[i] = new trimesh::point(xf * tri.v(i));
    }
    xf_tri = Triangle(tri_v);
}


void cal_trans(const trimesh::point** q, const trimesh::point** p, int num, trimesh::xform& ret)
{
    Eigen::Vector3f vq[num], vp[num];
    Eigen::Vector3f sum_q(0, 0, 0), sum_p(0, 0, 0);
    for (int i = 0; i < num; i++) {
        vq[i] << (*q[i])[0], (*q[i])[1], (*q[i])[2];
        vp[i] << (*p[i])[0], (*p[i])[1], (*p[i])[2];
        sum_q += vq[i];
        sum_p += vp[i];
    }
    Eigen::Vector3f vq_bar = sum_q / num;
    Eigen::Vector3f vp_bar = sum_p / num;

    // x = vq - vq_bar, y = vp - vp_bar
    Eigen::MatrixXf x(3, num), y(3, num);
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < 3; j++) {
            x(j, i) = vq[i](j) - vq_bar(j);
            y(j, i) = vp[i](j) - vp_bar(j);
        }
    }

    // S = X * Y^T, SVD(S) = U * D * V^T
    Eigen::Matrix3f s = x * y.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3f> svd(s, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // R = V * U^T, t = vp_bar - R * vq_bar
    Eigen::Matrix3f r = svd.matrixV() * svd.matrixU().transpose();
    Eigen::Vector3f t = vp_bar - r * vq_bar;

    // ret = {
    //     r(0, 0), r(0, 1), r(0, 2),
    //     r(1, 0), r(1, 1), r(1, 2),
    //     r(2, 0), r(2, 1), r(2, 2),
    //     t(0),    t(1),    t(2)
    // };

    ret = trimesh::xform(
        r(0, 0), r(1, 0), r(2, 0), 0,
        r(0, 1), r(1, 1), r(2, 1), 0,
        r(0, 2), r(1, 2), r(2, 2), 0,
        t(0),    t(1),    t(2),    1 
    );
}



inline const trimesh::xform angle_axis_to_xform(float rx, float ry, float rz, float tx, float ty, float tz)
{
    float t = sqrt(rx*rx + ry*ry + rz*rz);
    if (t > 0) {
        rx /= t;
        ry /= t;
        rz /= t;

        float ct = cos(t);
        float ct2 = 1 - ct;
        float st = sin(t);
        float st2 = 1 - st;

        float tmp121 = rx*ry*ct2; float tmp122 = rz*st;
        float tmp131 = rx*rz*ct2; float tmp132 = ry*st;
        float tmp231 = ry*rz*ct2; float tmp232 = rx*st;

        float R11 = ct + rx*rx*ct2;  float R12 = tmp121 - tmp122; float R13 = tmp131 + tmp132;
        float R21 = tmp121 + tmp122; float R22 = ct + ry*ry*ct2;  float R23 = tmp231 - tmp232;
        float R31 = tmp131 - tmp132; float R32 = tmp231 + tmp232; float R33 = ct + rz*rz*ct2;

        return trimesh::xform(
            R11, R12, R13, tx,
            R21, R22, R23, ty,
            R31, R32, R33, tz,
            0.0, 0.0, 0.0, 1.0
        );
    }
    // If t == 0, the rotation angle is 0 and no rotation is required
    else {
        return trimesh::xform(
            1.0, 0.0, 0.0, tx,
            0.0, 1.0, 0.0, ty,
            0.0, 0.0, 1.0, tz,
            0.0, 0.0, 0.0, 1.0
        ); // TODO: incorrect?
    }
}

inline const trimesh::xform angle_axis_to_xform(const trimesh::point& r, const trimesh::point& t)
{
    return angle_axis_to_xform(r[0], r[1], r[2], t[0], t[1], t[2]);
}

inline const trimesh::xform to_xform_matrix(float yaw, float pitch, float roll, float x, float y, float z)
{
    float cy = cos(yaw);
    float sy = sin(yaw);
    float cp = cos(pitch);
    float sp = sin(pitch);
    float cr = cos(roll);
    float sr = sin(roll);

    return trimesh::xform(
        cy*cp, cy*sp*sr-sy*cr, cy*sp*cr+sy*sr,   x,
        sy*cp, sy*sp*sr+cy*cr, sy*sp*cr-cy*sr,   y,
          -sp,          cp*sr,          cp*cr,   z,
          0.0,            0.0,            0.0, 1.0
    ); // TODO: incorrect?
}

void show_yprxyz(trimesh::xform& xf)
{

}


struct Quaternion
{
    float w, x, y, z;
};

Quaternion to_quaternion(float yaw, float pitch, float roll) // yaw (Z), pitch (Y), roll (X)
{
    // Abbreviations for the various angular functions
    float cy = cos(yaw * 0.5);
    float sy = sin(yaw * 0.5);
    float cp = cos(pitch * 0.5);
    float sp = sin(pitch * 0.5);
    float cr = cos(roll * 0.5);
    float sr = sin(roll * 0.5);

    Quaternion q;
    q.w = cr * cp * cy + sr * sp * sy;
    q.x = sr * cp * cy - cr * sp * sy;
    q.y = cr * sp * cy + sr * cp * sy;
    q.z = cr * cp * sy - sr * sp * cy;

    return q;
}

#endif