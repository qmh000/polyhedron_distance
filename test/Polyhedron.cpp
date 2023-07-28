#include "Polyhedron.hpp"
#include "Triangle.hpp"
#include <float.h>
#include <algorithm>

//--------------------------------------------------------------------------
// SegPoints() 
//
// Returns closest points between an segment pair.
// Implemented from an algorithm described in
//
// Vladimir J. Lumelsky,
// On fast computation of distance between line segments.
// In Information Processing Letters, no. 21, pages 55-61, 1985.   
//--------------------------------------------------------------------------

void SegPoints(Eigen::Vector3f& VEC,
    Eigen::Vector3f& X, Eigen::Vector3f& Y,             // closest points
    const Eigen::Vector3f P, const Eigen::Vector3f A, // seg 1 origin, vector
    const Eigen::Vector3f Q, const Eigen::Vector3f B) // seg 2 origin, vector
{
    float A_dot_A, B_dot_B, A_dot_B, A_dot_T, B_dot_T;
    Eigen::Vector3f T;
    Eigen::Vector3f TMP;

    T = Q - P;
    A_dot_A = A.dot(A);
    B_dot_B = B.dot(B);
    A_dot_B = A.dot(B);
    A_dot_T = A.dot(T);
    B_dot_T = B.dot(T);

    // t parameterizes ray P,A 
    // u parameterizes ray Q,B 

    float t, u;

    // compute t for the closest point on ray P,A to
    // ray Q,B

    float denom = A_dot_A * B_dot_B - A_dot_B * A_dot_B;

    t = (A_dot_T * B_dot_B - B_dot_T * A_dot_B) / denom;

    // clamp result so t is on the segment P,A

    if ((t < 0) || isnan(t)) t = 0; else if (t > 1) t = 1;

    // find u for point on ray Q,B closest to point at t

    u = (t * A_dot_B - B_dot_T) / B_dot_B;

    // if u is on segment Q,B, t and u correspond to 
    // closest points, otherwise, clamp u, recompute and
    // clamp t 

    if ((u <= 0) || isnan(u)) {

        Y = Q;

        t = A_dot_T / A_dot_A;

        if ((t <= 0) || isnan(t)) {
            X = P;
            VEC = Q - P;
        }
        else if (t >= 1) {
            X = P + A;
            VEC = Q - X;
        }
        else {
            X = P + A * t;
            TMP = T.cross(A);
            VEC = A.cross(TMP);
        }
    }
    else if (u >= 1) {

        Y = Q + B;

        t = (A_dot_B + A_dot_T) / A_dot_A;

        if ((t <= 0) || isnan(t)) {
            X = P;
            VEC = Y - P;
        }
        else if (t >= 1) {
            X = P + A;
            VEC = Y - X;
        }
        else {
            X = P + A * t;
            T = Y - P;
            TMP = T.cross(A);
            VEC = A.cross(TMP);
        }
    }
    else {

        Y = Q + B * u;

        if ((t <= 0) || isnan(t)) {
            X = P;
            TMP = T.cross(B);
            VEC = B.cross(TMP);
        }
        else if (t >= 1) {
            X = P + A;
            T = Q - X;
            TMP = T.cross(B);
            VEC = B.cross(TMP);
        }
        else {
            X = P + A * t;
            VEC = A.cross(B);
            if (VEC.dot(T) < 0) {
                VEC = VEC * -1;
            }
        }
    }
}

//--------------------------------------------------------------------------
// TriDist() 
//
// Computes the closest points on two triangles, and returns the 
// distance between them.
// 
// S and T are the triangles, stored tri[point][dimension].
//
// If the triangles are disjoint, P and Q give the closest points of 
// S and T respectively. However, if the triangles overlap, P and Q 
// are basically a random pair of points from the triangles, not 
// coincident points on the intersection of the triangles, as might 
// be expected.
//--------------------------------------------------------------------------

float TriDist(Eigen::Vector3f& P, Eigen::Vector3f& Q,
    const Eigen::Vector3f S[3], const Eigen::Vector3f T[3])
{
    // Compute vectors along the 6 sides

    Eigen::Vector3f Sv[3], Tv[3];
    Eigen::Vector3f VEC;

    Sv[0] = S[1] - S[0];
    Sv[1] = S[2] - S[1];
    Sv[2] = S[0] - S[2];

    Tv[0] = T[1] - T[0];
    Tv[1] = T[2] - T[1];
    Tv[2] = T[0] - T[2];

    // For each edge pair, the vector connecting the closest points 
    // of the edges defines a slab (parallel planes at head and tail
    // enclose the slab). If we can show that the off-edge vertex of 
    // each triangle is outside of the slab, then the closest points
    // of the edges are the closest points for the triangles.
    // Even if these tests fail, it may be helpful to know the closest
    // points found, and whether the triangles were shown disjoint

    Eigen::Vector3f V;
    Eigen::Vector3f Z;
    Eigen::Vector3f minP, minQ;
    float mindd;
    int shown_disjoint = 0;

    mindd = (S[0] - T[0]).squaredNorm() + 1;  // Set first minimum safely high

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            // Find closest points on edges i & j, plus the 
            // vector (and distance squared) between these points

            SegPoints(VEC, P, Q, S[i], Sv[i], T[j], Tv[j]);

            V = Q - P;
            float dd = V.dot(V);

            // Verify this closest point pair only if the distance 
            // squared is less than the minimum found thus far.

            if (dd <= mindd)
            {
                minP = P;
                minQ = Q;
                mindd = dd;

                Z = S[(i + 2) % 3] - P;
                float a = Z.dot(VEC);
                Z = T[(j + 2) % 3] - Q;
                float b = Z.dot(VEC);

                if ((a <= 0) && (b >= 0)) return sqrt(dd);

                float p = V.dot(VEC);

                if (a < 0) a = 0;
                if (b > 0) b = 0;
                if ((p - a + b) > 0) shown_disjoint = 1;
            }
        }
    }

    // No edge pairs contained the closest points.  
    // either:
    // 1. one of the closest points is a vertex, and the
    //    other point is interior to a face.
    // 2. the triangles are overlapping.
    // 3. an edge of one triangle is parallel to the other's face. If
    //    cases 1 and 2 are not true, then the closest points from the 9
    //    edge pairs checks above can be taken as closest points for the
    //    triangles.
    // 4. possibly, the triangles were degenerate.  When the 
    //    triangle points are nearly colinear or coincident, one 
    //    of above tests might fail even though the edges tested
    //    contain the closest points.

    // First check for case 1

    Eigen::Vector3f Sn;
    float Snl;

    Sn = Sv[0].cross(Sv[1]); // Compute normal to S triangle
    Snl = Sn.dot(Sn);      // Compute square of length of normal

    // If cross product is long enough,

    if (Snl > 1e-15)
    {
        // Get projection lengths of T points

        Eigen::Vector3f Tp;

        V = S[0] - T[0];
        Tp[0] = V.dot(Sn);

        V = S[0] - T[1];
        Tp[1] = V.dot(Sn);

        V = S[0] - T[2];
        Tp[2] = V.dot(Sn);

        // If Sn is a separating direction,
        // find point with smallest projection

        int point = -1;
        if ((Tp[0] > 0) && (Tp[1] > 0) && (Tp[2] > 0))
        {
            if (Tp[0] < Tp[1]) point = 0; else point = 1;
            if (Tp[2] < Tp[point]) point = 2;
        }
        else if ((Tp[0] < 0) && (Tp[1] < 0) && (Tp[2] < 0))
        {
            if (Tp[0] > Tp[1]) point = 0; else point = 1;
            if (Tp[2] > Tp[point]) point = 2;
        }

        // If Sn is a separating direction, 

        if (point >= 0)
        {
            shown_disjoint = 1;

            // Test whether the point found, when projected onto the 
            // other triangle, lies within the face.

            V = T[point] - S[0];
            Z = Sn.cross(Sv[0]);
            if (V.dot(Z) > 0)
            {
                V = T[point] - S[1];
                Z = Sn.cross(Sv[1]);
                if (V.dot(Z) > 0)
                {
                    V = T[point] - S[2];
                    Z = Sn.cross(Sv[2]);
                    if (V.dot(Z) > 0)
                    {
                        // T[point] passed the test - it's a closest point for 
                        // the T triangle; the other point is on the face of S

                        P = T[point] + Sn * Tp[point] / Snl;
                        Q = T[point];
                        return sqrt((P - Q).squaredNorm());
                    }
                }
            }
        }
    }

    float Tnl;
    Eigen::Vector3f Tn;

    Tn = Tv[0].cross(Tv[1]);
    Tnl = Tn.dot(Tn);

    if (Tnl > 1e-15)
    {
        float Sp[3];

        V = T[0] - S[0];
        Sp[0] = V.dot(Tn);

        V = T[0] - S[1];
        Sp[1] = V.dot(Tn);

        V = T[0] - S[2];
        Sp[2] = V.dot(Tn);

        int point = -1;
        if ((Sp[0] > 0) && (Sp[1] > 0) && (Sp[2] > 0))
        {
            if (Sp[0] < Sp[1]) point = 0; else point = 1;
            if (Sp[2] < Sp[point]) point = 2;
        }
        else if ((Sp[0] < 0) && (Sp[1] < 0) && (Sp[2] < 0))
        {
            if (Sp[0] > Sp[1]) point = 0; else point = 1;
            if (Sp[2] > Sp[point]) point = 2;
        }

        if (point >= 0)
        {
            shown_disjoint = 1;

            V = S[point] - T[0];
            Z = Tn.cross(Tv[0]);
            if (V.dot(Z) > 0)
            {
                V = S[point] - T[1];
                Z = Tn.cross(Tv[1]);
                if (V.dot(Z) > 0)
                {
                    V = S[point] - T[2];
                    Z = Tn.cross(Tv[2]);
                    if (V.dot(Z) > 0)
                    {
                        P = S[point];
                        Q = S[point] + Tn * Sp[point] / Tnl;
                        return sqrt((P - Q).squaredNorm());
                    }
                }
            }
        }
    }

    // Case 1 can't be shown.
    // If one of these tests showed the triangles disjoint,
    // we assume case 3 or 4, otherwise we conclude case 2, 
    // that the triangles overlap.

    if (shown_disjoint)
    {
        P = minP;
        Q = minQ;
        return sqrt(mindd);
    }
    else return 0;
}


float Polyhedron::PolyDist(Polyhedron B) {
	float dist = FLT_MAX;
	Eigen::Vector3f p, q;
	for (auto tri_1 : TriangleList) {
		for (auto tri_2 : B.TriangleList) {
            dist = std::min(dist, TriDist(p, q, tri_1->v, tri_2->v));
		}
	}
    return dist;
}
