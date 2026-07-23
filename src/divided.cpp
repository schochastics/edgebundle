// Divided edge bundling for directed graphs.
// Port of the reference implementation kakearney/divedgebundle (MATLAB), which
// replicates: Selassie, Heller, Heer (2011) "Divided Edge Bundling for
// Directional Network Data", IEEE TVCG 17:2354-2363.
#include <Rcpp.h>
#include "bundle_geom.h"
using namespace Rcpp;

typedef std::vector<std::vector<double>> Mat;

static double connectivity_compat(int ip, int iq, IntegerMatrix el, NumericMatrix node_dist) {
    int a1 = el(ip, 0) - 1, a2 = el(ip, 1) - 1;
    int b1 = el(iq, 0) - 1, b2 = el(iq, 1) - 1;
    double d = std::min(std::min(node_dist(a1, b1), node_dist(a1, b2)),
                        std::min(node_dist(a2, b1), node_dist(a2, b2)));
    return 1.0 / (1.0 + d);
}

// Accelerations on the interior control points; endpoints are held fixed.
static void compute_forces(const Mat &X, const Mat &Y, Mat &AX, Mat &AY,
                           const Mat &compat, const std::vector<std::vector<char>> &issame,
                           const NumericVector &weight, double ks, double kC,
                           double l, double s, int nc, int m) {
    // lane-offset target locations for every interior point (depend on current geometry)
    Mat MJx(m, std::vector<double>(nc, 0.0)), MJy(m, std::vector<double>(nc, 0.0));
    for (int e = 0; e < m; ++e) {
        for (int i = 1; i < nc - 1; ++i) {
            double Tx = X[e][i + 1] - X[e][i - 1];
            double Ty = Y[e][i + 1] - Y[e][i - 1];
            double Tn = std::sqrt(Tx * Tx + Ty * Ty);
            if (Tn > 0) {
                MJx[e][i] = X[e][i] + l * (-Ty / Tn);
                MJy[e][i] = Y[e][i] + l * (Tx / Tn);
            } else {
                MJx[e][i] = X[e][i];
                MJy[e][i] = Y[e][i];
            }
        }
    }

    for (int e = 0; e < m; ++e) {
        AX[e][0] = AY[e][0] = 0.0;
        AX[e][nc - 1] = AY[e][nc - 1] = 0.0;
        for (int i = 1; i < nc - 1; ++i) {
            double px = X[e][i], py = Y[e][i];
            double fsx = ks * nc * ((X[e][i - 1] - px) + (X[e][i + 1] - px)) * weight[e];
            double fsy = ks * nc * ((Y[e][i - 1] - py) + (Y[e][i + 1] - py)) * weight[e];
            double fcx = 0.0, fcy = 0.0;
            for (int q = 0; q < m; ++q) {
                double tx, ty;
                if (issame[e][q]) {
                    tx = X[q][i];
                    ty = Y[q][i];
                } else {
                    int mir = nc - 1 - i;
                    tx = MJx[q][mir];
                    ty = MJy[q][mir];
                }
                double vx = px - tx, vy = py - ty;
                double r = std::sqrt(vx * vx + vy * vy);
                if (r == 0) continue;
                double fcr = -s * kC * r / (M_PI * nc * (s * s + r * r) * (s * s + r * r));
                double wgt = compat[e][q] * weight[q];
                fcx += (vx / r) * fcr * wgt;
                fcy += (vy / r) * fcr * wgt;
            }
            AX[e][i] = fsx + fcx;
            AY[e][i] = fsy + fcy;
        }
    }
}

// [[Rcpp::export]]
List divided_bundle_iter(NumericMatrix edges_xy, IntegerMatrix el, NumericMatrix node_dist,
                         NumericVector weight, double ks, double kC, double l, double s,
                         double f, double dt, int npass, int nstep, bool use_connectivity) {
    int m = edges_xy.rows();
    kC = kC / std::sqrt((double)m);

    Mat compat(m, std::vector<double>(m, 0.0));
    std::vector<std::vector<char>> issame(m, std::vector<char>(m, 0));
    for (int ip = 0; ip < m; ++ip) {
        NumericVector P = edges_xy(ip, _);
        NumericVector pv = edge_as_vector(P);
        for (int iq = 0; iq < m; ++iq) {
            NumericVector Q = edges_xy(iq, _);
            NumericVector qv = edge_as_vector(Q);
            issame[ip][iq] = (pv[0] * qv[0] + pv[1] * qv[1]) > 0 ? 1 : 0;
            double c = angle_compatibility(P, Q) * scale_compatibility(P, Q) *
                position_compatibility(P, Q) * visibility_compatibility(P, Q);
            if (use_connectivity) c *= connectivity_compat(ip, iq, el, node_dist);
            compat[ip][iq] = c;
        }
    }

    int nc = 2;
    Mat X(m, std::vector<double>(nc)), Y(m, std::vector<double>(nc));
    for (int e = 0; e < m; ++e) {
        X[e][0] = edges_xy(e, 0);
        Y[e][0] = edges_xy(e, 1);
        X[e][1] = edges_xy(e, 2);
        Y[e][1] = edges_xy(e, 3);
    }

    for (int pass = 0; pass < npass; ++pass) {
        int nc_new = nc * 2 - 1;
        Mat Xn(m, std::vector<double>(nc_new)), Yn(m, std::vector<double>(nc_new));
        for (int e = 0; e < m; ++e) {
            for (int i = 0; i < nc; ++i) {
                Xn[e][2 * i] = X[e][i];
                Yn[e][2 * i] = Y[e][i];
            }
            for (int i = 0; i < nc - 1; ++i) {
                Xn[e][2 * i + 1] = (X[e][i] + X[e][i + 1]) / 2.0;
                Yn[e][2 * i + 1] = (Y[e][i] + Y[e][i + 1]) / 2.0;
            }
        }
        X = Xn;
        Y = Yn;
        nc = nc_new;

        Mat VX(m, std::vector<double>(nc, 0.0)), VY(m, std::vector<double>(nc, 0.0));
        Mat AX(m, std::vector<double>(nc, 0.0)), AY(m, std::vector<double>(nc, 0.0));
        for (int step = 0; step < nstep - 1; ++step) {
            for (int e = 0; e < m; ++e) {
                for (int i = 1; i < nc - 1; ++i) {
                    VX[e][i] = (VX[e][i] + AX[e][i] * dt / 2.0) * f;
                    VY[e][i] = (VY[e][i] + AY[e][i] * dt / 2.0) * f;
                    X[e][i] += VX[e][i] * dt;
                    Y[e][i] += VY[e][i] * dt;
                }
            }
            compute_forces(X, Y, AX, AY, compat, issame, weight, ks, kC, l, s, nc, m);
            for (int e = 0; e < m; ++e) {
                for (int i = 1; i < nc - 1; ++i) {
                    VX[e][i] += AX[e][i] * dt / 2.0;
                    VY[e][i] += AY[e][i] * dt / 2.0;
                }
            }
        }
        dt /= 2.0;
    }

    List out(m);
    for (int e = 0; e < m; ++e) {
        NumericMatrix mat(nc, 2);
        for (int i = 0; i < nc; ++i) {
            mat(i, 0) = X[e][i];
            mat(i, 1) = Y[e][i];
        }
        out[e] = mat;
    }
    return out;
}
