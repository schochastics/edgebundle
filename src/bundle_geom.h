// Shared geometry, compatibility and subdivision helpers for the force-directed
// bundlers (edge_bundle_force and the divided variant). Header-only; every
// function is inline so multiple translation units can include it.
#ifndef EDGEBUNDLE_BUNDLE_GEOM_H
#define EDGEBUNDLE_BUNDLE_GEOM_H

#include <Rcpp.h>
using namespace Rcpp;

inline double euclidean_distance(NumericVector P, NumericVector Q) {
    return sqrt((P[0] - Q[0]) * (P[0] - Q[0]) + (P[1] - Q[1]) * (P[1] - Q[1]));
}

inline double edge_length(NumericVector P, NumericVector Q, double eps) {
    if ((std::abs(P[0] - Q[0]) < eps) && (std::abs(P[1] - Q[1]) < eps)) {
        return eps;
    }
    return euclidean_distance(P, Q);
}

inline double vector_dot_product(NumericVector P, NumericVector Q) {
    return P[0] * Q[0] + P[1] * Q[1];
}

inline NumericVector edge_as_vector(NumericVector P) {
    NumericVector vec = {P[2] - P[0], P[3] - P[1]};
    return vec;
}

inline NumericVector project_point_on_line(NumericVector p, NumericVector Q) {
    double L = sqrt((Q[2] - Q[0]) * (Q[2] - Q[0]) + (Q[3] - Q[1]) * (Q[3] - Q[1]));
    double r = ((Q[1] - p[1]) * (Q[1] - Q[3]) - (Q[0] - p[0]) * (Q[2] - Q[0])) / (L * L);
    NumericVector vec = {Q[0] + r * (Q[2] - Q[0]), Q[1] + r * (Q[3] - Q[1])};
    return vec;
}

inline double edge_visibility(NumericVector P, NumericVector Q) {
    NumericVector qs = {Q[0], Q[1]};
    NumericVector qt = {Q[2], Q[3]};
    NumericVector I0 = project_point_on_line(qs, P);
    NumericVector I1 = project_point_on_line(qt, P);
    NumericVector midI = {(I0[0] + I1[0]) / 2.0, (I0[1] + I1[1]) / 2.0};
    NumericVector midP = {(P[0] + P[2]) / 2.0, (P[1] + P[3]) / 2.0};
    double tmp = 1.0 - 2.0 * euclidean_distance(midP, midI) / euclidean_distance(I0, I1);
    return tmp > 0 ? tmp : 0;
}

inline double compute_divided_edge_length(NumericMatrix emat) {
    int segments = emat.rows() - 1;
    double length = 0.0;
    for (int i = 0; i < segments; ++i) {
        length += euclidean_distance(emat(i, _), emat(i + 1, _));
    }
    return length;
}

// Resample every edge to P interior subdivision points (P + 2 points total).
inline List update_edge_divisions(List elist, int P) {
    for (int e_idx = 0; e_idx < elist.length(); ++e_idx) {
        NumericMatrix emat = elist[e_idx];
        if (P == 1) {
            NumericMatrix emat_new(3, 2);
            emat_new(0, 0) = emat(0, 0);
            emat_new(0, 1) = emat(0, 1);
            emat_new(1, 0) = (emat(0, 0) + emat(1, 0)) / 2.0;
            emat_new(1, 1) = (emat(0, 1) + emat(1, 1)) / 2.0;
            emat_new(2, 0) = emat(1, 0);
            emat_new(2, 1) = emat(1, 1);
            elist[e_idx] = emat_new;
        } else {
            double divided_edge_length = compute_divided_edge_length(emat);
            double segment_length = divided_edge_length / (P + 1);
            double current_segment_length = segment_length;
            NumericMatrix emat_new(P + 2, 2);
            emat_new(0, _) = emat(0, _);
            emat_new((emat_new.rows() - 1), _) = emat((emat.rows() - 1), _);
            int cur = 1;
            for (int i = 1; i < emat.rows(); ++i) {
                double old_segment_length = euclidean_distance(emat(i - 1, _), emat(i, _));
                while (old_segment_length > current_segment_length) {
                    double percent_position = current_segment_length / old_segment_length;
                    double new_x = emat(i - 1, 0) + percent_position * (emat(i, 0) - emat(i - 1, 0));
                    double new_y = emat(i - 1, 1) + percent_position * (emat(i, 1) - emat(i - 1, 1));
                    emat_new(cur, 0) = new_x;
                    emat_new(cur, 1) = new_y;
                    cur += 1;
                    old_segment_length -= current_segment_length;
                    current_segment_length = segment_length;
                }
                current_segment_length -= old_segment_length;
            }
            elist[e_idx] = emat_new;
        }
    }
    return elist;
}

inline double angle_compatibility(NumericVector P, NumericVector Q) {
    NumericVector P_source = {P[0], P[1]}, P_target = {P[2], P[3]};
    NumericVector Q_source = {Q[0], Q[1]}, Q_target = {Q[2], Q[3]};
    double dot_PQ = vector_dot_product(edge_as_vector(P), edge_as_vector(Q));
    double euc_PQ = euclidean_distance(P_source, P_target) * euclidean_distance(Q_source, Q_target);
    return std::abs(dot_PQ / euc_PQ);
}

inline double scale_compatibility(NumericVector P, NumericVector Q) {
    NumericVector P_source = {P[0], P[1]}, P_target = {P[2], P[3]};
    NumericVector Q_source = {Q[0], Q[1]}, Q_target = {Q[2], Q[3]};
    double euc_P = euclidean_distance(P_source, P_target);
    double euc_Q = euclidean_distance(Q_source, Q_target);
    double lavg = (euc_P + euc_Q) / 2.0;
    return 2.0 / (lavg / std::min(euc_P, euc_Q) + std::max(euc_P, euc_Q) / lavg);
}

inline double position_compatibility(NumericVector P, NumericVector Q) {
    NumericVector P_source = {P[0], P[1]}, P_target = {P[2], P[3]};
    NumericVector Q_source = {Q[0], Q[1]}, Q_target = {Q[2], Q[3]};
    double euc_P = euclidean_distance(P_source, P_target);
    double euc_Q = euclidean_distance(Q_source, Q_target);
    double lavg = (euc_P + euc_Q) / 2.0;
    NumericVector midP = {(P_source[0] + P_target[0]) / 2.0, (P_source[1] + P_target[1]) / 2.0};
    NumericVector midQ = {(Q_source[0] + Q_target[0]) / 2.0, (Q_source[1] + Q_target[1]) / 2.0};
    return lavg / (lavg + euclidean_distance(midP, midQ));
}

inline double visibility_compatibility(NumericVector P, NumericVector Q) {
    return std::min(edge_visibility(P, Q), edge_visibility(Q, P));
}

inline double compatibility_score(NumericVector P, NumericVector Q) {
    return angle_compatibility(P, Q) * scale_compatibility(P, Q) *
        position_compatibility(P, Q) * visibility_compatibility(P, Q);
}

inline bool are_compatible(NumericVector P, NumericVector Q, double compatibility_threshold) {
    return compatibility_score(P, Q) >= compatibility_threshold;
}

// Symmetric list of compatible edge indices per edge.
inline List compute_compatibility_lists(NumericMatrix edges_xy, double compatibility_threshold) {
    int m = edges_xy.rows();
    List elist_comp(m);
    for (int e = 0; e < (m - 1); ++e) {
        NumericVector P = edges_xy(e, _);
        for (int oe = (e + 1); oe < m; ++oe) {
            NumericVector Q = edges_xy(oe, _);
            if (are_compatible(P, Q, compatibility_threshold)) {
                if (elist_comp[e] == R_NilValue) {
                    IntegerVector ecomp = {oe};
                    elist_comp[e] = ecomp;
                } else {
                    IntegerVector ecomp = elist_comp[e];
                    ecomp.push_back(oe);
                    elist_comp[e] = ecomp;
                }
                if (elist_comp[oe] == R_NilValue) {
                    IntegerVector oecomp = {e};
                    elist_comp[oe] = oecomp;
                } else {
                    IntegerVector oecomp = elist_comp[oe];
                    oecomp.push_back(e);
                    elist_comp[oe] = oecomp;
                }
            }
        }
    }
    return elist_comp;
}

inline NumericVector apply_spring_force(List elist, int e_idx, int i, double kP) {
    NumericMatrix emat = elist[e_idx];
    NumericVector prec = emat(i - 1, _);
    NumericVector succ = emat(i + 1, _);
    NumericVector crnt = emat(i, _);
    double x = (prec[0] - crnt[0] + succ[0] - crnt[0]) * kP;
    double y = (prec[1] - crnt[1] + succ[1] - crnt[1]) * kP;
    return {x, y};
}

// Build the initial subdivision list (source + target point) from endpoint rows.
inline List init_edge_list(NumericMatrix edges_xy) {
    int m = edges_xy.rows();
    List elist(m);
    for (int e = 0; e < m; ++e) {
        NumericMatrix emat(2, 2);
        emat(0, 0) = edges_xy(e, 0);
        emat(0, 1) = edges_xy(e, 1);
        emat(1, 0) = edges_xy(e, 2);
        emat(1, 1) = edges_xy(e, 3);
        elist[e] = emat;
    }
    return elist;
}

#endif
