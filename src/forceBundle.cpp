// Re-implementation of https://github.com/upphiminn/d3.ForceBundle which implements
// Holten, Danny, and Jarke J. Van Wijk. "Force-Directed Edge Bundling for Graph Visualization."
// Computer Graphics Forum (Blackwell Publishing Ltd) 28, no. 3 (2009): 983-990.
#include <Rcpp.h>
#include "bundle_geom.h"
using namespace Rcpp;

NumericVector apply_electrostatic_force(List elist, List elist_comp, int e_idx, int i, double eps) {
    NumericVector sum_of_forces(2);
    if (elist_comp[e_idx] == R_NilValue) {
        return sum_of_forces;
    }
    IntegerVector ecomps = elist_comp[e_idx];
    NumericMatrix emat = elist[e_idx];

    for (int oe = 0; oe < ecomps.length(); ++oe) {
        NumericMatrix oemat = elist[ecomps[oe]];
        NumericVector force = {oemat(i, 0) - emat(i, 0), oemat(i, 1) - emat(i, 1)};
        if ((std::abs(force[0]) > eps) || (std::abs(force[1]) > eps)) {
            double diff = std::pow(euclidean_distance(oemat(i, _), emat(i, _)), -1.0);
            sum_of_forces[0] += force[0] * diff;
            sum_of_forces[1] += force[1] * diff;
        }
    }
    return sum_of_forces;
}

NumericMatrix apply_resulting_forces_on_subdivision_points(List elist, List elist_comp,
                                                           int e_idx, int P, double S, double K, double eps) {
    NumericMatrix emat = elist[e_idx];
    double kP = K / (edge_length(emat(0, _), emat(P + 1, _), eps) * (P + 1));

    NumericMatrix resulting_forces_for_subdivision_points(P + 2, 2);
    for (int i = 1; i < (P + 1); ++i) {
        NumericVector spring_force = apply_spring_force(elist, e_idx, i, kP);
        NumericVector electrostatic_force = apply_electrostatic_force(elist, elist_comp, e_idx, i, eps);
        resulting_forces_for_subdivision_points(i, 0) = S * (spring_force[0] + electrostatic_force[0]);
        resulting_forces_for_subdivision_points(i, 1) = S * (spring_force[1] + electrostatic_force[1]);
    }
    return resulting_forces_for_subdivision_points;
}

// [[Rcpp::export]]
List force_bundle_iter(NumericMatrix edges_xy,
                       double K, int C, int P, int P_rate,
                       double S, int I, double I_rate,
                       double compatibility_threshold, double eps) {
    int m = edges_xy.rows();

    List elist = init_edge_list(edges_xy);
    elist = update_edge_divisions(elist, P);
    List elist_comp = compute_compatibility_lists(edges_xy, compatibility_threshold);

    for (int cycle = 0; cycle < C; ++cycle) {
        for (int iteration = 0; iteration < I; ++iteration) {
            List forces(m);
            for (int e = 0; e < m; ++e) {
                forces[e] = apply_resulting_forces_on_subdivision_points(elist, elist_comp, e, P, S, K, eps);
            }
            for (int e = 0; e < m; ++e) {
                NumericMatrix emat = elist[e];
                NumericMatrix fmat = forces[e];
                for (int i = 0; i < (P + 1); ++i) {
                    emat(i, 0) += fmat(i, 0);
                    emat(i, 1) += fmat(i, 1);
                }
                elist[e] = emat;
            }
        }
        if (cycle != (C - 1)) {
            S = S / 2.0;
            P = P * P_rate;
            I = I * I_rate;
            elist = update_edge_divisions(elist, P);
        }
    }
    return elist;
}
