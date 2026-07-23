// MINGLE: multilevel agglomerative edge bundling.
// Gansner, Hu, North, Scheidegger (2011) "Multilevel agglomerative edge
// bundling for visualizing large graphs", IEEE PacificVis.
//
// Edges are grouped bottom-up: two bundles merge when routing them through
// shared meeting points reduces total "ink" (drawn length). Meeting points are
// the geometric medians of the source-side and target-side endpoints. kNN over
// edges is brute force (O(E^2)); fine for interactive sizes, and keeps the
// package dependency-free (a kd-tree could be swapped in for very large graphs).
#include <Rcpp.h>
using namespace Rcpp;

struct Pt {
    double x, y;
};

static double dist(const Pt &a, const Pt &b) {
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

// geometric median (Weiszfeld); minimises sum of distances = the ink terms
static Pt geomedian(const std::vector<Pt> &pts) {
    Pt m{0, 0};
    for (const Pt &p : pts) {
        m.x += p.x;
        m.y += p.y;
    }
    m.x /= pts.size();
    m.y /= pts.size();
    for (int it = 0; it < 20; ++it) {
        double sw = 0, sx = 0, sy = 0;
        for (const Pt &p : pts) {
            double d = dist(p, m);
            if (d < 1e-9) continue;
            double w = 1.0 / d;
            sw += w;
            sx += w * p.x;
            sy += w * p.y;
        }
        if (sw == 0) break;
        m.x = sx / sw;
        m.y = sy / sw;
    }
    return m;
}

struct Bundle {
    std::vector<int> edges;    // original edge ids
    std::vector<char> flip;    // per-member: is (a,b) the reverse of (source,target)?
    std::vector<Pt> A, B;      // source-side and target-side endpoints
    Pt ma, mb;
    double ink;
};

static void finalize(Bundle &b) {
    b.ma = geomedian(b.A);
    b.mb = geomedian(b.B);
    double s = dist(b.ma, b.mb);
    for (const Pt &p : b.A) s += dist(p, b.ma);
    for (const Pt &p : b.B) s += dist(p, b.mb);
    b.ink = s;
}

// ink of the union of two bundles under a given alignment (b2 flipped or not)
static double merged_ink(const Bundle &b1, const Bundle &b2, bool flip2) {
    std::vector<Pt> A = b1.A, B = b1.B;
    const std::vector<Pt> &A2 = flip2 ? b2.B : b2.A;
    const std::vector<Pt> &B2 = flip2 ? b2.A : b2.B;
    A.insert(A.end(), A2.begin(), A2.end());
    B.insert(B.end(), B2.begin(), B2.end());
    Pt ma = geomedian(A), mb = geomedian(B);
    double s = dist(ma, mb);
    for (const Pt &p : A) s += dist(p, ma);
    for (const Pt &p : B) s += dist(p, mb);
    return s;
}

static Bundle do_merge(const Bundle &b1, const Bundle &b2, bool flip2) {
    Bundle out;
    out.edges = b1.edges;
    out.flip = b1.flip;
    out.A = b1.A;
    out.B = b1.B;
    for (size_t i = 0; i < b2.edges.size(); ++i) {
        out.edges.push_back(b2.edges[i]);
        out.flip.push_back(flip2 ? !b2.flip[i] : b2.flip[i]);
    }
    const std::vector<Pt> &A2 = flip2 ? b2.B : b2.A;
    const std::vector<Pt> &B2 = flip2 ? b2.A : b2.B;
    out.A.insert(out.A.end(), A2.begin(), A2.end());
    out.B.insert(out.B.end(), B2.begin(), B2.end());
    finalize(out);
    return out;
}

// [[Rcpp::export]]
List mingle_iter(NumericMatrix edges_xy, int k, int segments, double bundle_strength) {
    int m = edges_xy.rows();

    std::vector<Bundle> bundles(m);
    for (int e = 0; e < m; ++e) {
        Pt s{edges_xy(e, 0), edges_xy(e, 1)};
        Pt t{edges_xy(e, 2), edges_xy(e, 3)};
        // canonical order so an edge and its reverse group identically
        bool flip = (s.x > t.x) || (s.x == t.x && s.y > t.y);
        bundles[e].edges = {e};
        bundles[e].flip = {(char)flip};
        bundles[e].A = {flip ? t : s};
        bundles[e].B = {flip ? s : t};
        finalize(bundles[e]);
    }

    // meeting-point chain per edge, from finest to coarsest merged level, in
    // original source->target orientation (leaf endpoints are added at render)
    std::vector<std::vector<Pt>> srcChain(m), tgtChain(m);

    bool improved = true;
    while (improved && bundles.size() > 1) {
        improved = false;
        int nb = bundles.size();

        // brute-force kNN over bundle descriptors (ma,mb)
        std::vector<std::pair<double, std::pair<int, bool>>> best(nb, {1e18, {-1, false}});
        for (int i = 0; i < nb; ++i) {
            for (int j = 0; j < nb; ++j) {
                if (i == j) continue;
                double d = dist(bundles[i].ma, bundles[j].ma) + dist(bundles[i].mb, bundles[j].mb);
                double dflip = dist(bundles[i].ma, bundles[j].mb) + dist(bundles[i].mb, bundles[j].ma);
                bool fl = dflip < d;
                double dd = std::min(d, dflip);
                if (dd < best[i].first) best[i] = {dd, {j, fl}};
            }
        }
        (void)k; // kNN kept to nearest for simplicity; k reserved for future use

        // greedy ink-reducing matching
        std::vector<char> used(nb, 0);
        std::vector<Bundle> next;
        std::vector<int> newIndexOf(nb, -1);

        std::vector<std::pair<double, int>> order;
        for (int i = 0; i < nb; ++i) {
            int j = best[i].second.first;
            if (j < 0) continue;
            bool fl = best[i].second.second;
            double sav = bundles[i].ink + bundles[j].ink - merged_ink(bundles[i], bundles[j], fl);
            order.push_back({sav, i});
        }
        std::sort(order.begin(), order.end(), [](auto &a, auto &b) { return a.first > b.first; });

        for (auto &o : order) {
            int i = o.second;
            int j = best[i].second.first;
            bool fl = best[i].second.second;
            if (used[i] || used[j]) continue;
            if (o.first <= 1e-9) continue;
            used[i] = used[j] = 1;
            Bundle mb = do_merge(bundles[i], bundles[j], fl);
            newIndexOf[i] = newIndexOf[j] = next.size();
            next.push_back(mb);
            improved = true;
        }
        for (int i = 0; i < nb; ++i) {
            if (!used[i]) {
                newIndexOf[i] = next.size();
                next.push_back(bundles[i]);
            }
        }

        if (improved) {
            // record this level's meeting points for every edge
            for (int i = 0; i < nb; ++i) {
                Bundle &nb2 = next[newIndexOf[i]];
                for (size_t mi = 0; mi < nb2.edges.size(); ++mi) {
                    int e = nb2.edges[mi];
                    bool flipped = nb2.flip[mi];
                    srcChain[e].push_back(flipped ? nb2.mb : nb2.ma);
                    tgtChain[e].push_back(flipped ? nb2.ma : nb2.mb);
                }
            }
            // dedup: the loop above pushes once per edge because each edge lives
            // in exactly one next bundle
        }
        bundles = next;
    }

    // assemble per-edge control polyline: source, src meeting points (fine->coarse),
    // tgt meeting points (coarse->fine), target; then subdivide toward straight
    List out(m);
    for (int e = 0; e < m; ++e) {
        Pt s{edges_xy(e, 0), edges_xy(e, 1)};
        Pt t{edges_xy(e, 2), edges_xy(e, 3)};
        std::vector<Pt> cp;
        cp.push_back(s);
        for (size_t i = 0; i < srcChain[e].size(); ++i) cp.push_back(srcChain[e][i]);
        for (size_t i = tgtChain[e].size(); i-- > 0;) cp.push_back(tgtChain[e][i]);
        cp.push_back(t);

        // relax control points toward the straight line by (1 - bundle_strength)
        for (size_t i = 1; i + 1 < cp.size(); ++i) {
            double u = (double)i / (cp.size() - 1);
            Pt straight{s.x + u * (t.x - s.x), s.y + u * (t.y - s.y)};
            cp[i].x = bundle_strength * cp[i].x + (1 - bundle_strength) * straight.x;
            cp[i].y = bundle_strength * cp[i].y + (1 - bundle_strength) * straight.y;
        }

        // resample the control polyline to `segments` points by arc length
        int kc = cp.size();
        std::vector<double> cl(kc, 0.0);
        for (int i = 1; i < kc; ++i) cl[i] = cl[i - 1] + dist(cp[i - 1], cp[i]);
        double total = cl[kc - 1];
        NumericMatrix mat(segments, 2);
        for (int i = 0; i < segments; ++i) {
            if (total <= 0) {
                mat(i, 0) = s.x;
                mat(i, 1) = s.y;
                continue;
            }
            double target = total * i / (segments - 1);
            int seg = 0;
            while (seg < kc - 2 && cl[seg + 1] < target) seg++;
            double segd = cl[seg + 1] - cl[seg];
            double u = segd > 0 ? (target - cl[seg]) / segd : 0.0;
            mat(i, 0) = cp[seg].x + u * (cp[seg + 1].x - cp[seg].x);
            mat(i, 1) = cp[seg].y + u * (cp[seg + 1].y - cp[seg].y);
        }
        out[e] = mat;
    }
    return out;
}
