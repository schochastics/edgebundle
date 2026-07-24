// MINGLE: multilevel agglomerative edge bundling.
// Gansner, Hu, North, Scheidegger (2011) "Multilevel agglomerative edge
// bundling for visualizing large graphs", IEEE PacificVis.
//
// Edges are grouped bottom-up: two bundles merge when routing them through
// shared meeting points reduces total "ink" (drawn length). Meeting points are
// the geometric medians of the source-side and target-side endpoints.
//
// A kNN proximity graph over the edges (points in canonical 4-D endpoint space)
// is built once with a kd-tree (vendored nanoflann), giving O(E log E); the
// graph is then coarsened combinatorially level by level.
#include <Rcpp.h>
#include "nanoflann.hpp"
#include <set>
#include <array>
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

// nanoflann point-cloud adaptor over canonical 4-D edge coordinates
struct EdgeCloud {
    std::vector<std::array<double, 4>> pts;
    inline size_t kdtree_get_point_count() const { return pts.size(); }
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const { return pts[idx][dim]; }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX &) const { return false; }
};

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

    // one-time kNN proximity graph in canonical 4-D endpoint space (kd-tree)
    std::vector<std::vector<int>> adj(m);
    if (m > 1) {
        EdgeCloud cloud;
        cloud.pts.resize(m);
        for (int e = 0; e < m; ++e) {
            cloud.pts[e] = {bundles[e].A[0].x, bundles[e].A[0].y,
                            bundles[e].B[0].x, bundles[e].B[0].y};
        }
        typedef nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<double, EdgeCloud>, EdgeCloud, 4>
            KDTree;
        KDTree tree(4, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
        tree.buildIndex();

        size_t kq = std::min((size_t)k + 1, (size_t)m);
        std::vector<std::set<int>> adjset(m);
        std::vector<size_t> nn_idx(kq);
        std::vector<double> nn_d2(kq);
        for (int e = 0; e < m; ++e) {
            nanoflann::KNNResultSet<double> res(kq);
            res.init(nn_idx.data(), nn_d2.data());
            tree.findNeighbors(res, cloud.pts[e].data(), nanoflann::SearchParameters());
            for (size_t t = 0; t < kq; ++t) {
                int j = (int)nn_idx[t];
                if (j == e) continue;
                adjset[e].insert(j);
                adjset[j].insert(e); // keep the candidate graph symmetric
            }
        }
        for (int e = 0; e < m; ++e) adj[e].assign(adjset[e].begin(), adjset[e].end());
    }

    // meeting-point chain per edge, from finest to coarsest merged level, in
    // original source->target orientation (leaf endpoints are added at render)
    std::vector<std::vector<Pt>> srcChain(m), tgtChain(m);

    bool improved = true;
    while (improved && bundles.size() > 1) {
        improved = false;
        int nb = bundles.size();

        // nearest candidate neighbour (by meeting-point distance) per bundle
        std::vector<int> bestj(nb, -1);
        std::vector<char> bestflip(nb, 0);
        std::vector<double> bestsav(nb, 0.0);
        for (int i = 0; i < nb; ++i) {
            double bd = 1e18;
            int bj = -1;
            bool bf = false;
            for (int j : adj[i]) {
                double d = dist(bundles[i].ma, bundles[j].ma) + dist(bundles[i].mb, bundles[j].mb);
                double dflip = dist(bundles[i].ma, bundles[j].mb) + dist(bundles[i].mb, bundles[j].ma);
                bool fl = dflip < d;
                double dd = std::min(d, dflip);
                if (dd < bd) {
                    bd = dd;
                    bj = j;
                    bf = fl;
                }
            }
            if (bj >= 0) {
                bestj[i] = bj;
                bestflip[i] = bf;
                bestsav[i] = bundles[i].ink + bundles[bj].ink - merged_ink(bundles[i], bundles[bj], bf);
            }
        }

        // greedy ink-reducing matching, largest saving first
        std::vector<int> order;
        for (int i = 0; i < nb; ++i) {
            if (bestj[i] >= 0) order.push_back(i);
        }
        std::sort(order.begin(), order.end(), [&](int a, int b) { return bestsav[a] > bestsav[b]; });

        std::vector<char> used(nb, 0);
        std::vector<Bundle> next;
        std::vector<int> newIndexOf(nb, -1);
        for (int i : order) {
            int j = bestj[i];
            if (used[i] || used[j]) continue;
            if (bestsav[i] <= 1e-9) continue;
            used[i] = used[j] = 1;
            newIndexOf[i] = newIndexOf[j] = next.size();
            next.push_back(do_merge(bundles[i], bundles[j], bestflip[i]));
            improved = true;
        }
        for (int i = 0; i < nb; ++i) {
            if (!used[i]) {
                newIndexOf[i] = next.size();
                next.push_back(bundles[i]);
            }
        }

        if (improved) {
            for (int i = 0; i < nb; ++i) {
                Bundle &nb2 = next[newIndexOf[i]];
                for (size_t mi = 0; mi < nb2.edges.size(); ++mi) {
                    int e = nb2.edges[mi];
                    bool flipped = nb2.flip[mi];
                    srcChain[e].push_back(flipped ? nb2.mb : nb2.ma);
                    tgtChain[e].push_back(flipped ? nb2.ma : nb2.mb);
                }
            }
            // contract the candidate graph onto the coarse bundles
            std::vector<std::set<int>> nadj(next.size());
            for (int i = 0; i < nb; ++i) {
                for (int j : adj[i]) {
                    int a = newIndexOf[i], b = newIndexOf[j];
                    if (a != b) {
                        nadj[a].insert(b);
                        nadj[b].insert(a);
                    }
                }
            }
            adj.assign(next.size(), {});
            for (size_t i = 0; i < next.size(); ++i) adj[i].assign(nadj[i].begin(), nadj[i].end());
        }
        bundles = next;
    }

    // assemble per-edge control polyline: source, src meeting points (fine->coarse),
    // tgt meeting points (coarse->fine), target; then relax toward the straight line
    List out(m);
    for (int e = 0; e < m; ++e) {
        Pt s{edges_xy(e, 0), edges_xy(e, 1)};
        Pt t{edges_xy(e, 2), edges_xy(e, 3)};
        std::vector<Pt> cp;
        cp.push_back(s);
        for (size_t i = 0; i < srcChain[e].size(); ++i) cp.push_back(srcChain[e][i]);
        for (size_t i = tgtChain[e].size(); i-- > 0;) cp.push_back(tgtChain[e][i]);
        cp.push_back(t);

        for (size_t i = 1; i + 1 < cp.size(); ++i) {
            double u = (double)i / (cp.size() - 1);
            Pt straight{s.x + u * (t.x - s.x), s.y + u * (t.y - s.y)};
            cp[i].x = bundle_strength * cp[i].x + (1 - bundle_strength) * straight.x;
            cp[i].y = bundle_strength * cp[i].y + (1 - bundle_strength) * straight.y;
        }

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
