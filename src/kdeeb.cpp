// Kernel-density-estimation edge bundling (KDEEB).
// Hurter, Ersoy, Telea (2012) "Graph Bundling by Kernel Density Estimation".
// This is the algorithm behind datashader's "hammer" bundling; implemented
// natively here so no Python/datashader dependency is needed.
#include <Rcpp.h>
using namespace Rcpp;

typedef std::vector<std::vector<double>> Mat;

static void gaussian_blur(std::vector<double> &grid, int n, double sigma) {
    if (sigma < 0.5) return;
    int r = (int)std::ceil(3.0 * sigma);
    std::vector<double> kern(2 * r + 1);
    double ksum = 0.0;
    for (int i = -r; i <= r; ++i) {
        kern[i + r] = std::exp(-(i * i) / (2.0 * sigma * sigma));
        ksum += kern[i + r];
    }
    for (double &k : kern) k /= ksum;

    std::vector<double> tmp(n * n, 0.0);
    // horizontal
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            double acc = 0.0;
            for (int i = -r; i <= r; ++i) {
                int xx = std::min(std::max(x + i, 0), n - 1);
                acc += grid[y * n + xx] * kern[i + r];
            }
            tmp[y * n + x] = acc;
        }
    }
    // vertical
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            double acc = 0.0;
            for (int i = -r; i <= r; ++i) {
                int yy = std::min(std::max(y + i, 0), n - 1);
                acc += tmp[yy * n + x] * kern[i + r];
            }
            grid[y * n + x] = acc;
        }
    }
}

// bilinear sample of a grid at continuous cell coordinates
static double sample(const std::vector<double> &grid, int n, double cx, double cy) {
    int x0 = (int)std::floor(cx), y0 = (int)std::floor(cy);
    double fx = cx - x0, fy = cy - y0;
    int x1 = std::min(x0 + 1, n - 1), y1 = std::min(y0 + 1, n - 1);
    x0 = std::min(std::max(x0, 0), n - 1);
    y0 = std::min(std::max(y0, 0), n - 1);
    return grid[y0 * n + x0] * (1 - fx) * (1 - fy) + grid[y0 * n + x1] * fx * (1 - fy) +
        grid[y1 * n + x0] * (1 - fx) * fy + grid[y1 * n + x1] * fx * fy;
}

// resample a polyline to `npoints` equally spaced points (endpoints preserved)
static void resample(std::vector<double> &px, std::vector<double> &py, int npoints) {
    int k = px.size();
    std::vector<double> cl(k, 0.0);
    for (int i = 1; i < k; ++i) {
        double dx = px[i] - px[i - 1], dy = py[i] - py[i - 1];
        cl[i] = cl[i - 1] + std::sqrt(dx * dx + dy * dy);
    }
    double total = cl[k - 1];
    std::vector<double> nx(npoints), ny(npoints);
    if (total <= 0) {
        for (int i = 0; i < npoints; ++i) {
            nx[i] = px[0];
            ny[i] = py[0];
        }
    } else {
        int seg = 0;
        for (int i = 0; i < npoints; ++i) {
            double target = total * i / (npoints - 1);
            while (seg < k - 2 && cl[seg + 1] < target) seg++;
            double segd = cl[seg + 1] - cl[seg];
            double t = segd > 0 ? (target - cl[seg]) / segd : 0.0;
            nx[i] = px[seg] + t * (px[seg + 1] - px[seg]);
            ny[i] = py[seg] + t * (py[seg + 1] - py[seg]);
        }
    }
    px = nx;
    py = ny;
}

// [[Rcpp::export]]
List kdeeb_iter(NumericMatrix edges_xy, int npoints, int niter, double bw,
                double decay, int grid, double step, int smooth_passes) {
    int m = edges_xy.rows();

    double xmin = edges_xy(0, 0), xmax = xmin, ymin = edges_xy(0, 1), ymax = ymin;
    for (int e = 0; e < m; ++e) {
        for (int c = 0; c < 4; c += 2) {
            xmin = std::min(xmin, edges_xy(e, c));
            xmax = std::max(xmax, edges_xy(e, c));
            ymin = std::min(ymin, edges_xy(e, c + 1));
            ymax = std::max(ymax, edges_xy(e, c + 1));
        }
    }
    double span = std::max(xmax - xmin, ymax - ymin);
    if (span <= 0) span = 1.0;
    double pad = 0.05;

    // initial polylines in normalized [pad, 1-pad] coordinates
    Mat X(m), Y(m);
    for (int e = 0; e < m; ++e) {
        double x0 = (edges_xy(e, 0) - xmin) / span * (1 - 2 * pad) + pad;
        double y0 = (edges_xy(e, 1) - ymin) / span * (1 - 2 * pad) + pad;
        double x1 = (edges_xy(e, 2) - xmin) / span * (1 - 2 * pad) + pad;
        double y1 = (edges_xy(e, 3) - ymin) / span * (1 - 2 * pad) + pad;
        X[e].resize(npoints);
        Y[e].resize(npoints);
        for (int i = 0; i < npoints; ++i) {
            double t = (double)i / (npoints - 1);
            X[e][i] = x0 + t * (x1 - x0);
            Y[e][i] = y0 + t * (y1 - y0);
        }
    }

    for (int it = 0; it < niter; ++it) {
        double h = bw * std::pow(decay, it);
        double sigma_px = h * grid;

        std::vector<double> dens(grid * grid, 0.0);
        for (int e = 0; e < m; ++e) {
            for (int i = 0; i < npoints; ++i) {
                int gx = std::min(std::max((int)std::floor(X[e][i] * (grid - 1)), 0), grid - 1);
                int gy = std::min(std::max((int)std::floor(Y[e][i] * (grid - 1)), 0), grid - 1);
                dens[gy * grid + gx] += 1.0;
            }
        }
        gaussian_blur(dens, grid, sigma_px);

        for (int e = 0; e < m; ++e) {
            for (int i = 1; i < npoints - 1; ++i) {
                double cx = X[e][i] * (grid - 1), cy = Y[e][i] * (grid - 1);
                double gxp = sample(dens, grid, cx + 1, cy) - sample(dens, grid, cx - 1, cy);
                double gyp = sample(dens, grid, cx, cy + 1) - sample(dens, grid, cx, cy - 1);
                double gn = std::sqrt(gxp * gxp + gyp * gyp);
                if (gn > 1e-12) {
                    X[e][i] += step * h * gxp / gn;
                    Y[e][i] += step * h * gyp / gn;
                }
            }
        }

        for (int e = 0; e < m; ++e) {
            for (int pass = 0; pass < smooth_passes; ++pass) {
                std::vector<double> sx = X[e], sy = Y[e];
                for (int i = 1; i < npoints - 1; ++i) {
                    X[e][i] = 0.5 * sx[i] + 0.25 * (sx[i - 1] + sx[i + 1]);
                    Y[e][i] = 0.5 * sy[i] + 0.25 * (sy[i - 1] + sy[i + 1]);
                }
            }
            resample(X[e], Y[e], npoints);
        }
    }

    List out(m);
    for (int e = 0; e < m; ++e) {
        NumericMatrix mat(npoints, 2);
        for (int i = 0; i < npoints; ++i) {
            mat(i, 0) = (X[e][i] - pad) / (1 - 2 * pad) * span + xmin;
            mat(i, 1) = (Y[e][i] - pad) / (1 - 2 * pad) * span + ymin;
        }
        out[e] = mat;
    }
    return out;
}
