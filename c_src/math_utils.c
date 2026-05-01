/**
 * math_utils.c
 * =============
 * C implementations of numerically critical geometric operations for H2 (Poincaré disk).
 * These functions are called from Ada SPARK via bindings in math_utils.ads.
 *
 * Compile with: gcc -c -O2 -g math_utils.c -o math_utils.o
 */

#include <math.h>
#include <stdio.h>

// --- Type Definitions ---
typedef struct {
    double x;
    double y;
} Point2D;

// --- Hyperbolic Geometry (H2) ---

/**
 * Möbius addition for Poincaré disk (Eq. 3 from the paper).
 * z ⊕ w = [(1 + 2<z,w> + ||w||²)z + (1 - ||z||²)w] / [1 + 2<z,w> + ||z||²||w||²]
 */
Point2D mobius_add(Point2D z, Point2D w) {
    double zx = z.x, zy = z.y;
    double wx = w.x, wy = w.y;

    double z_norm_sq = zx*zx + zy*zy;
    double w_norm_sq = wx*wx + wy*wy;
    double dot = zx*wx + zy*wy;

    double denom = 1.0 + 2.0*dot + z_norm_sq * w_norm_sq;
    if (fabs(denom) < 1e-12) {
        // Avoid division by zero (should not happen for valid inputs)
        fprintf(stderr, "Warning: Division by zero in mobius_add (denom=%.2e)\n", denom);
        denom = 1e-12;
    }

    Point2D result;
    result.x = ((1.0 + 2.0*dot + w_norm_sq) * zx + (1.0 - z_norm_sq) * wx) / denom;
    result.y = ((1.0 + 2.0*dot + w_norm_sq) * zy + (1.0 - z_norm_sq) * wy) / denom;

    return result;
}

/**
 * Hyperbolic distance in Poincaré disk (Eq. 3 from the paper).
 * d_D(z, w) = 2 * arctanh(||(-z) ⊕ w||)
 */
double hyperbolic_distance(Point2D z, Point2D w) {
    Point2D diff = mobius_add((Point2D){-z.x, -z.y}, w);
    double norm = sqrt(diff.x*diff.x + diff.y*diff.y);

    // Clamp norm to [0, 1) to avoid NaN in arctanh
    if (norm >= 1.0) {
        norm = 1.0 - 1e-12;
    } else if (norm < 0.0) {
        norm = 0.0;
    }

    return 2.0 * atanh(norm);
}

/**
 * Unit hyperbolic tangent direction (for angle-based insertion, Eq. 8).
 * Returns the unit vector in the direction from z to w in H2.
 */
Point2D hyperbolic_tangent(Point2D z, Point2D w) {
    Point2D diff = mobius_add((Point2D){-z.x, -z.y}, w);
    double norm = sqrt(diff.x*diff.x + diff.y*diff.y);

    if (norm < 1e-12) {
        // Return zero vector if points are identical
        return (Point2D){0.0, 0.0};
    }

    return (Point2D){diff.x / norm, diff.y / norm};
}

/**
 * Rotation of a 2D vector by angle alpha (in radians).
 */
Point2D rotate_vector(Point2D v, double alpha) {
    double cos_alpha = cos(alpha);
    double sin_alpha = sin(alpha);
    return (Point2D){
        v.x * cos_alpha - v.y * sin_alpha,
        v.x * sin_alpha + v.y * cos_alpha
    };
}

/**
 * Clamp a point to the Poincaré disk (||z|| < 0.999).
 */
Point2D clamp_to_disk(Point2D z) {
    double norm_sq = z.x*z.x + z.y*z.y;
    if (norm_sq >= 0.999*0.999) {
        double norm = sqrt(norm_sq);
        double scale = 0.999 / norm;
        return (Point2D){z.x * scale, z.y * scale};
    }
    return z;
}

// --- Euclidean Geometry (E2) ---

/**
 * Euclidean distance between two points.
 */
double euclidean_distance(Point2D a, Point2D b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx*dx + dy*dy);
}

/**
 * Normalize a 2D vector.
 */
Point2D normalize_vector(Point2D v) {
    double norm = sqrt(v.x*v.x + v.y*v.y);
    if (norm < 1e-12) {
        return (Point2D){0.0, 0.0};
    }
    return (Point2D){v.x / norm, v.y / norm};
}

// --- Spherical Geometry (S2) ---

/**
 * Spherical distance (great-circle distance) between two unit vectors.
 */
double spherical_distance(double ax, double ay, double az,
                           double bx, double by, double bz) {
    double dot = ax*bx + ay*by + az*bz;
    // Clamp dot product to [-1, 1] to avoid NaN in acos
    if (dot > 1.0) dot = 1.0;
    if (dot < -1.0) dot = -1.0;
    return acos(dot);
}

/**
 * Spherical tangent vector from point p to q (Eq. 2 from the paper).
 * Returns the unit tangent vector in the direction from p to q on S2.
 */
void spherical_tangent(double px, double py, double pz,
                       double qx, double qy, double qz,
                       double *tx, double *ty, double *tz) {
    // v = q - <p, q> * p
    double dot = px*qx + py*qy + pz*qz;
    double vx = qx - dot * px;
    double vy = qy - dot * py;
    double vz = qz - dot * pz;

    // Normalize v
    double norm = sqrt(vx*vx + vy*vy + vz*vz);
    if (norm < 1e-12) {
        *tx = *ty = *tz = 0.0;
        return;
    }

    *tx = vx / norm;
    *ty = vy / norm;
    *tz = vz / norm;
}
