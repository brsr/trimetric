#define PJ_LIB__

#include <errno.h>
#include <math.h>
#include <stdio.h>

#include "proj.h"
#include "proj_internal.h"

#define EPS4 1.e-4
#define TOL EPS4
#define MAXITER 10

namespace { // anonymous namespace
struct pj_opaque {
    struct {             // control point data
        double phi, lam; // sphere coords of control point
        double cosphi, sinphi;
        double s; // side length
        PJ_XY xy; // projected coords of control point
    } c[3];
    // matrices
    double Mmat[3][2];
    double Mpinv[3][2];
    double Vinv[3][3];
    double Amat[3][3];
};
} // anonymous namespace

PROJ_HEAD(mattri, "Matrix Trimetric")
"\n\tMisc Sph"
    "\n\tlat_1= lon_1= lat_2= lon_2= lat_3= lon_3=";

/* (spherical approximation to) distance from point 1 to point 2 */
static double dist(PJ_CONTEXT *ctx, double dphi, double c1, double s1,
                   double c2, double s2, double dlam) {
    double r, cdl, dp, dl;

    cdl = cos(dlam);
    if (fabs(dphi) > 1. || fabs(dlam) > 1.)
        r = aacos(ctx, s1 * s2 + c1 * c2 * cdl);
    else { /* more accurate for smaller distances */
        dp = sin(.5 * dphi);
        dl = sin(.5 * dlam);
        r = 2. * aasin(ctx, sqrt(dp * dp + c1 * c2 * dl * dl));
    }
    return r;
}

/* law of cosines */
static double lc(PJ_CONTEXT *ctx, double b, double c, double a) {
    return aacos(ctx, .5 * (b * b + c * c - a * a) / (b * c));
}

static PJ_XY mattri_s_forward(PJ_LP lp, PJ *P) { /* Spherical, forward */
    PJ_XY xy;
    struct pj_opaque *Q = static_cast<struct pj_opaque *>(P->opaque);
    double rsq;
    int i;
    double sinphi = sin(lp.phi);
    double cosphi = cos(lp.phi);

    xy.x = 0;
    xy.y = 0;

    for (i = 0; i < 3; ++i) { /* dist from control pts */
        rsq = dist(P->ctx, lp.phi - Q->c[i].phi, Q->c[i].cosphi, Q->c[i].sinphi,
                   cosphi, sinphi, lp.lam - Q->c[i].lam);
        rsq *= rsq;
        xy.x += Q->Mmat[i][0] * rsq;
        xy.y += Q->Mmat[i][1] * rsq;
    }
    return xy;
}

static PJ_LP mattri_s_inverse(PJ_XY xy, PJ *P) {
    PJ_LP lp;
    double k[3];
    double r[3];
    double cosr[3];
    double sincr[3];
    double Ac;
    double v[3] = {0.0, 0.0, 0.0};
    double f, fprime;
    int i, j, m;
    struct pj_opaque *Q = static_cast<struct pj_opaque *>(P->opaque);
    double h = 0.0;

    // sum k[i] = 0, can use that to skip part of the matrix product
    k[2] = 0.0;
    for (i = 0; i < 2; ++i) {
        k[i] = Q->Mpinv[i][0] * xy.x + Q->Mpinv[i][1] * xy.y;
        k[2] -= k[i];
        if (-k[i] > h) {
            h = -k[i];
        }
    }
    if (-k[2] > h) {
        h = -k[2];
    }

    for (m = 0; m < MAXITER; ++m) { // newton's method iteration
        for (i = 0; i < 3; ++i) {
            r[i] = k[i] + h;
            if (r[i] > EPS4) {
                r[i] = sqrt(r[i]);
                cosr[i] = cos(r[i]);
                sincr[i] = sin(r[i]) / r[i];
            } else {
                /*linearization of cos(sqrt(x)) and sinc(sqrt(x)) near 0.
                this is exactly equal for doubles where abs(r[i])) < EPS4*/
                cosr[i] = 1 - r[i] / 2.0;
                sincr[i] = 1 - r[i] / 6.0;
            }
        }
        f = -1.0;
        fprime = 0.0;
        for (i = 0; i < 3; ++i) {
            Ac = 0.0;
            for (j = 0; j < 3; ++j) {
                Ac += Q->Amat[i][j] * cosr[j];
            }
            f += Ac * cosr[i];
            fprime += Ac * sincr[i];
        }
        if (fprime <= 0) {
            break; // only happens when control triangle is almost a hemisphere
        }
        h += f / fprime;
        if (abs(f) < TOL) {
            break;
        }
    }

    for (i = 0; i < 3; ++i) { // FIXME DRY
        r[i] = k[i] + h;
        if (r[i] < 0) {
            // if it's still negative by now, assume it's really 0
            r[i] = 0;
        } else {
            r[i] = sqrt(r[i]);
        }
        cosr[i] = cos(r[i]);
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            v[i] += Q->Vinv[i][j] * cosr[j];
        }
    }
    lp.phi = atan(v[2] / sqrt(v[0] * v[0] + v[1] * v[1]));
    lp.lam = atan2(v[1], v[0]);
    return lp;
}

PJ *PROJECTION(mattri) {

    int i, j, k;
    double phi0, phi1, temp;
    double V[3][3];
    char line[10];
    struct pj_opaque *Q =
        static_cast<struct pj_opaque *>(calloc(1, sizeof(struct pj_opaque)));
    if (nullptr == Q)
        return pj_default_destructor(P, PROJ_ERR_OTHER /*ENOMEM*/);
    P->opaque = Q;
    for (i = 0; i < 3; ++i) { /* get control point locations */
        (void)sprintf(line, "rlat_%d", i + 1);
        Q->c[i].phi = pj_param(P->ctx, P->params, line).f;
        (void)sprintf(line, "rlon_%d", i + 1);
        Q->c[i].lam = pj_param(P->ctx, P->params, line).f;
        Q->c[i].lam = adjlon(Q->c[i].lam - P->lam0);
        Q->c[i].cosphi = cos(Q->c[i].phi);
        Q->c[i].sinphi = sin(Q->c[i].phi);
    }
    for (i = 0; i < 3; ++i) { /* inter-ctl pt. distances */
        j = (i + 1) % 3;
        k = (i + 2) % 3;
        Q->c[i].s = dist(P->ctx, Q->c[j].phi - Q->c[k].phi, Q->c[k].cosphi,
                         Q->c[k].sinphi, Q->c[j].cosphi, Q->c[j].sinphi,
                         Q->c[j].lam - Q->c[k].lam);
        if (Q->c[i].s == 0.0) {
            proj_log_error(
                P,
                _("Invalid value for control points: they should be distinct"));
            return pj_default_destructor(P,
                                         PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }
    }
    Q->c[2].xy.y = 0.0;
    Q->c[2].xy.x = Q->c[0].s * Q->c[1].s * Q->c[2].s; // circumradius
    Q->c[2].xy.x /= sqrt(2 * pow(Q->c[0].s * Q->c[1].s, 2) +
                         2 * pow(Q->c[1].s * Q->c[2].s, 2) +
                         2 * pow(Q->c[2].s * Q->c[0].s, 2) - pow(Q->c[0].s, 4) -
                         pow(Q->c[1].s, 4) - pow(Q->c[2].s, 4));
    // rest of the control points in the plane
    phi0 = lc(P->ctx, Q->c[1].s, Q->c[2].s, Q->c[0].s);
    phi1 = lc(P->ctx, Q->c[2].s, Q->c[0].s, Q->c[1].s);
    Q->c[0].xy.x = Q->c[2].xy.x * cos(2 * phi1);
    Q->c[1].xy.x = Q->c[2].xy.x * cos(2 * phi0);
    Q->c[0].xy.y = Q->c[2].xy.x * sin(2 * phi1);
    Q->c[1].xy.y = -Q->c[2].xy.x * sin(2 * phi0);
    /* matrices. some of these are a different transpose
    than in the text for efficiency. some get calculated in row-major but
    are used as column-major later, and this part only gets calculated once */
    temp = (Q->c[0].xy.x * Q->c[1].xy.y + Q->c[1].xy.x * Q->c[2].xy.y +
            Q->c[2].xy.x * Q->c[0].xy.y - Q->c[0].xy.x * Q->c[2].xy.y -
            Q->c[1].xy.x * Q->c[0].xy.y - Q->c[2].xy.x * Q->c[1].xy.y);
    if (temp < 0.) {
        proj_log_error(
            P, _("Invalid value for control points: wrong orientation"));
        return pj_default_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
    } else if (temp == 0.) {
        proj_log_error(P, _("Invalid value for control points: collinear"));
        return pj_default_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
    }
    for (i = 0; i < 3; ++i) {
        j = (i + 1) % 3;
        k = (i + 2) % 3;
        Q->Mmat[i][0] = (Q->c[k].xy.y - Q->c[j].xy.y) / 2 / temp;
        Q->Mmat[i][1] = (Q->c[j].xy.x - Q->c[k].xy.x) / 2 / temp;
        Q->Mpinv[i][0] =
            2.0 / 3.0 * (-2.0 * Q->c[i].xy.x + Q->c[j].xy.x + Q->c[k].xy.x);
        Q->Mpinv[i][1] =
            2.0 / 3.0 * (-2.0 * Q->c[i].xy.y + Q->c[j].xy.y + Q->c[k].xy.y);
        V[0][i] = Q->c[i].cosphi * cos(Q->c[i].lam);
        V[1][i] = Q->c[i].cosphi * sin(Q->c[i].lam);
        V[2][i] = Q->c[i].sinphi;
    }
    /* matrix inverse: not the most numerically sound way to do this but
    importing LAPACK would be overkill, and this is probably OK
    for this use case*/
    for (i = 0; i < 3; ++i) {
        j = (i + 1) % 3;
        k = (i + 2) % 3;
        Q->Vinv[i][0] = V[j][1] * V[k][2] - V[j][2] * V[k][1];
        Q->Vinv[i][1] = V[j][2] * V[k][0] - V[j][0] * V[k][2];
        Q->Vinv[i][2] = V[j][0] * V[k][1] - V[j][1] * V[k][0];
    }
    temp = Q->Vinv[0][0] * V[0][0] + Q->Vinv[0][1] * V[0][1] +
           Q->Vinv[0][2] * V[0][2];
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            Q->Vinv[i][j] /= temp;
        }
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            Q->Amat[i][j] = 0;
            for (k = 0; k < 3; ++k) {
                Q->Amat[i][j] += Q->Vinv[k][i] * Q->Vinv[k][j];
            }
        }
    }
    P->es = 0.0;
    P->fwd = mattri_s_forward;
    P->inv = mattri_s_inverse;

    return P;
}
