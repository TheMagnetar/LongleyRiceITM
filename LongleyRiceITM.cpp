#include "LongleyRiceITM.h"

#include <cmath>
#include <iostream>

#define set_warn(txt, err)

LongleyRiceITM::LongleyRiceITM() {}

LongleyRiceITM::~LongleyRiceITM() {}

float64_t LongleyRiceITM::FORTRAN_DIM(const float64_t x, const float64_t y) {
    if (x > y) {
        return x - y;
    } else {
        return 0.0;
    }
}

float64_t LongleyRiceITM::Fn(const float64_t v_square) {
    float64_t a = 0.0;

    // [Alg 6.1]
    if (v_square <= 5.76) { // this is the 2.40 from the text, but squared
        a = 6.02 + (9.11 * std::sqrt(v_square)) - (1.27 * (v_square));
    } else {
        a = 12.953 + (4.343 * std::log(v_square));
    }

    return a;
}

float64_t LongleyRiceITM::F(const float64_t x, const float64_t K) {
    float64_t fhtv = 0.0;

    if (x <= 200.0) {
        // F = F_2(x, L), which is defined in [Alg 6.6]

        float64_t w = -std::log(K);

        // XXX the text says "or x * w^3 > 450"
        if ((K < 1e-5) || ((x * w * w * w) > 5495.0)) {
            // F_2(x, k) = F_1(x), which is defined in [Alg 6.5]
            // XXX but this isn't the same as in itm_alg.pdf
            fhtv = -117.0;
            if (x > 1.0) {
                fhtv = (17.372 * std::log(x)) + fhtv;
            }
        } else {
            // [Alg 6.6], lower part
            fhtv = (2.5e-5 * x * x / K) - (8.686 * w) - 15.0;
        }
    } else {
        // [Alg 6.3] and [Alg 6.4], lower part, which is G(x)
        fhtv = (0.05751 * x) - (4.343 * std::log(x));

        // [Alg 6.4], middle part, but again XXX this doesn't match totally
        if (x < 2000.0) {
            float64_t w = 0.0134 * x * std::exp(-0.005 * x);
            fhtv = ((1.0 - w) * fhtv) + (w * (17.372 * std::log(x) - 117.0));
        }
    }

    return fhtv;
}

float64_t LongleyRiceITM::H_0(const float64_t r, const float64_t et) {
    // constants from [Alg 6.13]
    const float64_t a[5] = { 25.0, 80.0, 177.0, 395.0, 705.0 };
    const float64_t b[5] = { 24.0, 45.0,  68.0,  80.0, 105.0 };
    int32_t it = static_cast<int32_t>(et);

    float64_t q = 0.0;
    if (it <= 0) {
        it = 1;
    } else if (it >= 4) {
        it = 4;
    } else {
        q = et - static_cast<float64_t>(it);
    }

    float64_t x = 1.0 / r;
    x *= x;
    // [Alg 6.13], calculates something like H_01(r,j), but not really XXX
    float64_t h0fv = 4.343 * std::log((1.0 + (a[it - 1] * x) + b[it - 1]) * x); // TODO: Check the indexes. There is something funny

    // XXX not sure what this means
    if (q != 0.0) {
        h0fv = ((1.0 - q) * h0fv) + (q * 4.343 * std::log((((a[it] * x) + b[it]) * x) + 1.0)); // TODO: Check the indexes. There is something funny
    }

    return h0fv;
}

float64_t LongleyRiceITM::F_0(const float64_t td) {
    // [Alg 6.9]
    if (td <= 10e3) {
        // below 10 km
        return 133.4 + (104.6 * td) + (71.8 * std::log(td));
    } else if (td <= 70e3) {
        // between 10 km and 70 km
        return 0.332e-3 + (0.212e-3 * td) + (0.157e-3 * std::log(td));
    } else {
        // above 70 km
        return -4.343 + (-1.086 * td) + (2.171 * std::log(td));
    }
}

float64_t LongleyRiceITM::adiff(const float64_t s, prop_type &prop) {
    static float64_t wd1 = 0.0, xd1 = 0.0, A_fo = 0.0, qk = 0.0, aht = 0.0, xht = 0.0;
    const float64_t A = 151.03;      // dimensionles constant from [Alg 4.20]

    if (s == 0.0) {
        std::complex<float64_t> prop_zgnd(prop.Z_g_real, prop.Z_g_imag);

        // :11: Prepare initial diffraction constants, page 5
        float64_t q = prop.h_g[0] * prop.h_g[1];
        float64_t qk = (prop.h_e[0] * prop.h_e[1]) - q;

        if (prop.mdp == ControllingMode::PointToPoint) {
            q += 10.0;                                              // "C" from [Alg 4.9]
        }

        // wd1 and xd1 are parts of Q in [Alg 4.10], but I cannot find this there
        wd1 = std::sqrt(1.0 + (qk / q));
        xd1 = prop.d_L + (prop.theta_e / prop.gamma_e);             // [Alg 4.9] upper right

        const float64_t D = 50e3;                                   // 50 km from [Alg 3.9], scale distance for \delta_h(s)
        q = (1.0 - 0.8 * std::exp(-prop.d_Ls / D)) * prop.delta_h;  // \delta_h(s), [Alg 3.9]
        const float64_t H = 16;                                     // 16 m  from [Alg 3.10]
        q *= 0.78 * std::exp(-std::pow(q / H, 0.25));               // \sigma_h(s), [Alg 3.10]

        // A_fo is the "clutter factor"
        const float64_t ALPHA = 4.77e-4;                       // from [Alg 4.10]
        A_fo = std::min(15.0, 2.171 * std::log(1.0 + ALPHA * prop.h_g[0] * prop.h_g[1] * prop.k * q)); // [Alg 4.10]

        qk = 1.0 / std::fabs(prop_zgnd);                       // qk is part of the K_j calculation from [Alg 4.17]
        aht = 20.0;                                            // 20 dB approximation for C_1(K) from [Alg 6.7], see also [Alg 4.25]
        xht = 0.0;

        for (int32_t idx = 0; idx < 2; ++idx) {
            const float64_t gamma_j_recip = 0.5 * (prop.d_Lj[idx] * prop.d_Lj[idx]) / prop.h_e[idx]; // [Alg 4.15], a is reciproke of gamma_j
            const float64_t alpha = std::pow(gamma_j_recip * prop.k, third);                   // [Alg 4.16]
            const float64_t K = qk / alpha;                                                    // [Alg 4.17]
            q = A * (1.607 - K) * alpha * prop.d_Lj[idx] / gamma_j_recip;                      // [Alg 4.18 and 6.2]
            xht += q;                                                                          // [Alg 4.19, high gain part]
            aht += F(q, K);                                                                    // [Alg 4.20] ?,    F(x, k) is in [Alg 6.4]
        }
        return 0.0;
    }

    // :12: Compute diffraction attenuation, page 5
    const float64_t theta = prop.theta_e + s * prop.gamma_e;   // [Alg 4.12]
    const float64_t ds = s - prop.d_L;                         // XXX this is not [Alg 4.13]
    float64_t q = 0.0795775 * prop.k * ds * theta * theta;

    // [Alg 4.14], float64_t knife edge attenuation. Note that the arguments to Fn() are not v, but v^2
    const float64_t A_k = Fn(q * prop.d_Lj[0] / (ds + prop.d_Lj[0])) + Fn(q * prop.d_Lj[1] / (ds + prop.d_Lj[1]));

    const float64_t gamma_o_recip = ds / theta;                 // kinda [Alg 4.15], just so that gamma_o is 1/a
    const float64_t alpha = std::pow(gamma_o_recip * prop.k, third); // [Alg 4.16]
    const float64_t K = qk / alpha;                             // [Alg 4.17], note that qk is "1.0 / abs(prop_zgnd)" from above
    q = A * (1.607 - K) * alpha * theta + xht;                  // [Alg 4.19], q is now X_o

    // looks a bit like [Alg 4.20], rounded earth attenuation, or?? note that G(x) should be "0.05751 * x - 10 * std::log(q)"
    const float64_t A_r = 0.05751 * q - 4.343 * std::log(q) - aht;

    // I'm very unsure if this has anything to do with [Alg 4.9] or not
    q = (wd1 + xd1 / s) * std::min(((1.0 - 0.8 * std::exp(-s / 50e3)) * prop.delta_h * prop.k), 6283.2);

    // XXX this is NOT the same as the weighting factor from [Alg 4.9]
    const float64_t w = 25.1 / (25.1 + std::sqrt(q));

    return ((1.0 - w) * A_k) + (w * A_r) + A_fo;                // [Alg 4.11]
}

float64_t LongleyRiceITM::A_scat(const float64_t s, prop_type &prop) {
    static float64_t ad = 0.0, rr = 0.0, etq = 0.0, h0s = 0.0;

    if (s == 0.0) {
        // :23: Prepare initial scatter constants, page 10
        ad = prop.d_Lj[0] - prop.d_Lj[1];
        rr = prop.h_e[1] / prop.h_e[0];

        if (ad < 0.0) {
            ad = -ad;
            rr = 1.0 / rr;
        }

        etq = (5.67e-6 * prop.N_s - 2.32e-3) * prop.N_s + 0.031; // part of [Alg 4.67]
        h0s = -15.0;
        return 0.0;
    }

    float64_t h0 = 0.0;

    // :24: Compute scatter attenuation, page 11
    if (h0s > 15.0) {
        h0 = h0s;
    } else {
        const float64_t theta_tick = prop.theta_ej[0] + prop.theta_ej[1] + prop.gamma_e * s; // [Alg 4.61]
        float64_t r2 = 2.0 * prop.k * theta_tick;                                            // [Alg 4.62]
        float64_t r1 = r2 * prop.h_e[0];
        r2 *= prop.h_e[1];

        if ((r1 < 0.2) && (r2 < 0.2)) {
            // The function is undefined
            return 1001.0;
        }

        // XXX not like [Alg 4.65]
        float64_t ss = (s - ad) / (s + ad);
        float64_t q = rr / ss;
        ss = std::max(0.1, ss);
        q = std::min(std::max(0.1, q), 10.0);

        const float64_t z0 = (s - ad) * (s + ad) * theta_tick * 0.25 / s;                   // XXX not like [Alg 4.66]
        float64_t temp = std::min(1.7, z0 / 8.0e3);                                         // [Alg 4.67]
        temp = temp * temp * temp * temp * temp * temp;
        const float64_t et = (etq * std::exp(-temp) + 1.0) * z0 / 1.7556e3;

        const float64_t ett = std::max(et, 1.0);
        h0 = (H_0(r1, ett) + H_0(r2, ett)) * 0.5;                                           // [Alg 6.12]
        h0 += std::min(h0, (1.38 - std::log(ett)) * std::log(ss) * std::log(q) * 0.49);     // [Alg 6.10 and 6.11]
        h0 = FORTRAN_DIM(h0, 0.0);

        if (et < 1.0) {
            // [Alg 6.14]
            h0 = et * h0 + (1.0 - et) * 4.343 * std::log(std::pow((1.0 + 1.4142 / r1) * (1.0 + 1.4142 / r2), 2.0) * (r1 + r2) / (r1 + r2 + 2.8284));
        }
        if (h0 > 15.0 && h0s >= 0.0) {
            h0 = h0s;
        }
    }

    h0s = h0;
    const float64_t theta = prop.theta_e + s * prop.gamma_e;  // [Alg 4.60]

    const float64_t D_0 = 40e3;                               // 40 km from [Alg 6.8]
    const float64_t H = 47.7;                                 // 47.7 m from [Alg 4.63]
    return 4.343 * std::log(prop.k * H * theta * theta * theta * theta) + F_0(theta * s) - 0.1 * (prop.N_s - 301.0) * std::exp(-theta * s / D_0) + h0;
}

float64_t LongleyRiceITM::abq_alos(const std::complex<float64_t> &r) {
    return std::pow(r.real(), 2) + std::pow(r.imag(), 2);
}

float64_t LongleyRiceITM::A_los(const float64_t d, prop_type &prop) {
    static float64_t wls = 0.0;

    if (d == 0.0) {
        // :18: prepare initial line-of-sight constants, page 8
        const float64_t D1 = 47.7;                                           // 47.7 m from [Alg 4.43]
        const float64_t D2 = 10e3;                                           // 10 km  from [Alg 4.43]
        const float64_t D1R = 1.0 / D1;
        // weighting factor w
        wls = D1R / (D1R + prop.k * prop.delta_h / std::max(D2, prop.d_Ls)); // [Alg 4.43]
        return 0;
    }

    std::complex<float64_t> prop_zgnd(prop.Z_g_real, prop.Z_g_imag);

    // :19: compute line of sight attentuation, page 8
    const float64_t D = 50e3;                                                // 50 km from [Alg 3.9]
    const float64_t H = 16.0;                                                // 16 m  from [Alg 3.10]
    float64_t q = (1.0 - (0.8 * std::exp(-d / D))) * prop.delta_h;           // \Delta h(d), [Alg 3.9]
    const float64_t s = 0.78 * q * std::exp(-std::pow(q / H, 0.25));         // \sigma_h(d), [Alg 3.10]
    q = prop.h_e[0] + prop.h_e[1];
    const float64_t sps = q / std::sqrt((d * d) + (q * q));    // sin(\psi), [Alg 4.46]
    std::complex<float64_t> r = (sps - prop_zgnd) / (sps + prop_zgnd) * std::exp(-std::min(10.0, prop.k * s * sps)); // [Alg 4.47]
    q = abq_alos(r);

    if ((q < 0.25) || (q < sps)) {                      // [Alg 4.48]
        r = r * std::sqrt(sps / q);
    }

    const float64_t alosv = prop.emd * d + prop.aed;    // [Alg 4.45]
    q = prop.k * prop.h_e[0] * prop.h_e[1] * 2.0 / d;   // [Alg 4.49]

    // M_PI is pi, M_PI_2 is pi/2
    if (q > M_PI_2) {                                   // [Alg 4.50]
        q = M_PI - (M_PI_2 * M_PI_2) / q;
    }

    // [Alg 4.51 and 4.44]
    return (-4.343 * std::log(abq_alos(std::complex<float64_t>(std::cos(q), -std::sin(q)) + r)) - alosv) * wls + alosv;
}

void LongleyRiceITM::lrprop(const float64_t d, prop_type &prop) {
    static bool wlos = false, wscat = false;
    static float64_t dmin = 0.0, xae = 0.0;
    std::complex<float64_t> prop_zgnd(prop.Z_g_real, prop.Z_g_imag);

    float64_t q = 0.0;

    if (prop.mdp != ControllingMode::AreaContinuation) {
        // :6: Do secondary parameters, page 3
        // [Alg 3.5]
        for (int32_t idx = 0; idx < 2; idx++) {
            prop.d_Lsj[idx] = std::sqrt(2.0 * prop.h_e[idx] / prop.gamma_e);
        }

        prop.d_Ls = prop.d_Lsj[0] + prop.d_Lsj[1];                                                // [Alg 3.6]
        prop.d_L = prop.d_Lj[0] + prop.d_Lj[1];                                                   // [Alg 3.7]
        prop.theta_e = std::max(prop.theta_ej[0] + prop.theta_ej[1], -prop.d_L * prop.gamma_e);   // [Alg 3.8]
        wlos = false;
        wscat = false;

        /* :7: Check parameters range, page 3
         * kwx is some kind of error indicator. Setting kwx to 1 means "results are slightly bad", while setting it to 4
         * means "my calculations will be bogus"
         * */

        // Frequency must be between 40 MHz and 10 GHz
        // Wave number (wn) = 2 * M_PI / lambda, or f/f0, where f0 = 47.7 MHz*m;
        // 0.838 => 40 MHz, 210 => 10 GHz
        if ((prop.k < 0.838) || (prop.k > 210.0)) {
            set_warn("Frequency not optimal", Error::NearlyOutOfRange);
            prop.kwx = std::max(prop.kwx, Error::NearlyOutOfRange);
        }

        // Surface refractivity, see [Alg 1.2]
        if ((prop.N_s < 250.0) || (prop.N_s > 400.0)) {
            set_warn("Surface refractivity out-of-bounds", Error::Other);
            prop.kwx = Error::Other;
        } else {
            // Earth's effective curvature, see [Alg 1.3]
            if ((prop.gamma_e < 75e-9) || (prop.gamma_e > 250e-9)) {
                set_warn("Earth's curvature out-of-bounds", Error::Other);
                prop.kwx = Error::Other;
            } else {
                // Surface transfer impedance, see [Alg 1.4]
                if (prop_zgnd.real() <= std::abs(prop_zgnd.imag())) {
                    set_warn("Surface transfer impedance out-of-bounds", Error::Other);
                    prop.kwx = Error::Other;
                } else {
                    // Calculating outside of 20 MHz to 40 GHz is really bad
                    if ((prop.k < 0.419) || (prop.k > 420.0)) {
                        set_warn("Frequency out-of-bounds", Error::Other);
                        prop.kwx = Error::Other;
                    } else {
                        for (int32_t idx = 0; idx < 2; idx++) {
                            // Antenna structural height should be between 1 and 1000 m
                            if ((prop.h_g[idx] < 1.0) || (prop.h_g[idx] > 1000.0)) {
                                set_warn("Antenna height not optimal", Error::NearlyOutOfRange);
                                prop.kwx = std::max(prop.kwx, Error::NearlyOutOfRange);
                            }

                            // Horizon elevation angle
                            if (std::abs(prop.theta_ej[idx]) > 200e-3) {
                                set_warn("Horizon elevation weird", Error::ParameterCombinationOutOfRange);
                                prop.kwx = std::max(prop.kwx, Error::ParameterCombinationOutOfRange);
                            }

                            // Horizon distance dl,
                            // Smooth earth horizon distance dls
                            if ((prop.d_Lj[idx] < 0.1 * prop.d_Lsj[idx]) || (prop.d_Lj[idx] > 3.0 * prop.d_Lsj[idx])) {
                                set_warn("Horizon distance weird", Error::ParameterCombinationOutOfRange);
                                prop.kwx = std::max(prop.kwx, Error::ParameterCombinationOutOfRange);
                            }
                            // Antenna structural height must between  0.5 and 3000 m
                            if ((prop.h_g[idx] < 0.5) || (prop.h_g[idx] > 3000.0)) {
                                set_warn("Antenna heights out-of-bounds", Error::Other);
                                prop.kwx = Error::Other;
                            }
                        }
                    }
                }
            }
        }

        dmin = std::abs(prop.h_e[0] - prop.h_e[1]) / 200e-3;

        /*
         * :9: Diffraction coefficients, page 4
         *
         * This is the region beyond the smoot D_Lsa and short of where isotropic scatter takes over. It is a key to central
         * region and associated coefficients  must always be computed.
         */
        q = adiff(0.0, prop);
        xae = std::pow(prop.k * prop.gamma_e * prop.gamma_e, -third);            // [Alg 4.2]
        const float64_t d3 = std::max(prop.d_Ls, 1.3787 * xae + prop.d_L);  // [Alg 4.3]
        const float64_t d4 = d3 + 2.7574 * xae;                             // [Alg 4.4]
        const float64_t a3 = adiff(d3, prop);                               // [Alg 4.5]
        const float64_t a4 = adiff(d4, prop);                               // [Alg 4.7]

        prop.emd = (a4 - a3) / (d4 - d3);                                   // [Alg 4.8]
        prop.aed = a3 - prop.emd * d3;
    }

    if (prop.mdp >= ControllingMode::AreaContinuation) {
        prop.mdp = ControllingMode::AreaContinuation;
        prop.d = d;
    }

    if (prop.d > 0.0) {
        // :8: Check distance, page 3

        // Distance above 1000 km is guessing
        if (prop.d > 1000e3) {
            set_warn("Distance not optimal", Error::NearlyOutOfRange);
            prop.kwx = std::max(prop.kwx, Error::NearlyOutOfRange);
        }

        // Distance too small, use some indoor algorithm :-)
        if (prop.d < dmin) {
            set_warn("Distance too small", Error::ParameterCombinationOutOfRange);
            prop.kwx = std::max(prop.kwx, Error::ParameterCombinationOutOfRange);
        }

        // Distance above 2000 km is really bad, don't do that
        if ((prop.d < 1e3) || (prop.d > 2000e3)) {
            set_warn("Distance out-of-bounds", Error::Other);
            prop.kwx = Error::Other;
        }
    }

    if (prop.d < prop.d_Ls) {
        // :15: Line-of-sight calculations, page 7
        if (!wlos) {
            // :16: Line-of-sight coefficients, page 7
            q = A_los(0.0, prop);

            float64_t d0 = 1.908 * prop.k * prop.h_e[0] * prop.h_e[1]; // [Alg 4.38]
            float64_t d1 = 0.0;
            const float64_t d2 = prop.d_Ls;

            const float64_t a2 = prop.aed + d2 * prop.emd;

            if (prop.aed >= 0.0) {
                d0 = std::min(d0, 0.5 * prop.d_L);                     // [Alg 4.28]
                d1 = d0 + 0.25 * (prop.d_L - d0);                      // [Alg 4.29]
            } else {
                d1 = std::max(-prop.aed / prop.emd, 0.25 * prop.d_L);  // [Alg 4.39]
            }

            const float64_t a1 = A_los(d1, prop);                      // [Alg 4.31]
            bool wq = false;

            if (d0 < d1) {
                const float64_t a0 = A_los(d0, prop);                  // [Alg 4.30]
                q = std::log(d2 / d0);
                prop.ak2 = std::max(0.0, ((d2 - d0) * (a1 - a0) - (d1 - d0) * (a2 - a0)) / ((d2 - d0) * std::log(d1 / d0) - (d1 - d0) * q)); // [Alg 4.32]
                wq = prop.aed >= 0.0 || prop.ak2 > 0.0;

                if (wq) {
                    prop.ak1 = (a2 - a0 - prop.ak2 * q) / (d2 - d0);   // [Alg 4.33]

                    if (prop.ak1 < 0.0) {
                        prop.ak1 = 0.0;                                // [Alg 4.36]
                        prop.ak2 = FORTRAN_DIM(a2, a0) / q;            // [Alg 4.35]

                        if (prop.ak2 == 0.0) {                         // [Alg 4.37]
                            prop.ak1 = prop.emd;
                        }
                    }
                }
            }

            if (!wq) {
                prop.ak1 = FORTRAN_DIM(a2, a1) / (d2 - d1);   // [Alg 4.40]
                prop.ak2 = 0.0;                               // [Alg 4.41]

                if (prop.ak1 == 0.0) {                        // [Alg 4.37]
                    prop.ak1 = prop.emd;
                }
            }

            prop.ael = a2 - prop.ak1 * d2 - prop.ak2 * std::log(d2); // [Alg 4.42]
            wlos = true;
        }

        if (prop.d > 0.0) {
            // [Alg 4.1]
            /*
             * The reference attenuation is determined as a function of the distance
             * d from 3 piecewise formulatios. This is part one.
             */
            prop.A_ref = prop.ael + prop.ak1 * prop.d + prop.ak2 * std::log(prop.d);
        }
    }

    if ((prop.d <= 0.0) || (prop.d >= prop.d_Ls)) {
        // :20: Troposcatter calculations, page 9
        if (!wscat) {
            // :21: Troposcatter coefficients
            const float64_t DS = 200e3;             // 200 km from [Alg 4.53]

            q = A_scat(0.0, prop);
            const float64_t d5 = prop.d_L + DS;     // [Alg 4.52]
            const float64_t d6 = d5 + DS;           // [Alg 4.53]
            const float64_t a6 = A_scat(d6, prop);  // [Alg 4.54]
            const float64_t a5 = A_scat(d5, prop);  // [Alg 4.55]

            if (a5 < 1000.0) {
                const float64_t HS = 47.7;          // 47.7 m from [Alg 4.59]
                prop.ems = (a6 - a5) / DS;          // [Alg 4.57]
                prop.dx = std::max(prop.d_Ls,       // [Alg 4.58]
                        std::max(prop.d_L + 0.3 * xae * std::log(HS * prop.k),
                                (a5 - prop.aed - prop.ems * d5) / (prop.emd - prop.ems)));
                prop.aes = (prop.emd - prop.ems) * prop.dx + prop.aed; // [Alg 4.59]
            } else {
                prop.ems = prop.emd;
                prop.aes = prop.aed;
                prop.dx = 10.e6;                    // [Alg 4.56]
            }
            wscat = true;
        }

        // [Alg 4.1], part two and three.
        if (prop.d > prop.dx) {
            prop.A_ref = prop.aes + prop.ems * prop.d;  // scatter region
        } else {
            prop.A_ref = prop.aed + prop.emd * prop.d;  // diffraction region
        }
    }

    prop.A_ref = std::max(prop.A_ref, 0.0);
}

void LongleyRiceITM::qlra(const SiteCriteria *const kst, const RadioClimate klimx, const VariabilityMode mdvarx, prop_type &prop, propv_type &propv) {
    float64_t q = 0.0;

    for (int32_t idx = 0; idx < 2; ++idx) {
        if (kst[idx] <= SiteCriteria::Random) {
            prop.h_e[idx] = prop.h_g[idx];  // [Alg 3.1]
        } else {
            q = 4.0;

            if (kst[idx] != SiteCriteria::Careful) {
                q = 9.0;
            }

            if (prop.h_g[idx] < 5.0) {
                q *= std::sin(0.3141593 * prop.h_g[idx]);
            }

            prop.h_e[idx] = prop.h_g[idx] + (1.0 + q) * std::exp(-std::min(20.0, 2.0 * prop.h_g[idx] / std::max(1e-3, prop.delta_h)));
        }

        // [Alg 3.3], upper function. q is d_Ls_j
        const float64_t H_3 = 5; // 5m from [Alg 3.3]
        q = std::sqrt(2.0 * prop.h_e[idx] / prop.gamma_e);
        prop.d_Lj[idx] = q * std::exp(-0.07 * std::sqrt(prop.delta_h / std::max(prop.h_e[idx], H_3)));
        prop.theta_ej[idx] = (0.65 * prop.delta_h * (q / prop.d_Lj[idx] - 1.0) - 2.0 * prop.h_e[idx]) / q; // [Alg 3.4]
    }

    prop.mdp = ControllingMode::StartOfArea;
    propv.lvar = std::max(propv.lvar, ControlSwitch::FrequencyChanged);

    if (mdvarx >= VariabilityMode::Single) {
        propv.mdvar = static_cast<int32_t>(mdvarx);
        propv.lvar = std::max(propv.lvar, ControlSwitch::MdvarChanged);
    }

    if ((klimx >= RadioClimate::Equatorial) && (klimx <= RadioClimate::MaritimeTemperateOverSea)) {
        propv.klim = klimx;
        propv.lvar = ControlSwitch::ClimateChangedOrInitialise;
    }
}

float64_t LongleyRiceITM::qerfi(const float64_t q) {
    const float64_t c0 = 2.515516698;
    const float64_t c1 = 0.802853;
    const float64_t c2 = 0.010328;
    const float64_t d1 = 1.432788;
    const float64_t d2 = 0.189269;
    const float64_t d3 = 0.001308;

    const float64_t x = 0.5 - q;
    float64_t t = std::max(0.5 - std::fabs(x), 0.000001);
    t = std::sqrt(-2.0 * std::log(t));
    float64_t v = t - ((((c2 * t + c1)) * t) + c0) / ((((((d3 * t) + d2) * t) + d1) * t) + 1.0);

    if (x < 0.0) {
        v = -v;
    }

    return v;
}

void LongleyRiceITM::qlrps(const float64_t fmhz, const float64_t zsys, const float64_t en0, const Polarization ipol, const float64_t eps, const float64_t sgm, prop_type &prop) {
    const float64_t gma = 157e-9;                                     // 157e-9 1/m    from [Alg 1.3]
    const float64_t N_1 = 179.3;                                      // 179.3 N-units from [Alg 1.3]
    const float64_t Z_0 = 376.62;                                     // 376.62 Ohm    from [Alg 1.5]

    // Frequecy -> Wave number
    prop.k = fmhz / f_0;                                              // [Alg 1.1]

    // Surface refractivity
    prop.N_s = en0;
    if (zsys != 0.0) {
        const float64_t Z_1 = 9.46e3;                                 // 9.46 km       from [Alg 1.2]
        prop.N_s *= std::exp(-zsys / Z_1);                            // [Alg 1.2]
    }

    // Earths effective curvature
    prop.gamma_e = gma * (1.0 - 0.04665 * std::exp(prop.N_s / N_1));  // [Alg 1.3]

    // Surface transfer impedance
    std::complex<float64_t> prop_zgnd(prop.Z_g_real, prop.Z_g_imag);
    std::complex<float64_t> zq(eps, Z_0 * sgm / prop.k);              // [Alg 1.5]
    prop_zgnd = std::sqrt(zq - 1.0);

    // adjust surface transfer impedance for Polarization
    if (ipol != Polarization::Horizontal) {
        prop_zgnd = prop_zgnd / zq;                                   // [Alg 1.4]
    }

    prop.Z_g_real = prop_zgnd.real();
    prop.Z_g_imag = prop_zgnd.imag();
}

float64_t LongleyRiceITM::curve(float64_t const c1, float64_t const c2, float64_t const x1, float64_t const x2, float64_t const x3, float64_t const de) {
    float64_t temp1 = (de - x2) / x3;
    float64_t temp2 = de / x1;

    temp1 *= temp1;
    temp2 *= temp2;

    return (c1 + c2 / (1.0 + temp1)) * temp2 / (1.0 + temp2);
}

float64_t LongleyRiceITM::avar(const float64_t zzt, const float64_t zzl, const float64_t zzc, prop_type &prop, propv_type &propv) {
    static int32_t kdv = 0;
    static float64_t dexa = 0.0, de = 0.0, vmd = 0.0, vs0 = 0.0, sgl = 0.0, sgtm = 0.0, sgtp = 0.0, sgtd = 0.0, tgtd = 0.0;
    static float64_t gm = 0.0, gp = 0.0, cv1 = 0.0, cv2 = 0.0, yv1 = 0.0, yv2 = 0.0, yv3 = 0.0, csm1 = 0.0, csm2 = 0.0;
    static float64_t ysm1 = 0.0, ysm2 = 0.0, ysm3 = 0.0, csp1 = 0.0, csp2 = 0.0, ysp1 = 0.0, ysp2 = 0.0, ysp3 = 0.0, csd1 = 0.0;
    static float64_t zd = 0.0, cfm1 = 0.0, cfm2 = 0.0, cfm3 = 0.0, cfp1 = 0.0, cfp2 = 0.0, cfp3 = 0.0;

    static bool no_location_variability = false, no_situation_variability = false;

    if (propv.lvar > ControlSwitch::NoRecalculationNeeded) {
        float64_t q = 0.0;
        uint32_t temp_klim = 5;
        /* :29: Climatic constants, page 15
         * Indexes are:
         *   0: equatorial
         *   1: continental suptropical
         *   2: maritime subtropical
         *   3: desert
         *   4: continental temperature
         *   5: maritime over land
         *   6: maritime over sea
         *                          equator  contsup  maritsup  desert  conttemp  mariland  marisea */
        const float64_t bv1[7]  = { -9.67,   -0.62,    1.26,    -9.21,   -0.62,   -0.39,    3.15    };
        const float64_t bv2[7]  = { 12.7,    9.19,    15.5,     9.05,    9.19,    2.86,     857.9   };
        const float64_t xv1[7]  = { 144.9e3, 228.9e3, 262.6e3,  84.1e3, 228.9e3,  141.7e3,  2222.e3 };
        const float64_t xv2[7]  = { 190.3e3, 205.2e3, 185.2e3,  101.1e3, 205.2e3, 315.9e3,  164.8e3 };
        const float64_t xv3[7]  = { 133.8e3, 143.6e3, 99.8e3,   98.6e3,  143.6e3, 167.4e3,  116.3e3 };
        const float64_t bsm1[7] = { 2.13,    2.66,    6.11,     1.98,    2.68,    6.86,     8.51    };
        const float64_t bsm2[7] = { 159.5,   7.67,    6.65,     13.11,   7.16,    10.38,    169.8   };
        const float64_t xsm1[7] = { 762.2e3, 100.4e3, 138.2e3,  139.1e3, 93.7e3,  187.8e3,  609.8e3 };
        const float64_t xsm2[7] = { 123.6e3, 172.5e3, 242.2e3,  132.7e3, 186.8e3, 169.6e3,  119.9e3 };
        const float64_t xsm3[7] = { 94.5e3,  136.4e3, 178.6e3,  193.5e3, 133.5e3, 108.9e3,  106.6e3 };
        const float64_t bsp1[7] = { 2.11,    6.87,    10.08,    3.68,    4.75,    8.58,     8.43    };
        const float64_t bsp2[7] = { 102.3,   15.53,    9.60,    159.3,   8.12,    13.97,    8.19    };
        const float64_t xsp1[7] = { 636.9e3, 138.7e3, 165.3e3,  464.4e3, 93.2e3,  216.0e3,  136.2e3 };
        const float64_t xsp2[7] = { 134.8e3, 143.7e3, 225.7e3,  93.1e3,  135.9e3, 152.0e3,  188.5e3 };
        const float64_t xsp3[7] = { 95.6e3,  98.6e3,  129.7e3,  94.2e3,  113.4e3, 122.7e3,  122.9e3 };
        const float64_t bsd1[7] = { 1.224,   0.801,   1.380,    1.000,   1.224,   1.518,    1.518   }; // bds1 -> is similar to C_D from table 5.1 at [Alg 5.8]
        const float64_t bzd1[7] = { 1.282,   2.161,   1.282,    20.0,    1.282,   1.282,    1.282   }; // bzd1 -> is similar to z_D from table 5.1 at [Alg 5.8]
        const float64_t bfm1[7] = { 1.0,     1.0,     1.0,      1.0,     0.92,    1.0,      1.0     };
        const float64_t bfm2[7] = { 0.0,     0.0,     0.0,      0.0,     0.25,    0.0,      0.0     };
        const float64_t bfm3[7] = { 0.0,     0.0,     0.0,      0.0,     1.77,    0.0,      0.0     };
        const float64_t bfp1[7] = { 1.0,     0.93,    1.0,      0.93,    0.93,    1.0,      1.0     };
        const float64_t bfp2[7] = { 0.0,     0.31,    0.0,      0.19,    0.31,    0.0,      0.0     };
        const float64_t bfp3[7] = { 0.0,     2.00,    0.0,      1.79,    2.00,    0.0,      0.0     };

        // :31: Setup variablity constants, page 16
        switch (propv.lvar) {
            default:
                // Initialization or climate change

                // if climate is wrong, use some "continental temperature" as default
                // and set error indicator
                if ((propv.klim < RadioClimate::Equatorial) || (propv.klim > RadioClimate::MaritimeTemperateOverSea)) {
                    propv.klim = RadioClimate::ContinentalTemperate;
                    prop.kwx = std::max(prop.kwx, Error::DefaultParamsSubstituted);
                    set_warn("Climate index set to continental", Error::DefaultParamsSubstituted);
                }

                // convert climate number into index into the climate tables
                temp_klim = static_cast<int32_t>(propv.klim) - 1;

                // :32: Climatic coefficients, page 17
                cv1  = bv1[temp_klim];
                cv2  = bv2[temp_klim];
                yv1  = xv1[temp_klim];
                yv2  = xv2[temp_klim];
                yv3  = xv3[temp_klim];
                csm1 = bsm1[temp_klim];
                csm2 = bsm2[temp_klim];
                ysm1 = xsm1[temp_klim];
                ysm2 = xsm2[temp_klim];
                ysm3 = xsm3[temp_klim];
                csp1 = bsp1[temp_klim];
                csp2 = bsp2[temp_klim];
                ysp1 = xsp1[temp_klim];
                ysp2 = xsp2[temp_klim];
                ysp3 = xsp3[temp_klim];
                csd1 = bsd1[temp_klim];
                zd   = bzd1[temp_klim];
                cfm1 = bfm1[temp_klim];
                cfm2 = bfm2[temp_klim];
                cfm3 = bfm3[temp_klim];
                cfp1 = bfp1[temp_klim];
                cfp2 = bfp2[temp_klim];
                cfp3 = bfp3[temp_klim];
                // fall throught

            case ControlSwitch::MdvarChanged:
                // :33: Mode of variablity coefficients, page 17

                // This code means that propv.mdvar can be
                //  0 ..  3
                // 10 .. 13, then no_location_variability is set              (no location variability)
                // 20 .. 23, then no_situation_variability is set             (no situatian variability)
                // 30 .. 33, then no_location_variability and no_situation_variability is set
                kdv = propv.mdvar;
                no_situation_variability = kdv >= 20;
                if (no_situation_variability){
                    kdv -= 20;
                }

                no_location_variability = kdv >= 10;
                if (no_location_variability) {
                    kdv -= 10;
                }

                if ((kdv < 0) || (kdv > 3)) {
                    kdv = 0;
                    set_warn("kdv set to 0", Error::DefaultParamsSubstituted);
                    prop.kwx = std::max(prop.kwx, Error::DefaultParamsSubstituted);
                }

                // fall throught

            case ControlSwitch::FrequencyChanged:
                // :34: Frequency coefficients, page 18
                q = std::log(0.133 * prop.k);
                gm = cfm1 + cfm2 / ((cfm3 * q * cfm3 * q) + 1.0);
                gp = cfp1 + cfp2 / ((cfp3 * q * cfp3 * q) + 1.0);
                // fall throught

            case ControlSwitch::AntennaHeightsChanged: {
                // :35: System coefficients, page 18
                // [Alg 5.3], effective distance

                const float64_t a_1 = 9000e3; // 9000 km from [[Alg 5.3]
                //XXX I don't have any idea how they made up the third summand,
                //XXX text says    a_1 * std::pow(k * D_1, -THIRD)
                //XXX with const float64_t D_1 = 1266; // 1266 km
                dexa = std::sqrt(2 * a_1 * prop.h_e[0]) + std::sqrt(2 * a_1 * prop.h_e[1]) + std::pow((575.7e12 / prop.k), third);
            }
            // fall throught

            case ControlSwitch::DistanceChanged: {
                // :36: Distance coefficients, page 18
                // [Alg 5.4]
                const float64_t D_0 = 130e3; // 130 km from [Alg 5.4]
                if (prop.d < dexa) {
                    de = D_0 * prop.d / dexa;
                } else {
                    de = D_0 + prop.d - dexa;
                }
            }
        }
        /*
         * Quantiles of time variability are computed using a variation of the methods described in Section 10 and Annex III.7 of
         * NBS~TN101, and also in CCIR Report {238-3}. Those methods speak of eight or nine discrete radio climates, of which seven
         * have been documented with corresponding empirical curves. It is these empirical curves to which we refer below. They are
         * all curves of quantiles of deviations versus the effective distance @de.
         */
        vmd = curve(cv1, cv2, yv1, yv2, yv3, de);             // [Alg 5.5]
        // [Alg 5.7], the slopes or "pseudo-standard deviations":
        // sgtm -> \sigma T-
        // sgtp -> \sigma T+
        sgtm = curve(csm1, csm2, ysm1, ysm2, ysm3, de) * gm;
        sgtp = curve(csp1, csp2, ysp1, ysp2, ysp3, de) * gp;
        // [Alg 5.8], ducting, "sgtd" -> \sigma TD
        sgtd = sgtp * csd1;
        tgtd = (sgtp - sgtd) * zd;

        // Location variability
        if (no_location_variability) {
            sgl = 0.0;
        } else {
            // Alg [3.9]
            q = (1.0 - 0.8 * std::exp(-prop.d / 50e3)) * prop.delta_h * prop.k;
            // [Alg 5.9]
            sgl = 10.0 * q / (q + 13.0);
        }

        // Situation variability
        if (no_situation_variability) {
            vs0 = 0.0;
        } else {
            const float64_t D = 100e3;          // 100 km
            vs0 = (5.0 + 3.0 * std::exp(-de / D));   // [Alg 5.10]
            vs0 *= vs0;
        }

        // Mark all constants as initialized
        propv.lvar = ControlSwitch::NoRecalculationNeeded;
    }


    // :37: Correct normal deviates, page 19
    float64_t zt = zzt;
    float64_t zl = zzl;
    float64_t zc = zzc;
    // kdv is derived from prop.mdvar
    switch (kdv) {
        case 0:
            // Single message mode. Time, location and situation variability are combined to form a confidence level.
            zt = zc;
            zl = zc;
            break;

        case 1:
            // Individual mode. Reliability is given by time variability. Confidence. is a combination of location
            // and situation variability.
            zl = zc;
            break;

        case 2:
            // Mobile modes. Reliability is a combination of time and location variability. Confidence. is given by the
            // situation variability.
            zl = zt;
            // case 3: Broadcast mode. like avar(zt, zl, zc).
            // Reliability is given by the two-fold statement of at least qT of the time in qL of the locations. Confidence.
            // is given by the situation variability.
    }
    if ((std::fabs(zt) > 3.1) || (std::fabs(zl) > 3.1) || (std::fabs(zc) > 3.1)) {
        set_warn("Situations variables not optimal", Error::NearlyOutOfRange);
        prop.kwx = std::max(prop.kwx, Error::NearlyOutOfRange);
    }

    // :38: Resolve standard deviations, page 19
    float64_t sgt = 0.0;
    if (zt < 0.0) {
        sgt = sgtm;
    } else if (zt <= zd) {
        sgt = sgtp;
    } else {
        sgt = sgtd + tgtd / zt;
    }
    // [Alg 5.11], situation variability
    const float64_t rt = 7.8, rl = 24.0;
    const float64_t vs = vs0 + (sgt * zt * sgt * zt) / (rt + zc * zc) + (sgl * zl * sgl * zl) / (rl + zc * zc);

    // :39: Resolve deviations, page 19
    float64_t yr = 0.0;
    if (kdv == 0) {
        propv.sgc = std::sqrt(sgt * sgt + sgl * sgl + vs);
    } else if (kdv == 1) {
        yr = sgt * zt;
        propv.sgc = std::sqrt(sgl * sgl + vs);
    } else if (kdv == 2) {
        yr = std::sqrt(sgt * sgt + sgl * sgl) * zt;
        propv.sgc = std::sqrt(vs);
    } else {
        yr = sgt * zt + sgl * zl;
        propv.sgc = std::sqrt(vs);
    }

    // [Alg 5.1], area variability
    float64_t avarv = prop.A_ref - vmd - yr - (propv.sgc * zc);

    // [Alg 5.2]
    if (avarv < 0.0) {
        avarv = avarv * (29.0 - avarv) / (29.0 - 10.0 * avarv);
    }

    return avarv;
}

void LongleyRiceITM::hzns(const float64_t *const pfl, prop_type &prop) {

    int32_t np = static_cast<int32_t>(pfl[0]);
    float64_t xi = pfl[1];
    float64_t za = pfl[2] + prop.h_g[0];
    float64_t zb = pfl[np + 2] + prop.h_g[1];
    float64_t qc = 0.5 * prop.gamma_e;
    float64_t q = qc * prop.d;
    prop.theta_ej[1] = (zb - za) / prop.d;
    prop.theta_ej[0] = prop.theta_ej[1] - q;
    prop.theta_ej[1] = -prop.theta_ej[1] - q;
    prop.d_Lj[0] = prop.d;
    prop.d_Lj[1] = prop.d;

    if (np < 2) {
        return;
    }
    float64_t sa = 0.0;
    float64_t sb = prop.d;
    bool wq = true;

    for (int32_t idx = 1; idx < np; idx++) {
        sa += xi;
        sb -= xi;
        q = pfl[idx + 2] - (qc * sa + prop.theta_ej[0]) * sa - za;

        if (q > 0.0) {
            prop.theta_ej[0] += q / sa;
            prop.d_Lj[0] = sa;
            wq = false;
        }

        if (!wq) {
            q = pfl[idx + 2] - (qc * sb + prop.theta_ej[1]) * sb - zb;

            if (q > 0.0) {
                prop.theta_ej[1] += q / sb;
                prop.d_Lj[1] = sb;
            }
        }
    }
}

void LongleyRiceITM::zlsq1(const float64_t *const z, const float64_t x1, const float64_t x2, float64_t& z0, float64_t& zn) {
    float64_t temp = FORTRAN_DIM(x1 / z[1], 0.0);
    float64_t xa = (temp >= 0.0) ? floor(temp) : ceil(temp);

    const float64_t xn = z[0];
    temp = FORTRAN_DIM(xn, x2 / z[1]);
    temp = (temp >= 0.0) ? floor(temp) : ceil(temp);
    float64_t xb = xn - temp;

    if (xb <= xa) {
        xa = FORTRAN_DIM(xa, 1.0);
        xb = xn - FORTRAN_DIM(xn, xb + 1.0);
    }

    int32_t ja = static_cast<int32_t>(xa);
    int32_t jb = static_cast<int32_t>(xb);
    const int32_t n = jb - ja;
    xa = xb - xa;
    float64_t x = -0.5 * xa;
    xb += x;
    float64_t a = 0.5 * (z[ja + 2] + z[jb + 2]);
    float64_t b = 0.5 * (z[ja + 2] - z[jb + 2]) * x;

    for (int32_t idx = 2; idx <= n; ++idx) {
        ++ja;
        x += 1.0;
        a += z[ja + 2];
        b += z[ja + 2] * x;
    }

    a /= xa;
    b = b * 12.0 / ((xa * xa + 2.0) * xa);
    z0 = a - b * xb;
    zn = a + b * (xn - xb);
}

float64_t LongleyRiceITM::qtile(const int32_t nn, float64_t *const a, const int32_t ir) {
    float64_t q = 0.0;
    int32_t j1 = 0;
    int32_t i0 = 0;
    bool done = false;
    bool goto10 = true;

    int32_t m = 0;
    int32_t n = nn;
    int32_t k = std::min(std::max(0, ir), n);

    while (!done) {
        if (goto10) {
            q = a[k];
            i0 = m;
            j1 = n;
        }

        int32_t i = i0;

        while (i <= n && a[i] >= q) {
            i++;
        }

        if (i > n) {
            i = n;
        }

        int32_t j = j1;

        while (j >= m && a[j] <= q) {
            j--;
        }

        if (j < m) {
            j = m;
        }

        if (i < j) {
            const float64_t r = a[i];
            a[i] = a[j];
            a[j] = r;
            i0 = i + 1;
            j1 = j - 1;
            goto10 = false;
        } else if (i < k) {
            a[k] = a[i];
            a[i] = q;
            m = i + 1;
            goto10 = true;
        } else if (j > k) {
            a[k] = a[j];
            a[j] = q;
            n = j - 1;
            goto10 = true;
        } else {
            done = true;
        }
    }

    return q;
}

float64_t LongleyRiceITM::dlthx(const float64_t *const pfl, const float64_t x1, const float64_t x2) {
    float64_t *s = nullptr;

    const int32_t np = static_cast<int32_t>(pfl[0]);
    float64_t xa = x1 / pfl[1];
    float64_t xb = x2 / pfl[1];
    float64_t dlthxv = 0.0;

    if (xb - xa < 2.0) {// exit out
        return dlthxv;
    }

    int32_t ka = static_cast<int32_t>(0.1 * (xb - xa + 8.0));
    ka = std::min(std::max(4, ka), 25);
    const int32_t n = (10 * ka) - 5;
    const int32_t kb = n - ka + 1;
    const float64_t sn = static_cast<float64_t>(n - 1);
    if (n + 2 < 0) {
        std::cerr << "dlthx(): Negative array size" << std::endl;
        return 0.0;
    }
    s = new float64_t[n + 2];
    if (s == nullptr) {
        std::cerr << "dlthx(): Could not allocate array of size " << n + 2 << std::endl;
        return 0.0;
    }

    s[0] = sn;
    s[1] = 1.0;
    xb = (xb - xa) / sn;
    int32_t k = static_cast<int32_t>(xa + 1.0);
    xa -= static_cast<float64_t>(k);

    for (int32_t idx = 0; idx < n; idx++) {
        // Reduce
        while ((xa > 0.0) && (k < np)) {
            xa -= 1.0;
            ++k;
        }

        s[idx + 2] = pfl[k + 2] + (pfl[k + 2] - pfl[k + 1]) * xa;
        xa = xa + xb;
    }

    zlsq1(s, 0.0, sn, xa, xb);
    xb = (xb - xa) / sn;

    for (int32_t idx = 0; idx < n; idx++) {
        s[idx + 2] -= xa;
        xa = xa + xb;
    }

    dlthxv = qtile(n - 1, s + 2, ka - 1) - qtile(n - 1, s + 2, kb - 1);
    dlthxv /= 1.0 - 0.8 * std::exp(-(x2 - x1) / 50.0e3);
    delete[] s;

    return dlthxv;
}

void LongleyRiceITM::qlrpfl(const float64_t *const pfl, const RadioClimate klimx, const int32_t mdvarx, prop_type &prop, propv_type &propv) {
    prop.d = pfl[0] * pfl[1];
    const int32_t np = static_cast<int32_t>(pfl[0]);

    // :44: determine horizons and dh from pfl, page 23
    hzns(pfl, prop);
    float64_t xl[2] = {0.0, 0.0};
    for (int32_t idx = 0; idx < 2; idx++) {
        xl[idx] = std::min(15.0 * prop.h_g[idx], 0.1 * prop.d_Lj[idx]);
    }

    xl[1] = prop.d - xl[1];
    prop.delta_h = dlthx(pfl, xl[0], xl[1]);

    float64_t q = 0.0, za = 0.0, zb = 0.0;
    if (prop.d_Lj[0] + prop.d_Lj[1] > 1.5 * prop.d) {
        // :45: Redo line-of-sight horizons, page 23

        /*
         * If the path is line-of-sight, we still need to know where the horizons might have been, and so we turn to
         * techniques used in area prediction mode.
         */
        zlsq1(pfl, xl[0], xl[1], za, zb);
        prop.h_e[0] = prop.h_g[0] + FORTRAN_DIM(pfl[2], za);
        prop.h_e[1] = prop.h_g[1] + FORTRAN_DIM(pfl[np + 2], zb);

        for (int32_t idx = 0; idx < 2; idx++) {
            prop.d_Lj[idx] = std::sqrt(2.0 * prop.h_e[idx] / prop.gamma_e) * std::exp(-0.07 * std::sqrt(prop.delta_h / std::max(prop.h_e[idx], 5.0)));
        }

        q = prop.d_Lj[0] + prop.d_Lj[1];

        if (q <= prop.d) {
            q = ((prop.d / q) * (prop.d / q));

            for (int32_t idx = 0; idx < 2; idx++) {
                prop.h_e[idx] *= q;
                prop.d_Lj[idx] = std::sqrt(2.0 * prop.h_e[idx] / prop.gamma_e) * std::exp(-0.07 * std::sqrt(prop.delta_h / std::max(prop.h_e[idx], 5.0)));
            }
        }

        for (int32_t idx = 0; idx < 2; idx++) {
            q = std::sqrt(2.0 * prop.h_e[idx] / prop.gamma_e);
            prop.theta_ej[idx] = (0.65 * prop.delta_h * (q / prop.d_Lj[idx] - 1.0) - 2.0 * prop.h_e[idx]) / q;
        }
    } else {
        // :46: Transhorizon effective heights, page 23
        zlsq1(pfl, xl[0], 0.9 * prop.d_Lj[0], za, q);
        zlsq1(pfl, prop.d - 0.9 * prop.d_Lj[1], xl[1], q, zb);
        prop.h_e[0] = prop.h_g[0] + FORTRAN_DIM(pfl[2], za);
        prop.h_e[1] = prop.h_g[1] + FORTRAN_DIM(pfl[np + 2], zb);
    }

    prop.mdp = ControllingMode::PointToPoint;
    propv.lvar = std::max(propv.lvar, ControlSwitch::FrequencyChanged);

    if (mdvarx >= 0) {
        propv.mdvar = mdvarx;
        propv.lvar = std::max(propv.lvar, ControlSwitch::MdvarChanged);
    }

    if ((klimx >= RadioClimate::Equatorial) && (klimx <= RadioClimate::MaritimeTemperateOverSea)) {
        propv.klim = klimx;
        propv.lvar = ControlSwitch::ClimateChangedOrInitialise;
    }

    lrprop(0.0, prop);
}

void LongleyRiceITM::pointToPoint(const float64_t *const elev,
        const float64_t tht_m,
        const float64_t rht_m,
        const float64_t eps_dielect,
        const float64_t sgm_conductivity,
        const float64_t eno,
        const float64_t frq_mhz,
        const RadioClimate radio_climate,
        const Polarization pol,
        const float64_t conf,
        const float64_t rel,
        float64_t &dbloss,
        PropagationMode &propMode,
        int32_t &p_mode,
        float64_t(&horizons)[2],
        Error &errnum) {

    if (elev == nullptr) {
        dbloss = -992.0;
        propMode = PropagationMode::Undefined;
        p_mode = 0;
        horizons[0] = 0.0;
        horizons[1] = 0.0;
        errnum = Error::Other;
        return;
    }

    prop_type   prop;
    prop.h_g[0] = tht_m; // Tx height above ground level
    prop.h_g[1] = rht_m; // Rx height above ground level
    prop.kwx    = Error::NoError;
    prop.mdp    = ControllingMode::PointToPoint;

    propv_type  propv;
    propv.klim = radio_climate;
    propv.lvar = ControlSwitch::ClimateChangedOrInitialise; // initialize all constants
    propv.mdvar = 12;

    float64_t zc = qerfi(conf);
    float64_t zr = qerfi(rel);
    const int64_t np = static_cast<int64_t>(elev[0]);

    const int64_t ja = static_cast<int64_t>(3.0 + 0.1 * elev[0]);
    const int64_t jb = np - ja + 6;
    float64_t zsys = 0.0;
    for (int64_t idx = ja - 1; idx < jb; ++idx) {
        zsys += elev[idx];
    }
    zsys = static_cast<float64_t>(jb - ja + 1)/zsys;

    float64_t q = eno;
    qlrps(frq_mhz, zsys, q, pol, eps_dielect, sgm_conductivity, prop);
    qlrpfl(elev, propv.klim, propv.mdvar, prop, propv);
    const float64_t fs = 32.45 + (20.0 * std::log10(frq_mhz)) + (20.0 * std::log10(prop.d / 1000.0));
    q = prop.d - prop.d_L;

    horizons[0] = 0.0;
    horizons[1] = 0.0;
    if (static_cast<int32_t>(q) < 0) {
        propMode = PropagationMode::LineOfSight;
        p_mode = 0;
    } else {
        if (static_cast<int32_t>(q) == 0) {
            propMode = PropagationMode::SingleHorizon;
            horizons[0] = prop.d_Lj[0];
            p_mode = 1;
        } else if (static_cast<int32_t>(q) > 0) {
            propMode = PropagationMode::DoubleHorizon;

            horizons[0] = prop.d_Lj[0];
            horizons[1] = prop.d_Lj[1];
            p_mode = 1;
        }

        if ((prop.d <= prop.d_Ls) || (prop.d <= prop.dx)) {
            propMode = static_cast<PropagationMode>(static_cast<int32_t>(propMode) + 1); // Diffraction Dominant
            p_mode = 1;
        } else if (prop.d > prop.dx) {
            propMode = static_cast<PropagationMode>(static_cast<int32_t>(propMode) + 2); // Troposcatter Dominant
            p_mode = 2;
        }
    }

    dbloss = avar(zr, 0.0, zc, prop, propv) + fs;
    errnum = prop.kwx;
}

void LongleyRiceITM::pointToPointMDH(const float64_t *const elev,
        const float64_t tht_m,
        const float64_t rht_m,
        const float64_t eps_dielect,
        const float64_t sgm_conductivity,
        const float64_t eno,
        const float64_t frq_mhz,
        const RadioClimate radio_climate,
        const Polarization pol,
        const float64_t timepct,
        const float64_t locpct,
        const float64_t confpct,
        float64_t &dbloss,
        PropagationMode &propMode,
        float64_t &deltaH,
        Error &errnum) {

    if (elev == nullptr) {
        dbloss = -992.0;
        propMode = PropagationMode::Undefined;
        deltaH = 0.0;
        errnum = Error::Other;
        return;
    }

    prop_type   prop;
    prop.h_g[0] = tht_m;
    prop.h_g[1] = rht_m;
    prop.kwx    = Error::NoError;
    prop.mdp = ControllingMode::PointToPoint;

    propv_type  propv;
    propv.lvar = ControlSwitch::ClimateChangedOrInitialise;
    propv.klim = radio_climate;
    propv.mdvar = 12;

    const float64_t ztime = qerfi(timepct);
    const float64_t zloc  = qerfi(locpct);
    const float64_t zconf = qerfi(confpct);

    const int64_t np = static_cast<int64_t>(elev[0]);

    const int64_t ja = static_cast<int64_t>(3.0 + 0.1 * elev[0]);
    const int64_t jb = np - ja + 6;
    float64_t zsys = 0.0;
    for (int64_t i = ja - 1; i < jb; ++i) {
        zsys += elev[i];
    }
    zsys = static_cast<float64_t>(jb - ja + 1) / zsys;

    float64_t q = eno;
    qlrps(frq_mhz, zsys, q, pol, eps_dielect, sgm_conductivity, prop);
    qlrpfl(elev, propv.klim, propv.mdvar, prop, propv);
    const float64_t fs = 32.45 + 20.0 * std::log10(frq_mhz) + 20.0 * std::log10(prop.d / 1000.0);
    deltaH = prop.delta_h;
    q = prop.d - prop.d_L;

    propMode = PropagationMode::Undefined;
    if (static_cast<int32_t>(q) < 0) {
        propMode = PropagationMode::LineOfSight;
    } else {
        if (static_cast<int32_t>(q) == 0) {
            propMode = PropagationMode::SingleHorizon;
        } else if (static_cast<int32_t>(q) > 0) {
            propMode = PropagationMode::DoubleHorizon;
        }

        if ((prop.d <= prop.d_Ls) || (prop.d <= prop.dx)) {
            propMode = static_cast<PropagationMode>(static_cast<int32_t>(propMode) + 1); // Diffraction Dominant
        } else if (prop.d > prop.dx) {
            propMode = static_cast<PropagationMode>(static_cast<int32_t>(propMode) + 2); // Troposcatter Dominant
        }
    }

    dbloss = avar(ztime, zloc, zconf, prop, propv) + fs;  //avar(time,location,confidence)
    errnum = prop.kwx;
}

void LongleyRiceITM::pointToPointDH(const float64_t *const elev,
        const float64_t tht_m,
        const float64_t rht_m,
        const float64_t eps_dielect,
        const float64_t sgm_conductivity,
        const float64_t eno,
        const float64_t frq_mhz,
        const RadioClimate radio_climate,
        const Polarization pol,
        const float64_t conf,
        const float64_t rel,
        float64_t &dbloss,
        PropagationMode &propMode,
        float64_t &deltaH,
        Error &errnum) {

    if (elev == nullptr) {
        dbloss = -992.0;
        propMode = PropagationMode::Undefined;
        deltaH = 0.0;
        errnum = Error::Other;
        return;
    }

    prop_type prop;
    prop.h_g[0] = tht_m;
    prop.h_g[1] = rht_m;
    prop.kwx    = Error::NoError;
    prop.mdp = ControllingMode::PointToPoint;
    prop.kwx = Error::NoError;

    propv_type propv;
    propv.klim = radio_climate;
    propv.lvar = ControlSwitch::ClimateChangedOrInitialise;
    propv.mdvar = 12;

    const float64_t zc = qerfi(conf);
    const float64_t zr = qerfi(rel);
    const int64_t np = static_cast<int64_t>(elev[0]);

    const int64_t ja = static_cast<int64_t>(3.0 + 0.1 * elev[0]);
    const int64_t jb = np - ja + 6;
    float64_t zsys = 0;
    for (int64_t idx = ja - 1; idx < jb; ++idx) {
        zsys += elev[idx];
    }
    zsys = static_cast<float64_t>(jb - ja + 1)/zsys;

    float64_t q = eno;
    qlrps(frq_mhz, zsys, q, pol, eps_dielect, sgm_conductivity, prop);
    qlrpfl(elev, propv.klim, propv.mdvar, prop, propv);
    const float64_t fs = 32.45 + 20.0 * std::log10(frq_mhz) + 20.0 * std::log10(prop.d / 1000.0);
    deltaH = prop.delta_h;
    q = prop.d - prop.d_L;

    propMode = PropagationMode::Undefined;
    if (static_cast<int32_t>(q) < 0) {
        propMode = PropagationMode::LineOfSight;
    } else {
        if (static_cast<int32_t>(q) == 0) {
            propMode = PropagationMode::SingleHorizon;
        } else if (static_cast<int32_t>(q) > 0) {
            propMode = PropagationMode::DoubleHorizon;
        }

        if ((prop.d <= prop.d_Ls) || (prop.d <= prop.dx)) {
            propMode = static_cast<PropagationMode>(static_cast<int32_t>(propMode) + 1); // Diffraction Dominant
        } else if (prop.d > prop.dx) {
            propMode = static_cast<PropagationMode>(static_cast<int32_t>(propMode) + 2); // Troposcatter Dominant
        }
    }

    dbloss = avar(zr, 0.0, zc, prop, propv) + fs; //avar(time,location,confidence)
    errnum = prop.kwx;
}

void LongleyRiceITM::area(const VariabilityMode ModVar,
        const float64_t deltaH,
        const float64_t tht_m,
        const float64_t rht_m,
        const float64_t dist_km,
        const SiteCriteria TSiteCriteria,
        const SiteCriteria RSiteCriteria,
        const float64_t eps_dielect,
        const float64_t sgm_conductivity,
        const float64_t eno,
        const float64_t frq_mhz,
        const RadioClimate radio_climate,
        const Polarization pol,
        const float64_t pctTime,
        const float64_t pctLoc,
        const float64_t pctConf,
        float64_t &dbloss,
        Error &errnum) {

    SiteCriteria kst[2] = { TSiteCriteria, RSiteCriteria };

    const float64_t zt = qerfi(pctTime);
    const float64_t zl = qerfi(pctLoc);
    const float64_t zc = qerfi(pctConf);

    prop_type prop;
    prop.delta_h = deltaH;
    prop.h_g[0] = tht_m;
    prop.h_g[1] = rht_m;
    prop.N_s = eno;
    prop.kwx = Error::NoError;

    propv_type propv;
    propv.klim = radio_climate;

    qlrps(frq_mhz, 0.0, eno, pol, eps_dielect, sgm_conductivity, prop);
    qlra(kst, propv.klim, ModVar, prop, propv);

    if (propv.lvar < ControlSwitch::DistanceChanged) {
        propv.lvar = ControlSwitch::DistanceChanged;
    }

    lrprop(dist_km * 1000.0, prop);
    const float64_t fs = 32.45 + (20.0 * std::log10(frq_mhz)) + (20.0 * std::log10(prop.d / 1000.0));
    dbloss = fs + avar(zt, zl, zc, prop, propv);

    errnum = prop.kwx;
}
