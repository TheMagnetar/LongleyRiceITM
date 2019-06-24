/*****************************************************************************\
 *                                                                           *
 *  The following code was derived from public domain ITM code available     *
 *  at ftp://flattop.its.bldrdoc.gov/itm/ITMDLL.cpp that was released on     *
 *  June 26, 2007.  It was modified to remove Microsoft Windows "dll-isms",  *
 *  redundant and unnecessary #includes, redundant and unnecessary { }'s,    *
 *  to initialize uninitialized variables, type cast some variables,         *
 *  re-format the code for easier reading, and to replace pow() function     *
 *  calls with explicit multiplications wherever possible to increase        *
 *  execution speed and improve computational accuracy.                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Added comments that refer to itm.pdf and itmalg.pdf in a way to easy     *
 *  searching.                                                               *
 *                                                                           *
 * // :0: Blah, page 0     This is the implementation of the code from       *
 *                         itm.pdf, section "0". The description is          *
 *                         found on page 0.                                  *
 * [Alg 0.0]               please refer to algorithm 0.0 from itm_alg.pdf    *
 *                                                                           *
 * Holger Schurig, DH3HS                                                     *
 *****************************************************************************
 *                                                                           *
 *  More C++ improvements.                                                   *
 *                                                                           *
 *  Additional changes include using class enums for code clarity, use of    *
 *  const correctness, move to its own class, initialise all variables upon  *
 *  declaration, declare variables as local as possible, remove if/else      *
 *  FORTRAN style and use of fixed types for increased portability.          *
 *                                                                           *
 *  Documentation ported to Doxygen.                                         *
 *                                                                           *
 * Ferran Obón Santacana (Magnetar)                                          *
\*****************************************************************************/


// *************************************
// C++ routines for this program are taken from
// a translation of the FORTRAN code written by
// U.S. Department of Commerce NTIA/ITS
// Institute for Telecommunication Sciences
// *****************
// Irregular Terrain Model (ITM) (Longley-Rice)
// *************************************

#ifndef LONGLEYRICEITM_H_
#define LONGLEYRICEITM_H_

#include <cstdint>
#include <complex>
#include <string>

#include "Types.h"

class LongleyRiceITM {
public:
    LongleyRiceITM();
    virtual ~LongleyRiceITM();

    /**
     * Radio Climante
     */
    enum class RadioClimate : uint8_t {
        Equatorial = 1,             /*!< Equatorial */
                ContinentalSubtropical,     /*!< Continental Subtropical */
                MaritimeTropical,           /*!< Maritime Tropical */
                Desert,                     /*!< Desert */
                ContinentalTemperate,       /*!< Continental Temperate */
                MaritimeTemperateOverLand,  /*!< Maritime Temperate, Over Land */
                MaritimeTemperateOverSea,   /*!< Maritime Temperate, Over Sea */
    };

    /**
     * Error code of the Longley-Rice routines
     */
    enum class Error : uint8_t {
        NoError,                        /*!< No Error */
        NearlyOutOfRange,               /*!< Some parameters are nearly out of range. Results should be used with caution. */
        DefaultParamsSubstituted,       /*!< Default parameters have been substituted for impossible ones. */
        ParameterCombinationOutOfRange, /*!< Warning: A combination of parameters is out of range. Results are probably invalid. */
        Other                           /*!< Some parameters are out of range. Results are probably invalid. */
    };

    /**
     * Antenna polarization
     */
    enum class Polarization : uint8_t {
        Horizontal,
        Vertical
    };

    /**
     * Propagation mode
     */
    enum class PropagationMode : int8_t {
        Undefined  = static_cast<int8_t>(-1),
                LineOfSight = 0,
                SingleHorizon = 4,
                SingleHorizon_Diffraction = 5,
                SingleHorizon_Troposcatter = 6,
                DoubleHorizon = 8,
                DoubleHorizon_Diffraction = 9,
                DoubleHorizon_Troposcatter = 10
    };

    /**
     * Site criteria for area calculation mode
     */
    enum class SiteCriteria : uint8_t {
        Random,        /*!< Random */
        Careful,       /*!< Careful */
        VeryCareful    /*!< Very Careful */
    };

    enum class VariabilityMode : uint8_t {
        Single,
        Individual,
        Mobile,
        Broadcast
    };

    /**
     * Point to point calculation
     *
     * @params[in]   elev              [num points - 1], [delta dist(meters)], [height(meters) point 1], ..., [height(meters) point n]
     * @params[in]   tht_m             Transceiver above ground level.
     * @params[in]   rht_m             Receiver above ground level.
     * @params[in]   eps_dielect       Earth dielectric constant (rel. permittivity).
     * @params[in]   sgm_conductivity  Earth conductivity (Siemens/m).
     * @params[in]   eno               Atmospheric bending const, n-Units
     * @params[in]   frq_mhz           Frequency in MHz.
     * @params[in]   radio_climate     Radio Climate
     * @params[in]   pol               Polarization
     * @params[in]   conf              0.01 .. .99, Fractions of situations
     * @params[in]   rel               0.01 .. .99, Fractions of time
     * @params[out]  dblos             Loss in DB
     * @params[out]  propMode          Propagation Mode
     * @params[out]  p_mode            propagation mode selector
     * @params[out]  horizons          Horizon distances
     * @params[out]  errnum            Error code
     */
    void pointToPoint(const float64_t *const elev,
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
            Error &errnum);

    /**
     * Point to point calculation with MDH
     *
     * @params[in]   elev              [num points - 1], [delta dist(meters)], [height(meters) point 1], ..., [height(meters) point n]
     * @params[in]   tht_m             Transceiver above ground level.
     * @params[in]   rht_m             Receiver above ground level.
     * @params[in]   eps_dielect       Earth dielectric constant (rel. permittivity).
     * @params[in]   sgm_conductivity  Earth conductivity (Siemens/m).
     * @params[in]   eno               Atmospheric bending const, n-Units
     * @params[in]   frq_mhz           Frequency in MHz.
     * @params[in]   radio_climate     Radio Climate
     * @params[in]   pol               Polarization
     * @params[in]   timepct           0.01 .. .99
     * @params[in]   locpct            0.01 .. .99
     * @params[in]   confpct           0.01 .. .99
     * @params[out]  dblos             Loss in DB
     * @params[out]  propMode          Propagation Mode
     * @params[out]  deltaH
     * @params[out]  errnum            Error code
     */
    void pointToPointMDH(const float64_t *const elev,
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
            Error &errnum);

    /**
     * Point to point calculation with MDH
     *
     * @params[in]   elev              [num points - 1], [delta dist(meters)], [height(meters) point 1], ..., [height(meters) point n]
     * @params[in]   tht_m             Transceiver above ground level.
     * @params[in]   rht_m             Receiver above ground level.
     * @params[in]   eps_dielect       Earth dielectric constant (rel. permittivity).
     * @params[in]   sgm_conductivity  Earth conductivity (Siemens/m).
     * @params[in]   eno               Atmospheric bending const, n-Units
     * @params[in]   frq_mhz           Frequency in MHz.
     * @params[in]   radio_climate     Radio Climate
     * @params[in]   pol               Polarization
     * @params[in]   conf              0.01 .. .99
     * @params[in]   rel               0.01 .. .99
     * @params[out]  dbloss            propagation mode selector
     * @params[out]  deltaH
     * @params[out]  errnum            Error code
     */
    void pointToPointDH(const float64_t *const elev,
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
            Error &errnum);

    /**
     * Area calculation
     *
     * @params[in]   ModVar            0 - Single: pctConf is "Time/Situation/Location", pctTime, pctLoc not used
     *                                 1 - Individual: pctTime is "Situation/Location", pctConf is "Confidence", pctLoc not used
     *                                 2 - Mobile: pctTime is "Time/Locations (Reliability)", pctConf is "Confidence", pctLoc not used
     *                                 3 - Broadcast: pctTime is "Time", pctLoc is "Location", pctConf is "Confidence"
     * @params[in]   deltaH
     * @params[in]   tht_m             Transceiver above ground level.
     * @params[in]   rht_m             Receiver above ground level.
     * @params[in]   dist_km           Distance between transmitter and receiver
     * @params[in]   TSiteCriteria     Transmitting site criteria
     * @params[in]   RSiteCriteria     Receiving site criteria
     * @params[in]   eps_dielect       Earth dielectric constant (rel. permittivity).
     * @params[in]   sgm_conductivity  Earth conductivity (Siemens/m).
     * @params[in]   eno               Atmospheric bending const, n-Units
     * @params[in]   frq_mhz           Frequency in MHz
     * @params[in]   radio_climate     Radio Climate
     * @params[in]   pol               Antenna polarization
     * @params[in]   pctTime           0.01 .. .99
     * @params[in]   pctLoc            0.01 .. .99
     * @params[in]   pctConf           0.01 .. .99
     * @params[out]  dbloss            Loss in DB
     * @params[out]  errnum            Error code
     */
    void area(const VariabilityMode ModVar,
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
            Error &errnum);
private:

    enum class ControllingMode : int8_t {
        PointToPoint = static_cast<int8_t>(-1),   /*!< Point to point */
        AreaContinuation,                         /*!< Area Continuation */
         StartOfArea                               /*!< Start of Area */
    };

    enum class ControlSwitch : uint8_t {
        NoRecalculationNeeded,      /*!< 0: no recalculations needed */
        DistanceChanged,            /*!< 1: distance changed */
        AntennaHeightsChanged,      /*!< 2: antenna heights changed */
        FrequencyChanged,           /*!< 3: frequency changed */
        MdvarChanged,               /*!< 4: mdvar changed */
        ClimateChangedOrInitialise  /*!< 5: climate changed, or initialize everything */
    };

    static constexpr float64_t third = 1.0 / 3.0;  /*!< One third constant */
    static constexpr float64_t f_0   = 47.7;       /*!< 47.7 MHz from [Alg 1.1], to convert frequency into wavenumber and vica versa */

    struct prop_type {
        // General input
        float64_t d;           /*!< distance between the two terminals */
        float64_t h_g[2];      /*!< antenna structural heights (above ground) */
        float64_t k;           /*!< wave number (= radio frequency) */
        float64_t delta_h;     /*!< terrain irregularity parameter */
        float64_t N_s;         /*!< minimum monthly surface refractivity, n-Units */
        float64_t gamma_e;     /*!< earth's effective curvature */
        float64_t Z_g_real;    /*!< surface transfer impedance of the ground */
        float64_t Z_g_imag;

        // Additional input for point-to-point mode
        float64_t h_e[2];      /*!< antenna effective heights */
        float64_t d_Lj[2];     /*!< horizon distances */
        float64_t theta_ej[2]; /*!< horizontal elevation angles */
        ControllingMode mdp;   /*!< controlling mode: -1: point to point, 1 start of area, 0 area continuation */

        // Output
        Error kwx;           /*!< error indicator */
        float64_t A_ref;       /*!< reference attenuation */

        // used to be propa_type, computed in lrprop()
        float64_t dx;          /*!< scatter distance */
        float64_t ael;         /*!< line-of-sight coefficients */
        float64_t ak1;         /*!< dito */
        float64_t ak2;         /*!< dito */
        float64_t aed;         /*!< diffraction coefficients */
        float64_t emd;         /*!< dito */
        float64_t aes;         /*!< scatter coefficients */
        float64_t ems;         /*!< dito */
        float64_t d_Lsj[2];    /*!< smooth earth horizon distances */
        float64_t d_Ls;        /*!< d_Lsj[] accumulated */
        float64_t d_L;         /*!< d_Lj[] accumulated */
        float64_t theta_e;     /*!< theta_ej[] accumulated, total bending angle */
    };

    /**
     * :27: Variablility parameters, page 13
     */
    struct propv_type {
        // Input:
        ControlSwitch lvar;  /*!< control switch for initialisation and/or recalculation
         * 0: no recalculations needed
         * 1: distance changed
         * 2: antenna heights changed
         * 3: frequency changed
         * 4: mdvar changed
         * 5: climate changed, or initialize everything */
        int32_t mdvar;       /*!< desired mode of variability */
        RadioClimate klim;   /*!< climate indicator */
        // Output
        float64_t sgc;   /*!< standard deviation of situation variability (confidence) */
    };

    /**
     * This performs the FORTRAN DIM function.
     *
     * @return x-y if x is greater than y; otherwise result is 0.0
     */
    float64_t FORTRAN_DIM(const float64_t x, const float64_t y);

    /**
     * :13: single-knife attenuation, page 6
     *
     * The attenuation due to a single knife edge -- this is an approximation of a Fresnel integral as a function of v^2.
     * The non-approximated integral is documented as [Alg 4.21]
     *
     * Now, what is "single knife attenuation"? Googling found some paper http://www.its.bldrdoc.gov/pub/ntia-rpt/81-86/81-86.pdf,
     * which actually talks about multi-knife attenuation calculation. However, it says there that single-knife attenuation
     * models attenuation over the edge of one isolated hill.
     *
     * Note that the arguments to this function aren't v, but v^2
     */
    float64_t Fn(const float64_t v_square);

    /**
     * :14: page 6
     *
     * The heigh-gain over a smooth spherical earth -- to be used in the "three radii" mode. The approximation is that given in [Alg 6.4ff].
     */
    float64_t F(const float64_t x, const float64_t K);

    /**
     * :25: Tropospheric scatter frequency gain, [Alg 6.10ff], page 12
     */
    float64_t H_0(const float64_t r, const float64_t et);

    /**
     * :25: This is the F(\Theta d) function for scatter fields, page 12
     */
    float64_t F_0(const float64_t td);

    /**
     * :10: Diffraction attenuation, page 4
     *
     * The function adiff finds the "Diffraction attenuation" at the distance s. It uses a convex combination of smooth earth
     * diffraction and knife-edge diffraction.
     */
    float64_t adiff(const float64_t s, prop_type &prop);

    /**
     * :22: Scatter attenuation, page 9
     *
     * The function ascat finds the "scatter attenuation" at the distance d. It uses an approximation to the methods of NBS
     * TN101 with check for inadmissable situations. For proper operation, the larger distance (d = d6) must be the first
     * called. A call with d = 0 sets up initial constants.
     *
     * One needs to get TN101, especially chaper 9, to understand this function.
     */
    float64_t A_scat(const float64_t s, prop_type &prop);


    float64_t abq_alos(const std::complex<float64_t> &r);

    /**
     * :17: line-of-sight attenuation, page 8
     *
     * The function alos finds the "line-of-sight attenuation" at the distance d. It uses a convex combination of plane earth
     * fields and diffracted fields. A call with d=0 sets up initial constants.
     */
    float64_t A_los(const float64_t d, prop_type &prop);

    /**
     * :5: LRprop, page 2
     *
     * The value mdp controls some of the program flow. When it equals -1 we are in point-to-point mode, when 1 we are beginning
     * the area mode, and when 0 we are continuing the area mode. The assumption is that when one uses the area mode, one will
     * want a sequence of results for varying distances.
     */
    void lrprop(const float64_t d, prop_type &prop);

    void qlra(const SiteCriteria *const kst, const RadioClimate klimx, const VariabilityMode mdvarx, prop_type &prop, propv_type &propv);

    /**
     * :51: Inverse of standard normal complementary probability
     *
     * The approximation is due to C. Hastings, Jr. ("Approximations for digital computers," Princeton Univ. Press, 1955) and the
     * maximum error should be  4.5e-4.
     */
    float64_t qerfi(const float64_t q);

    /**
     * :41: preparatory routine, page 20
     *
     * This subroutine converts
     *   frequency @fmhz
     *   surface refractivity reduced to sea level @en0
     *   general system elevation @zsys
     *   polarization and ground constants @eps, @sgm
     * to
     *   wave number @wn,
     *   surface refractivity @ens
     *   effective earth curvature @gme
     *   surface impedance @zgnd
     *
     * It may be used with either the area prediction or the point-to-point mode.
     */
    void qlrps(const float64_t fmhz, const float64_t zsys, const float64_t en0, const Polarization ipol, const float64_t eps, const float64_t sgm, prop_type &prop);

    /**
     * :30: Function curv, page 15
     */
    float64_t curve(float64_t const c1, float64_t const c2, float64_t const x1, float64_t const x2, float64_t const x3, float64_t const de);


    /**
     *  :28: Area variablity, page 14
     */
    float64_t avar(const float64_t zzt, const float64_t zzl, const float64_t zzc, prop_type &prop, propv_type &propv);

    /*
     * :45: Find to horizons, page 24
     *
     * Here we use the terrain profile @pfl to find the two horizons. Output consists of the horizon distances @dl and the horizon
     * take-off angles @the. If the path is line-of-sight, the routine sets both horizon distances equal to @dist.
     */
    void hzns(const float64_t *const pfl, prop_type &prop);

    /**
     * :53: Linear least square fit, page 28
     */
    void zlsq1(const float64_t *const z, const float64_t x1, const float64_t x2, float64_t& z0, float64_t& zn);

    /**
     * :52: Provide a quantile and reorders array @a, page 27
     */
    float64_t qtile(const int32_t nn, float64_t *const a, const int32_t ir);

    /**
     * :41: Prepare model for point-to-point operation, page 22
     *
     * This mode requires the terrain profile lying between the terminals. This should be a sequence of surface elevations at
     * points along the great circle path joining the two points. It should start at the ground beneath the first terminal and
     * end at the ground beneath the second. In the present routine it is assumed that the elevations are equispaced along the
     * path. They are stored in the array @pfl along with two defining parameters.
     *
     * We will have:
     *   pfl[0] = np, the number of points in the path
     *   pfl[1] = xi, the length of each increment
     */
    void qlrpfl(const float64_t *const pfl, const RadioClimate klimx, const int32_t mdvarx, prop_type &prop, propv_type &propv);

    /**
     * :48: Find interdecile range of elevations, page 25
     *
     * Using the terrain profile @pfl we find \Delta h, the interdecile range of elevations between the two points @x1 and @x2.
     */
    float64_t dlthx(const float64_t *const pfl, const float64_t x1, const float64_t x2);
};

#endif /* LONGLEYRICEITM_H_ */
