/* Mathematical Constants */
#define LOG10 2.302585093
#define PI    3.141592653589793
#define RtD   57.295779513082325
/* Physical Constants */

/* Gauss Units */
#define CGS_c      2.997925E+10   /* cm/s */
#define CGS_me     9.109E-28      /* g */
#define CGS_me_erg 8.1760E-7     /* erg */
#define CGS_mp     1.6726E-24     /* g */
#define CGS_e      4.8066E-10     /* statcoul  */
#define CGS_k      1.3807E-16     /* erg/K */
#define CGS_re     2.818E-13      /* cm */
#define CGS_pc     3.0856E+18     /* cm */
#define CGS_h      6.6261E-27     /* erg s */
#define CGS_hbar   1.054589E-27  /* erg s */
#define CGS_tcs    6.653E-25      /* Thompson Cross section in cm**2 */
#define CGS_SBK    5.5670400e-5   // erg s-1 m-2 -4

#define CGS_MSUN   1.98892E+33      // Solar mass in grams
#define CGS_G      6.67300E-8       // Gravitational Constant in dyne cm^2 gm-2

/* SI Units */
#define SI_c    2.997925E+08   /* m/s */
#define SI_h    6.6261E-34     /* J s */
#define SI_hbar 1.05457657E-34 /* J s */
#define SI_me   9.109E-31      /* kg */
#define SI_mp   1.6726E-27     /* kg */
#define SI_e    1.6022E-19     /* C  */
#define SI_k    1.3807E-23     /* J/K */
#define SI_re   2.818E-15      /* m */
#define SI_e0   8.8542E-12     /* C**2 / m**2 N */
#define SI_pc   3.0856E+16     /*  m */
#define SI_tcs  6.653E-29      /* Thompson Cross section in m**2 */
#define SI_SBK  5.5670400e-8   // W m-2 K-4

#define SI_MSUN 1.98892E+30 // kg
#define SI_G	6.67300E-11 // m+3 kg-1 s-2

#define C_me_eV 511E+3         /* eV */
#define NU_eV   4.136E-15
#define NU_TeV  4.136E-27
#define eV_NU   2.41779E+14
#define NU_erg  6.6176E-27
#define J_eV    6.242197253e+18
#define eV_J    1.602E-19
#define J_erg   1E+7           /* Joule in ergs */
#define erg_J   1E-7
#define eV_erg  1.6E-12  
#define erg_eV  6.25E+11  

static double sqrarg;
#define sq(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double cubearg;
#define cube(a) ((cubearg=(a)) == 0.0 ? 0.0 : cubearg*cubearg*cubearg)

static double maxarg1,maxarg2;
#define maximum(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static double minarg1,minarg2;
#define minimum(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static double signA,signB;
#define SIGN(a,b) (signA=a,signB=b,(signB) >= 0.0 ? fabs(signA) : -fabs(signA))
