/* Structures */

#define DISK_HORIZON_MU 51
#define DISK_STRUCTURE 10001
#define DIM 4 
#define MaxscatD 100

typedef struct 
{
    double M,a,epsilon3,Mdot;
    
	// important for ray-tracing
    int n_horizon_mu;
    double horizon[DISK_HORIZON_MU]; // event horizon as function of mu [M]

    // important for modeling of emission and scattering
    double ISCO; // [M]
    double r1,r2,rmax;
    
    int n_disk_structure;

    double u     [DISK_STRUCTURE][DIM]; // 4-velocity of disk element
    double emui  [DISK_STRUCTURE][DIM][DIM];
    double eimu  [DISK_STRUCTURE][DIM][DIM];
    
    // radial integration

//    double r1      [DISK_SAMPLE]; // lower end of bin
//    double r2      [DISK_SAMPLE]; // upper end of bin

    double r       [DISK_STRUCTURE]; // center of bin
    double r_T     [DISK_STRUCTURE]; // T in bin
    double r_weight[DISK_STRUCTURE]; // weight of bin (proportional of flux from bin)
    
    double rCorona_ends;
    double rCorona_start;
    int tedad;
    int nn;
    
} Disk, *Disk_p;


typedef struct 
{
    /* emission position and angle */
    double r0;
    double T; // temperature in eV */
    double u0PF[DIM]; // Plasma Frame 4-velocity
    double u0BL[DIM]; // Boyer Lindquist 4-velocity
    double E0;
    double E,L,b;
    
    /* tracking variables */
    double xBL[DIM]; // BL 4-coord
    double uBL[DIM]; // BL 4-velocity
    double fBL[DIM]; // polarization 4-vector  
    double pol;

    double xSC[DIM];
    double uSC[DIM];
    
    double uCS[DIM];  // 4-wave vector in coordinate-stationary [CS] reference frame
    double fCS[DIM];  // polarization  in coordinate-stationary [CS] reference frame
    double NuCS,NfCS; // parameters to monitor the accuracy of the ray tracing
    double stokes[3]; // Stokes vector in CS frame.

    double weight; // weight
    int    nScatter; // number of scatterings
    int   timeout;
    
    double weight_ep;  //BB
    int timeout_corona;
    int nScatter_Corona;
    double E_b_a_scat_Disk_BL[MaxscatD][5];
    double weight_ep_nopol;
    
    int totalN;
} Photon, *Photon_p;


typedef struct 
{
    /* emission position and angle */
    float r0; // Boyer Lindquist radial coordinate r0 from which the photon was launched
    float T; // temperature of the disk where the photon was launched in eV */
    float u0PF[DIM]; // Plasma Frame 4-velocity with which the photon was launched
    float E0; // energy at infinity when the photon was launched
    float E; // energy at infinity when the photon arrives; may differ from E0 owing to scatterings
    
    /* tracking variables */
    float xBL[DIM]; // Boyer Lindquist 4-vector of the final position of the photon
    float uCS[DIM];  // 4-wave vector in coordinate-stationary [CS] reference frame
    
    float stokes[3]; // Stokes vector in CS frame.
    int nScatter; // number of scatterings
    int timeout;
    
    
    int nScatter_Corona;

    double weight_ep_nopol;

    
} PhotonRed, *PhotonRed_p;

struct ChristoffelParams {
    double a_2;
    double a_3;
    double a_4;
    double a_6;
    double a_8;
    double M_2;
    double M_3;
    double M_4;
    double M_6;
    double M_7;
    double hair_2;
    double y1;
    double y1_2;
    double y1_3;
    double y1_4;
    double y1_5;
    double y1_6;
    double y1_7;
    double y1_8;
    double y1_9;
    double y2;
    double sin_y2;
    double sin_y2_2;
    double sin_y2_3;
    double sin_y2_4;
    double sin_2y2;
    double sin_2y2_2;
    double cos_y2;
    double cos_y2_2;
    double cos_y2_3;
    double cos_y2_4;
    double cos_y2_6;
    double cos_y2_8;
    double cos_y2_10;
    double cos_2y2;
    double cos_4y2;
    double pow2_y1_2_p_a_2cos_y2_2;
    double pow3_y1_2_p_a_2cos_y2_2;
    double pow4_y1_2_p_a_2cos_y2_2;
    double pow5_y1_2_p_a_2cos_y2_2;
    double pow2_y2_2_p_a_2cos_y2_2;
};


/* === === === === === === === === === === === === === === */

void loadConfig(int nConfig, Disk *disk_p);
void generatePhoton(Disk *disk_p, int ir,Photon *res_p);

// void trackPhoton(Disk *disk_p,Photon *gamma_p, double accuracy1, double accuracy2,int flag);
void trackPhoton(Disk *disk_p,Photon *gamma_p, double accuracy1, double accuracy2,std::ofstream &os,int flag);
void tab24(double mu, double &Is, double &pol);
void tab25(double *F, double mu0, double phi0, double *I, double mu, double phi);


double ChristoffelKerr(int rho, int mu, int nu, double M, double a, double hair,const double *y, const ChristoffelParams &p);
double metric(int mu, int nu, double M, double a, double epsilon3,double *y);
void CS_emui(double M, double a, double epsilon3,double *y,double emui[4][4]);


void load_table(double *t_x,int &numbins, char *fname1,int checknum);
double get_value_log(double *t_x,double *t_v,int numbins,double x);
double get_value_lin(double *t_x,double *t_v,int numbins,double x);

double ran1(long *idum);

int linBin(double x1,double x2,int numbins,double x);
int logBin(double x1,double x2,int numbins,double x);
double binLog(double x1,double x2,int numbins,int nbin);
double cot(double x);

void Matrix_SMult(double c,double a1[][3], double a2[][3]);
void Matrix_Add(double a1[][3], double a2[][3], double a3[][3]);
void Matrix_Mult(double a1[][3], double a2[][3], double a3[][3]);
double dot3(double a1[][3], double a2[][3], int aRow, int bCol);
void Matrix_Print(double *F,int n, int l);
void MatrixVector_Mult(double a1[][3], double v[3], double r[3]);
double dot(double *a, double *b);
void normalize(double *v,double *n);
void perp(double *v,double *k,double *p);
double angle(double *aNorm,double *bNorm);
void cross(double *a, double *b, double *c);
double getChi(double q, double u);

//BB

double RootFinder(double THETA,double lambda0); //BB
void BZ_emunu(double M, double a, double epsilon3, double *y ,double emunu[4][4], double enumu[4][4]); //BB
void Lorentz(double gamma_e, double beta_e, double *n ,double LAMBDA[4][4], double LAMBDAi[4][4]); //BB
double Gamma_MB(double Tcorona);  //BB
void Scatter(Photon *gamma_p, double u[4], double f_corona[4],double Tcorona, double emunu[4][4], double enumu[4][4], double fac);
void PF_eimu(double M, double a, double epsilon3, double *y ,double eimu[4][4], double emui[4][4]); //BB // transforms BL to PF and vsv
void BP_eimu(double M, double a, double epsilon3, double *y ,double eimuBP[4][4], double emuiPB[4][4]); //BB // emunu transforms BL to ZAMO, enumu, inverse metrix, transfroms ZAMO to BL

// Planck functions used:

double planck_sample ( double a, double b, int &seed );
int zipf_sample ( double a, int &seed );
double gamma_sample ( double a, double b, double c, int &seed );
double r8_uniform_01 ( int &seed );
double exponential_01_sample ( int &seed );
double normal_01_sample ( int &seed );
double r8_abs ( double x );
double r8_sign ( double x );


#define NMAX 4
void migs (double a[][NMAX],int n,double x[][NMAX]);
void elgs (double a[][NMAX],int n,int indx[NMAX]);







