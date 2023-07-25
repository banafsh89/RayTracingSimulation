// This code is based of the raytracing code created by Henric Krawczynski on 5/07/12.
//  Copyright (c) 2014 Washington University in St. Louis. All rights reserved.

//
// Compile with g++ --std=c++11 -lm -O3 spAGN.cpp -o spAGN
// spAGN [config-#] [events-per-bin-#] [seed] [output-file-name]
//

// Description:
// The code simulates photons emerging from a Novikov/Thorne accretion disk with
// a polarization given by the classical Chandrasekhar result.
//
// The real work is done in the function "trackPhoton", where the photon trajectories are
// calculated, and the photon parameters (time, position, energy, momentum, polarization vector)
// are updated, and photons are scattered if they collide with the accretion disk or with Corona.

// Corona is extended from ISCO to rcorona and the accretion disk is extended after the corona radius. The scattering in corona is simulated in function "scatter". Here we assume that corona is not emitting seed photons and only up-scatter and energize the photons.
//Size of the corona is set in the LoadConfig function. You can change the size manually in that function. You can also change the tau0 (optical depth of the corona) by searching for tau0 and manually change it. ( I have not yet chance to put it on the global format to change it easier).
//The source photons are all generated from the disk. 

/* Here is what the most important functions do:
 
 loadConfig: loads the information about the parameters
 (black hole mass, accretion rate, black hole spin)
 and the radial disk emissivity.
 
 generatePhoton(Disk *disk_p, int ir,Photon *res_p):
 emits a photon in a radial bin into a random direction, makign use of Chandrasekhar's
 results concenring the mean polarization and elevation-dependent intensity of the emission
 
 trackPhoton: does the real work of propagating the photon and scattering it.
 Scatter: calculate the scattering, transform the frames, and calculate the scattered photon polarisation and energy.
 tab24, tab25: have Chandrasekhar's information about the initial polarization
 of photons emerging from the disk (tab24) and how scattering changes the polarization (tab25).
 
 Scatter: It does the scattering of photon inside corona and transorm the frames from BL to ZAMO and electron rest frame.
 
 Christoffel: the code can do all kinds of metrics.
 The parallel transport requires the Christoffel symbols.
 
 */


// C++
#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <iomanip>   // format manipulation
#include <string>
#include <cmath>
#include <vector>

// C
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>    

using namespace std;

#include "jpsch.h"
#include "macro.h"

#define NO_DEBUG
#define NO_DEBUG_SCAT
#define NO_DEBUG_FINAL
#define NO_PRINT

inline double sqr(double x) { return x*x; }

inline int make_christoffel_index(int rho, int mu, int nu)
{
    return nu+(4*(mu+4*rho));
}

static long seed = -90;
static int trial=0;
int Pseed=10;

#define MAXT 100000000
int main(int argc, char** argv)
{
    
    /* Aim: produce num events per radial bin and store them in a file */
    if (argc<=1) 
	{
        /* Help */
        cout << "Compile with g++ --std=c++11 -lm -O5 spAGN.cpp -o spAGN"<<endl;
        cout << "Run     with ./spAGN [config-#] [events-per-bin-#] [seed] [output-file-name]"<<endl;
    }
    
    /* load disk information */ 
    Disk *disk_p=new Disk;
    
    int nConfig=1;
    if (argc>=2) nConfig = (double)atoi(argv[1]);
    
    loadConfig(nConfig,disk_p);
    
    int nBin=1;
    if (argc>=3) nBin = (double)atoi(argv[2]);
    cout <<endl<<"Simulating "<<nBin<<" events for all "<<DISK_STRUCTURE<<" radial bins."<<endl; 
    
    if (argc>=4) seed = -atoi(argv[3]);
    cout <<"Seed "<<seed<<endl;
    /* Loop over radial bins */
    
    cout<<endl<<" ISCO "<<disk_p->ISCO;
    cout<<" r1 "<<disk_p->r1;
    cout<<" r2 "<<disk_p->r2;
    cout<<" rmax "<<disk_p->rmax<<endl<<endl;
    
    int np=0;
    char  fname[100];
    
    if (argc>=5) sprintf(fname,"%s.dat",argv[4]);
    else sprintf(fname,"out.dat");
    
    ofstream o_file (fname,ios::out | ios::binary);   //open file for writing
    if ( !o_file.is_open()) {cout << "couldn't open output file" << fname<<" "<<endl; exit(1);}
    cout <<"Opened output file "<<fname<<endl;
    
    Photon    *gamma_p    = new Photon;
    PhotonRed *gammaRed_p = new PhotonRed;
    int size = sizeof(PhotonRed);
    
    cout <<"Writing data in chunks of "<<size<<" units."<<endl;
    int total=0;
    int flag;
    
    char track[100];
    sprintf(track, "track.dat");
    ofstream os(track);
    
    if( !os )
    {
        cout << "error opening" <<track<< endl;
        exit(-1);
    }
    
    disk_p->tedad=0;
    disk_p->nn=0;
    
    int ircorona=logBin(disk_p->r1,disk_p->r2,disk_p->n_disk_structure,disk_p->rCorona_ends);

    
    /* Loop over radial bins of the accretion disk */
    for (int ir=0; ir<ircorona; ir++) {
        /* In each radial bin, generate nBin photons */
        for (int ip=0; ip<nBin; ip++) {
            
            gamma_p->totalN=total;

            generatePhoton(disk_p,ir,gamma_p);
            
            
            flag=-1;

#ifdef PRINT
            if (total<100) {
                flag=total;
                cout <<"Writing track!"<<total<<endl;
            }
#endif
            
            trackPhoton(disk_p,gamma_p,1E-5,1E-5,os,flag); /* includes scattering */
            
            /* save data t binary file */
            gammaRed_p->r0       = gamma_p->r0;
            gammaRed_p->T        = gamma_p->T;
            gammaRed_p->E0       = gamma_p->E0;
            gammaRed_p->E        = gamma_p->E;
            gammaRed_p->nScatter = gamma_p->nScatter;
            gammaRed_p->timeout  = gamma_p->timeout;
            gammaRed_p->nScatter_Corona=gamma_p->nScatter_Corona;
            gammaRed_p->weight_ep_nopol=gamma_p->weight_ep_nopol;

            

            for (int i=0;i<3;i++) gammaRed_p->stokes[i]=gamma_p->stokes[i];
            
            for (int i=0;i<4;i++){
                gammaRed_p->u0PF[i] = gamma_p->u0PF[i];
                gammaRed_p->xBL[i]  = gamma_p->xBL[i];
                gammaRed_p->uCS[i]  = gamma_p->uCS[i];

            }
            
            o_file.write ((char*)gammaRed_p,size);
            if (total%5000==0)
            cout <<"#" <<total<<" ir "<<ir<<" r0 "<<gammaRed_p->r0<<" T "<<gammaRed_p->T<<" E "<<gammaRed_p->E <<" nS "<<gammaRed_p->nScatter<<" xBL "<<gammaRed_p->xBL[0]<<" "<<gammaRed_p->xBL[1]<<" "<<gammaRed_p->xBL[2]<<" "<<gammaRed_p->xBL[3]<<" "<<" nto "<<gammaRed_p->timeout<<" Nscat "<<gammaRed_p->nScatter<<" NCoronaScat "<<gammaRed_p->nScatter_Corona<<endl;
           
 
            if (++total>MAXT) break; 
        }
        if (total>MAXT) break; 
    }
    
    o_file.close();
    os.close();
    cout<<" dist2>0.1 "<<disk_p->tedad<<" u<0 "<<disk_p->nn<<endl;
    cout <<"Done."<<endl;
}

/* custom functions */

void loadConfig(int nConfig, Disk *disk_p)
{
    char root[20][30]={
    "Parm_1_0p0_0",
    "Parm_1_0p5_0",    
    "Parm_1_0p9_0",
    "Parm_1_0p99_0", // 3    
    "Parm_1_0p5_0a25",
    "Parm_1_0p5_0a4",
    "Parm_1_0p5-2p5",
    "Parm_1_0p5+2p5", 
    "Parm_1_0p5-5",
    "Parm_1_0p5+5", // 9
    "Parm_1_0p99-5", // 10 
    "Parm_1_0p99-2p5",
    "Parm_1_0p5+6p32915",// 12
    "Parm_1_0p5-30p611"}; //13
    
    char emui[6][30]={"lmi00","lmi03","lmi11","lmi22","lmi30","lmi33"};
    char eimu[6][30]={"lim00","lim03","lim11","lim22","lim30","lim33"};
    
    int a[4][4]={  
        0,-1,-1,1,
        -1,2,-1,-1,
        -1,-1,3,-1,
        4,-1,-1,5};
    
    
    char fname[50]; 
    int dummy;
    
    sprintf(fname,"Data/%s/horizon.dat",root[nConfig]);
    load_table(disk_p->horizon,disk_p->n_horizon_mu,fname,DISK_HORIZON_MU);
    
    double in[10];
    sprintf(fname,"Data/%s/isco.dat",root[nConfig]);
    load_table(in,dummy,fname,9);
    
    disk_p->M       =in[0];
    disk_p->a       =in[1];
    disk_p->epsilon3=in[2];
    
    disk_p->ISCO=in[3];
    disk_p->r1  =in[4];
    disk_p->r2  =in[5];
    disk_p->rmax=in[6];
        
    disk_p->n_disk_structure=in[7];
    disk_p->Mdot=in[8];
    cout <<"M "<<disk_p->M<<" a "<<disk_p->a<<" e3 "<<disk_p->epsilon3<<" Mdot "<<disk_p->Mdot<<endl;
    
    //radius of corona should be set here
    disk_p->rCorona_ends=15; //in r_g unit
    disk_p->rCorona_start=disk_p->ISCO;    //starts farther away from ISCO to prevent many photons traped in BH
    
    sprintf(fname,"Data/%s/Tlist.dat",root[nConfig]);
    load_table(disk_p->r_T,dummy,fname,DISK_STRUCTURE);
    
    sprintf(fname,"Data/%s/Wlist.dat",root[nConfig]);
    load_table(disk_p->r_weight,dummy,fname,DISK_STRUCTURE);
    
    double inp1[6][DISK_STRUCTURE];
    double inp2[6][DISK_STRUCTURE];
    
    for (int i=0;i<6;i++){
        sprintf(fname,"Data/%s/%s.dat",root[nConfig],emui[i]);
        load_table(inp1[i],dummy,fname,DISK_STRUCTURE);
        sprintf(fname,"Data/%s/%s.dat",root[nConfig],eimu[i]);
        load_table(inp2[i],dummy,fname,DISK_STRUCTURE);
    };
    
    for (int n=0; n<DISK_STRUCTURE; n++)
        for (int i=0;i<4;i++)
            for (int j=0;j<4;j++)
                if (a[i][j]>=0) {
                    disk_p->emui[n][i][j]=inp1[a[i][j]][n];
                    disk_p->eimu[n][i][j]=inp2[a[i][j]][n];
                } else {
                    disk_p->emui[n][i][j]=0.;
                    disk_p->eimu[n][i][j]=0.;
                };
}

void generatePhoton(Disk *disk_p, int ir,Photon *gamma_p)
{
    double y[4];
       
    gamma_p->r0=binLog(disk_p->r1,disk_p->r2,disk_p->n_disk_structure,ir);
    gamma_p->T = disk_p->r_T[ir];
    
    gamma_p->timeout=0;
    double amu;
    
    if(gamma_p->r0<=disk_p->rCorona_ends)
        amu=1-2*ran1(&seed);
    else
        amu=ran1(&seed);

    
    double phi=2. * PI * ran1(&seed);
    
    double hardening=1.8;
    double A =SI_h*J_eV/(1000*gamma_p->T*hardening);  // Temprature is in unit of keV, So I multiplied by 1000 to make it in unit of eV
    
    
    double freq=planck_sample(A,2,Pseed);
    double Ei = SI_h*J_eV*freq/1000; //keV  h*\nu \nu is from BB radiation
    
    gamma_p->u0PF[0] =  Ei*1.;
    gamma_p->u0PF[1] =  Ei*sqrt(1.-sq(amu))*cos(phi); // r
    gamma_p->u0PF[2] =  Ei*(-amu); // theta ==> emission into the upper hemisphere!
    gamma_p->u0PF[3] =  Ei*sqrt(1.-sq(amu))*sin(phi); // phi
    
    
   
    
    if (gamma_p->r0<=disk_p->rCorona_ends) {
        y[0]=0.;
        y[1]=gamma_p->r0;
        double mu=1-2*ran1(&seed);
        y[2]=acos(mu);
        y[3]=0;
    }
    else{
        y[0]=0.;
        y[1]=gamma_p->r0;
        y[2]=PI/2.;
        y[3]=0.;
    }

    double emunu[4][4];
    double enumu[4][4];
    
    BZ_emunu(disk_p->M,disk_p->a,disk_p->epsilon3,y,emunu,enumu); //emunu transforms BL to Zamo and enumu transforms back
    
    if (gamma_p->r0<=disk_p->rCorona_ends) {
        for (int mu=0;mu<4;mu++){
            gamma_p->u0BL[mu]=0.;
            for (int i=0; i<4; i++) gamma_p->u0BL[mu]+=emunu[mu][i]*gamma_p->u0PF[i];
            gamma_p->uBL[mu]=gamma_p->u0BL[mu];
        }

    }
    else{
        for (int mu=0;mu<4;mu++){
            gamma_p->u0BL[mu]=0.;
            for (int i=0; i<4; i++) gamma_p->u0BL[mu]+=disk_p->emui[ir][mu][i]*gamma_p->u0PF[i];
            gamma_p->uBL[mu]=gamma_p->u0BL[mu];
        }
    }
    
#ifdef DEBUG
    cout << "gamma_p->u0BL ";
    for (int i=0;i<4;i++) cout<<gamma_p->u0BL[i]<<" ";
    cout <<endl<<endl;
    
    cout << "emui "<<endl;
    for (int i=0;i<4;i++) 
    {
        for (int j=0;j<4;j++) 
            cout<< disk_p->emui[ir][i][j] <<" ";
        cout <<endl;
    }

    double dist1=0.;
    
    for (int i=0; i<4; i++) 
        for (int j=0; j<4; j++) 
            dist1 += metric(i,j,disk_p->M,disk_p->a,disk_p->epsilon3,y)*gamma_p->uBL[i]*gamma_p->uBL[j];
    
    cout <<"norm of velocity:"<<dist1<<endl<<endl;

#endif


    gamma_p->xBL[0]=y[0];
    gamma_p->xBL[1]=y[1]; // r
    gamma_p->xBL[2]=y[2]; // theta
    gamma_p->xBL[3]=y[3]; // phi
    
    double g00=metric(0,0,disk_p->M,disk_p->a,disk_p->epsilon3,y);
    double g03=metric(0,3,disk_p->M,disk_p->a,disk_p->epsilon3,y);
    double g11=metric(1,1,disk_p->M,disk_p->a,disk_p->epsilon3,y);
    double g22=metric(2,2,disk_p->M,disk_p->a,disk_p->epsilon3,y);
    double g33=metric(3,3,disk_p->M,disk_p->a,disk_p->epsilon3,y);
    
    gamma_p->E = gamma_p->E0 = -1.*(g00*gamma_p->u0BL[0]+g03*gamma_p->u0BL[3]);
    gamma_p->L =     (g33*gamma_p->u0BL[3]+g03*gamma_p->u0BL[0]);
    gamma_p->b = gamma_p->L/gamma_p->E;
    
    
    double fPF[DIM];// f should point along the y-direction for phi=0.
    fPF[0]=0.;
    fPF[1]=-sin(phi); // e_r     == x
    fPF[2]=0.;        // e_theta == z
    fPF[3]=cos(phi);  // e_phi   == y
    
    
    
    
    if (gamma_p->r0<=disk_p->rCorona_ends) {
        for (int mu=0;mu<4;mu++){
            gamma_p->fBL[mu]=0.;
            for (int i=0; i<4; i++) gamma_p->fBL[mu]+=emunu[mu][i]*fPF[i];
        }
    }
    else{
        for (int mu=0;mu<4;mu++){
            gamma_p->fBL[mu]=0.;
            for (int i=0; i<4; i++) gamma_p->fBL[mu]+=disk_p->emui[ir][mu][i]*fPF[i];
        }
    }
    
    double weight, pol;
    tab24(amu,weight,pol);
    
    gamma_p->weight  =weight*disk_p->r_weight[ir];
    gamma_p->pol     =pol;
    
    if (gamma_p->r0<=disk_p->rCorona_ends) {
        gamma_p->weight=1*disk_p->r_weight[ir];
        gamma_p->pol=0;
    }
    
#ifdef DEBUG_SCAT
    cout <<endl<<"Initialization :"<<" E "<<gamma_p->E<<" L "<<gamma_p->L<<" b "<<gamma_p->b<<" polarization "<<gamma_p->pol<<endl;
    cout <<"uPF :"<<gamma_p->u0PF[0]<<" "<<gamma_p->u0PF[1]<<" "<<gamma_p->u0PF[2]<<" "<<gamma_p->u0PF[3]<<endl;
    cout <<"uBL :"<<gamma_p->uBL[0]<<" "<<gamma_p->uBL[1]<<" "<<gamma_p->uBL[2]<<" "<<gamma_p->uBL[3]<<endl<<endl;
#endif
    
    gamma_p->nScatter=0;
    gamma_p->nScatter_Corona=0;

}

void trackPhoton(Disk *disk_p,Photon *gamma_p, double accuracy1, double accuracy2,std::ofstream &os,int flag)
{
#define NUM_RK 10
    
    double y   [4][NUM_RK];
    double ydot[4][NUM_RK];
    double k   [4][NUM_RK];
    double g00,g03,g11,g22,g33,den;
    double dist1,dist2,dist3,dist4,dist5,ff,fff,fac1,fac2,base;
    
    double delta[NUM_RK],vel[4];
    
    double u[4],uPF[4],fPF[4],vPF[4],hPF[4];
    
    double eta[4]={-1.,1.,1.,1.};
    
    
    int    timeout=0;
    int    timeout_corona=0;
    gamma_p->timeout_corona=0;
    
    double w1[3]={0.5,0.5,1.};
    double w2[4]={6.,3.,3.,6.};
    
    double dh=0.1;
    double rh=0.;
    
    
    y[0][0] =gamma_p->xBL[0]; // t
    y[0][1] =gamma_p->xBL[1]; // r
    y[0][2] =gamma_p->xBL[2]; // theta
    y[0][3] =gamma_p->xBL[3]; // phi
    
    y[0][4] =gamma_p->uBL[1]/gamma_p->E; // d_r/d_lambda' = d_r/d_lambda * d_lambda/d_lambda' (with lambda=lambda'/E)
    y[0][5] =gamma_p->uBL[2]/gamma_p->E; // d_theta/d_lambda'
    
    y[0][6] =gamma_p->fBL[0]; // fBL[0]
    y[0][7] =gamma_p->fBL[1]; // fBL[1]
    y[0][8] =gamma_p->fBL[2]; // fBL[2]
    y[0][9] =gamma_p->fBL[3]; // fBL[3]
    
    for (int i=0;i<4;i++){
        gamma_p->xSC[i] = 0.;
        gamma_p->uSC[i] = 0.;
    }
    
    
    double t_prev=y[0][0];
    double r_prev=y[0][1];
    double theta_prev=y[0][2];
    double phi_prev=y[0][3];
    double y04_prev=y[0][4];
    double y05_prev=y[0][5];
    double dh_prev=0;
    int flagtest=0;
    int flagdisktest=0;
    double adjustment=1;
    double f_prev[4]={y[0][6],y[0][7],y[0][8],y[0][9]};

        
    do {
        
        double correct=1.;
        
        //Here while integrate the geodesic equation, I check to make sure if dtau in the corona after one step of photon is less than
        //0.1, otherwise we change the step size to make smaller steps in the corona.
        
        while (adjustment==1) {
            // cout<<" xposition before one step "<<y[0][0]<<" "<<y[0][1]<<" "<<y[0][2]<<" "<<y[0][3]<<endl;
            
            for (int iter=0; iter<4; iter++) {
                
                g00 = metric(0,0,disk_p->M,disk_p->a,disk_p->epsilon3,y[iter]);
                g03 = metric(0,3,disk_p->M,disk_p->a,disk_p->epsilon3,y[iter]);
                g33 = metric(3,3,disk_p->M,disk_p->a,disk_p->epsilon3,y[iter]);
                
                // (1) compute change of coordinates and velocities
                
                den = g33*g00-sq(g03);
                ydot[iter][0]=(-g33 - gamma_p->b*g03)/den; // dt/dlambda'
                ydot[iter][1]=y[iter][4]; // dr/dlambda'
                ydot[iter][2]=y[iter][5]; // dtheta/dlambda'
                ydot[iter][3]=(gamma_p->b*g00+g03)/den; // dphi/dlambda'
                
                
                //Fabian add
                double christoffel[64];
                const double M = disk_p->M;
                const double a = disk_p->a;
                const double hair = disk_p->epsilon3;
                const double *yiter = y[iter];
                ChristoffelParams p;
                p.a_2 = a*a;
                p.a_3 = p.a_2*a;
                p.a_4 = p.a_3*a;
                p.a_6 = p.a_4*p.a_2;
                p.a_8 = p.a_6*p.a_2;
                p.M_2 = M*M;
                p.M_3 = p.M_2*M;
                p.M_4 = p.M_3*M;
                p.M_6 = p.M_3*p.M_3;
                p.M_7 = p.M_6*M;
                p.hair_2 = hair*hair;
                p.y1 = yiter[1];
                p.y1_2 = p.y1*p.y1;
                p.y1_3 = p.y1_2*p.y1;
                p.y1_4 = p.y1_3*p.y1;
                p.y1_5 = p.y1_4*p.y1;
                p.y1_6 = p.y1_5*p.y1;
                p.y1_7 = p.y1_6*p.y1;
                p.y1_8 = p.y1_7*p.y1;
                p.y1_9 = p.y1_8*p.y1;
                p.y2 = yiter[2];
                p.sin_y2 = sin(p.y2);
                p.sin_y2_2 = p.sin_y2*p.sin_y2;
                p.sin_y2_3 = p.sin_y2_2*p.sin_y2;
                p.sin_y2_4 = p.sin_y2_3*p.sin_y2;
                p.cos_y2 = cos(p.y2);
                p.cos_y2_2 = p.cos_y2*p.cos_y2;
                p.cos_y2_3 = p.cos_y2_2*p.cos_y2;
                p.cos_y2_4 = p.cos_y2_3*p.cos_y2;
                p.cos_y2_6 = p.cos_y2_4*p.cos_y2_2;
                p.cos_y2_8 = p.cos_y2_6*p.cos_y2_2;
                p.cos_y2_10 = p.cos_y2_8*p.cos_y2_2;
                p.cos_2y2 = p.cos_y2_2 - p.sin_y2_2;     // cos(2x) = cos^2(x) - sin^2(x)
                p.cos_4y2 = 2*p.cos_2y2*p.cos_2y2 - 1;   // cos(4x) = 2cos^2(2x) - 1
                p.sin_2y2 = 2*p.sin_y2*p.cos_y2;         // sin(2x) = 2*sin(x)*cos(x)
                p.sin_2y2_2 = sqr(p.sin_2y2);
                p.pow2_y2_2_p_a_2cos_y2_2 = sqr(p.y1_2 + p.a_2*p.cos_y2_2);
                const double y1_2_p_a_2cos_y2_2 = p.y1_2 + p.a_2*p.cos_y2_2;
                p.pow2_y1_2_p_a_2cos_y2_2 = sqr(y1_2_p_a_2cos_y2_2);
                p.pow3_y1_2_p_a_2cos_y2_2 = p.pow2_y1_2_p_a_2cos_y2_2*y1_2_p_a_2cos_y2_2;
                p.pow4_y1_2_p_a_2cos_y2_2 = sqr(p.pow2_y1_2_p_a_2cos_y2_2);
                p.pow5_y1_2_p_a_2cos_y2_2 = p.pow4_y1_2_p_a_2cos_y2_2*y1_2_p_a_2cos_y2_2;
                
                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        for (int k = 0; k < 4; ++k) {
                            int x = make_christoffel_index(k, i, j);
                            christoffel[x] = ChristoffelKerr(k, i, j, M, a, hair, yiter, p);
                        }
                    }
                }
                
                
                ydot[iter][4]=
                -christoffel[make_christoffel_index(1, 0, 0)]    * sq(ydot[iter][0])
                -christoffel[make_christoffel_index(1, 1, 1)]    * sq(ydot[iter][1])
                -christoffel[make_christoffel_index(1, 2, 2)]    * sq(ydot[iter][2])
                -christoffel[make_christoffel_index(1, 3, 3)]    * sq(ydot[iter][3])
                -2.*christoffel[make_christoffel_index(1, 3, 0)]    * ydot[iter][3] * ydot[iter][0]
                -2.*christoffel[make_christoffel_index(1, 2, 1)]    * ydot[iter][2] * ydot[iter][1];
                
                ydot[iter][5]=
                -christoffel[make_christoffel_index(2, 0, 0)]    * sq(ydot[iter][0])
                -christoffel[make_christoffel_index(2, 1, 1)]    * sq(ydot[iter][1])
                -christoffel[make_christoffel_index(2, 2, 2)]    * sq(ydot[iter][2])
                -christoffel[make_christoffel_index(2, 3, 3)]    * sq(ydot[iter][3])
                -2.*christoffel[make_christoffel_index(2, 3, 0)]    * ydot[iter][3] * ydot[iter][0]
                -2.*christoffel[make_christoffel_index(2, 2, 1)]    * ydot[iter][2] * ydot[iter][1];
                
                // (2) compute change of polarization vector
                for (int i=0;i<4;i++){
                    
                    ydot[iter][6+i] =0;
                    
                    for (int gamma=0;gamma<4;gamma++)
                    for (int beta=0;beta<4;beta++)
                    {
                        ydot[iter][6+i]-= christoffel[make_christoffel_index(i, gamma, beta)]*y[iter][6+gamma]*ydot[iter][beta];
                    }
                }
                
                for (int i=0; i<NUM_RK; i++) {
                    k[iter][i] = dh * ydot[iter][i];
                    if (iter<3) y[iter+1][i]=y[0][i]+w1[iter]*k[iter][i];
                }
            }
            
            delta[1] = 0.;
            for (int iter=0; iter<4; iter++)
                delta[1] += k[iter][1]/w2[iter];
            if (y[0][1]+delta[1]>disk_p->rmax) correct = (disk_p->rmax+0.001-y[0][1])/delta[1];
            
            delta[3] = 0.;
            for (int iter=0; iter<4; iter++)
                delta[3] += k[iter][3]/w2[iter];
            
            if (delta[3]>0.01) correct = 0.01/delta[3];
            
            
            for (int j=0; j<NUM_RK; j++) {
                delta[j] = 0.;
                for (int iter=0; iter<4; iter++)
                    delta[j] += k[iter][j]/w2[iter];
                y[0][j] += correct*delta[j];
            }
            
            
            // in the ray-tracing we have to make sure that the photons have integration step which makes d\tau<<1 as is said in Schnittman& Krolik (2010-2013) paper. To insure that in each integration step we calculate this and if it is larger than one, we manually make the integration step for that step smaller. Here is how I do this.
            
            // check if photon is in corona s
            if ((y[0][1]>=disk_p->rCorona_start)&&(y[0][1]<=disk_p->rCorona_ends)&&(timeout_corona<2000)&&(flagtest==0)) // Here I am assuming that corona can have the Max scat of timeout_corona
            {
                
                double emunu[4][4];
                double enumu[4][4];
                
                BZ_emunu(disk_p->M,disk_p->a,disk_p->epsilon3,y[0],emunu,enumu);  //emunu transforms BL to Zamo and enumu transforms back
                
                
                //calculate BL displacement 4-vector
                double dispBL[4];
                dispBL[0] = y[0][0]-t_prev;
                dispBL[1] = y[0][1]-r_prev;
                dispBL[2] = y[0][2]-theta_prev;
                dispBL[3] = y[0][3]-phi_prev;
                
                double dispZAMO[4];
                for (int i=0;i<4;i++)
                {
                    dispZAMO[i]=0;
                    for (int mu=0; mu<4; mu++)
                    {
                        dispZAMO[i] += enumu[i][mu]*dispBL[mu];
                    }
                }
                
                
                //find spacial magnitude dl
                
                double dl = sqrt(sq(dispZAMO[1])+sq(dispZAMO[2])+sq(dispZAMO[3]));
                
                
                //calculate density and probability of scattering (Schnittman&Krolik 2010, pg.910)
                
                
                double Tcorona= 100; //Corona temp in kev
                double tau0 = 1.5; //vertical optical depth
                double alpha=tau0/(disk_p->rCorona_ends-disk_p->rCorona_start);
                if (alpha<0) {
                    cout<<"Wrong slope for tau"<<endl;
                    exit(1);
                }
                double dtau=dl*alpha;

                
                int dtaunew=dtau*100;
                
                if (dtaunew>11){
                    double scale= dtau/0.1;
                    
                    dh_prev=dh;
                    dh/=(scale);
                    
                    adjustment=1;
                    
                    //return all posison to previous one
                    y[0][0] =t_prev; // t
                    y[0][1] =r_prev; // r
                    y[0][2] =theta_prev; // theta
                    y[0][3] =phi_prev; // phi
                    y[0][4]=y04_prev;
                    y[0][5]=y05_prev;
                    y[0][6]=f_prev[0];
                    y[0][7]=f_prev[1];
                    y[0][8]=f_prev[2];
                    y[0][9]=f_prev[3];
                    
                }
                else
                    adjustment=0;
                
            }
            else
                adjustment=0;
            
        }
        
        
        
        // location and velocity: check accuracy and adapt step-size
        
        g00 = metric(0,0,disk_p->M,disk_p->a,disk_p->epsilon3,y[0]);
        g03 = metric(0,3,disk_p->M,disk_p->a,disk_p->epsilon3,y[0]);
        g11 = metric(1,1,disk_p->M,disk_p->a,disk_p->epsilon3,y[0]);
        g22 = metric(2,2,disk_p->M,disk_p->a,disk_p->epsilon3,y[0]);
        g33 = metric(3,3,disk_p->M,disk_p->a,disk_p->epsilon3,y[0]);
        
        den = g33*g00-sq(g03);
        u[0]  = (-g33 - gamma_p->b*g03)/den; // dt/dlambda'
        u[1]  = y[0][4];
        u[2]  = y[0][5];
        u[3]  = (gamma_p->b*g00+g03)/den; // dphi/dlambda'        
        dist2 = g00*sq(u[0]) + g33*sq(u[3]) + 2.*g03*u[0]*u[3] + g11*sq(u[1]) + g22*sq(u[2]); // velocities
        
        if (fabs(dist2)>0.1) {
            disk_p->tedad++;
        }
        
    
        // f.f
        dist5 = g00*sq(y[0][6+0])+ 2.*g03*y[0][6+0]*y[0][6+3]+
        g11* sq(y[0][6+1]) + g22 * sq(y[0][6+2]) + g33*sq(y[0][6+3]);
                        
        double ratio=y[0][1]/disk_p->ISCO;
        dh=correct*ratio/5.;

        if (y[0][1]<10.) dh/=10.;

        if (dh>200) dh=200.;
        if ((y[0][1]>9000.)&&(dh>100)) dh=100.;
        
        for (int i=2;i<4;i++)
        {
            double Change=dh*u[i];
            if (Change>0.01) dh=dh*0.01/Change;
        }

        if (dh>fabs(den)/10.) dh=fabs(den)/10.;
        
        if (dh<0.00001) dh=0.00001;
  
        flagtest=0;
        
        double tp=theta_prev;
        double tc=y[0][2];
        
        if (tp<0.)
        {
            tp=-tp;
            tc=-tc;
        }
        
        if (tp>PI)
        {
            tp=tp-(int)(tp/PI)*PI;
            tc=tc-(int)(tc/PI)*PI;
        }
        


        /*________________________________________________________________________________________________________________*/
        
        // check if the photon hits the disk from above and bellow
        if (((y[0][1]>disk_p->rCorona_ends)&&(y[0][1]<disk_p->r2)&&(timeout>1))&&(((tp<PI/2.)&&(tc>PI/2.))||((tp>PI/2.)&&(tc<PI/2.)))) // The disk starts after ISCO
        {
            // scatter photon!
            
            
            //get the Energy before scatering in disk in BL frame
            if(gamma_p->nScatter>(MaxscatD-1))
            {
                cout << " change size of the scatter array " <<gamma_p->nScatter <<endl;
                exit(1);
            }
            
            
            double eimu[4][4];
            double emui[4][4];
            
            PF_eimu(disk_p->M,disk_p->a,disk_p->epsilon3,y[0],eimu,emui);  //transfrom BL to PF and reverse
            
            // first transform u and f into the PF:
            int ir=logBin(disk_p->r1,disk_p->r2,disk_p->n_disk_structure,y[0][1]);
            // get eimu!

#ifdef DEBUG_SCAT
            cout <<endl<<"Scattering "<<"y[0][1-3] "<<y[0][1]<<" "<<y[0][2]<<" "<<y[0][3]<<" "<<" ir "<<ir<<endl;
#endif
            
            for (int i=0;i<4;i++){
                uPF[i]=0.;
                fPF[i]=0;
                for (int mu=0; mu<4; mu++) {
                    uPF[i] += eimu[i][mu]*u[mu];
                    fPF[i] += eimu[i][mu]*y[0][6+mu];   
                }
            }
            
            
            if (uPF[0]<0) {
                if(gamma_p->E>0){
                    disk_p->nn+=1;
                    cout<<" photons number "<< gamma_p->totalN<<endl;
                }
                else{
                    uPF[0]*=-1;
                    uPF[1]*=-1;
                    uPF[2]*=-1;
                    uPF[3]*=-1;
                }
                
            }
          

            
            // Scatter photon back into the upper hemisphere:
            
            // get PF k vector of scattered photon
            double amu=ran1(&seed);
            double phi=2. * PI * ran1(&seed);
            
            gamma_p->nScatter++;
            
            // the photon is emitted into the upper hemisphere; phi starts at 0 at e_r, and is positive towards e_phi.  
            vPF[0] =  uPF[0] ; 
            vPF[1] =  -uPF[0]* sqrt(1.-sq(amu))*cos(phi); // r ; phi=0 ==> the phi of the direction is 0. 
            vPF[2] =  -uPF[0]* (amu); // theta ==> emission into the upper hemisphere!
            vPF[3] =  -uPF[0]* sqrt(1.-sq(amu))*sin(phi); // phi

            double nvPF[3];
            normalize(&vPF[1],nvPF);

            // get parameters for Chandrasekhar's scattering equation (p. 260)   
            
            double kV[3];
            normalize(&uPF[1],kV);
            
            double ethetaN[3]={0.,-1.,0.}; 
            double etheta_perp[3];
            perp(ethetaN,kV,etheta_perp); 
            
            double f[3];
            normalize(&fPF[1],f);
            
            double f_perp[3];
            perp(f,kV,f_perp);

            double help[3];
            cross(ethetaN,kV,help);
            int sign=(dot(help,f_perp)>0);
            
            // this gives only positive u's!
            double chi= angle(etheta_perp,f_perp);
            
            double mu0=fabs(kV[1]);
            
            double kVn[3];
            kVn[0]=-kV[0];
            kVn[1]=-kV[1];
            kVn[2]=-kV[2];
            
            double phi0;
            phi0 = atan2(kVn[2],kVn[0]); // phi=0 for kVN=e_r, and PI/2. for kVN=e_phi
            if (phi0<0.) phi0+=2.*PI;
            
            double stokes_i = 1.;
            double stokes_q = gamma_p->pol * cos(2.*chi);
            double stokes_u = -sign*gamma_p->pol * sin(2.*chi); // the minus is required because Chandra defined U seen from the 'far end' of the incoming rays
            double pol_in=gamma_p->pol;
            
            double chandra[3];
            
            chandra[0] = (stokes_i+stokes_q)/2.;
            chandra[1] = (stokes_i-stokes_q)/2.;
            chandra[2] = stokes_u;
                        
            double I[3];
            tab25(chandra,mu0,phi0,I,amu,phi); 

            stokes_i=I[0]+I[1];
            stokes_q=I[0]-I[1];
            stokes_u=I[2];
            
            gamma_p->pol = sqrt(sq(stokes_q)+sq(stokes_u))/stokes_i;
            gamma_p->weight *= stokes_i;

            // use I to get "f after scattering" 
            chi=getChi(stokes_q,stokes_u);
            double e_l[3],e_r[3];
            
            perp(ethetaN,nvPF,e_l);
            cross(e_l,nvPF,e_r);
            
            fPF[0]=0.;
            for(int i=0;i<3;i++)
                fPF[i+1]=cos(chi)*e_l[i]+sin(chi)*e_r[i];
  
            double check=0.; 
            for (int i=0; i<4; i++)
                check += eta[i]*sq(vPF[i]);

//            cout <<"vPF * vPF "<<check<<endl;
            // transform everything back into BL coordinates and fill y[0][4] ... y[0][9] 
            for (int mu=0;mu<4;mu++){
                u[mu]=0.;
                for (int i=0; i<4; i++) {
                    u[mu] += emui[mu][i]*vPF[i]*abs(gamma_p->E);
                }                
            }

            // now re-stuff y[0][4] and y[0][5]
            // note that E, L, and b need to be re-calculated!
            
            gamma_p->E = -1.*(g00*u[0] + g03 * u[3]);
            gamma_p->L =     (g33*u[3] + g03 * u[0]);
            gamma_p->b = gamma_p->L/gamma_p->E;
            
            y[0][4]=u[1]/gamma_p->E;
            y[0][5]=u[2]/gamma_p->E;
            
#ifdef DEBUG_SCAT
            cout <<"vPF "<<vPF[0]<<" "<<vPF[1]<<" "<<vPF[2]<<" "<<vPF[3]<<endl;
            cout <<"uBL "<<  u[0]<<" "<<  u[1]<<" "<<  u[2]<<" "<<  u[3]<<endl;
            cout <<"E   "<<gamma_p->E<<" L "<<gamma_p->L<<" b "<<gamma_p->b<<endl;
#endif 
            
            for (int i=0;i<4;i++){
                gamma_p->xSC[i] = y[0][i];
                gamma_p->uSC[i] = u[i];
            }
            
            for (int mu=0;mu<4;mu++){
                y[0][6+mu]=0.;
                for (int i=0; i<4; i++) y[0][6+mu] += emui[mu][i]*fPF[i];
            }
            
            flagtest=1;
            
            if (fabs(u[0])>1) {
                dh/=fabs(u[0]);
                
            }
            
            y[0][2]=PI/2;
            
        };
        
        /*____________________________________________________________________________________________________________________*/
        
        // check if photon is in corona
        
        // flagtest here make sure that if the photon scattered off the disk, do not scatter in in the same step of the corona.
        if ((y[0][1]>=disk_p->rCorona_start)&&(y[0][1]<=disk_p->rCorona_ends)&&(timeout_corona<2000)&&(flagtest==0)) //Here I am assuming that corona can have the Max scat of timeout_corona
        {
            
            // get transformation elements for BL-ZAMO
            double emunu[4][4];
            double enumu[4][4];
            
            BZ_emunu(disk_p->M,disk_p->a,disk_p->epsilon3,y[0],emunu,enumu);   //emunu transforms BL to Zamo and enumu transforms back
            
            //get bin info (DISK INFO NOT CORONA NEED TO IMPROVE)
            
            //calculate BL displacement 4-vector
            double dispBL[4];
            dispBL[0] = y[0][0]-t_prev;
            dispBL[1] = y[0][1]-r_prev;
            dispBL[2] = y[0][2]-theta_prev;
            dispBL[3] = y[0][3]-phi_prev;
            
            double dispZAMO[4];
            for (int i=0;i<4;i++)
            {
                dispZAMO[i]=0;
                for (int mu=0; mu<4; mu++)
                {
                    dispZAMO[i] += enumu[i][mu]*dispBL[mu];
                }
            }
            
            
            //find spacial magnitude dl
            
            double dl = sqrt(sq(dispZAMO[1])+sq(dispZAMO[2])+sq(dispZAMO[3]));
            
            
            //calculate density and probability of scattering (Schnittman&Krolik 2010, pg.910)
            
            double Tcorona= 100; //Corona temp in kev
            double tau0 = 1.5; //vertical optical depth
            double alpha=tau0/(disk_p->rCorona_ends-disk_p->rCorona_start);
            if (alpha<0) {
                cout<<"Wrong slope for tau"<<endl;
                exit(1);
            }
            
            double dtau=dl*alpha;
            double Pscat = 1.-exp(-dtau);
            
            if(dtau>0.2)
                cout<<" dtae<< "<<dtau<<endl;
            //cout<<">>> Pscat = "<<Pscat<<endl;
            
            //generate random number
            double p = ran1(&seed);
            
            double f_corona[4];
            for (int i=0; i<4; i++) {
                f_corona[i]=y[0][6+i];
            }
            
            // check if photon scatters!!!
            if (p<Pscat)
            {
                
                Scatter(gamma_p, u,f_corona, Tcorona, emunu,enumu,p/Pscat);
                
                
                gamma_p->E = -1.*(g00*u[0] + g03 * u[3]);
                gamma_p->L =     (g33*u[3] + g03 * u[0]);
                gamma_p->b = gamma_p->L/gamma_p->E;
                
                y[0][4]=u[1]/gamma_p->E;
                y[0][5]=u[2]/gamma_p->E;
                
                
                for (int i=0; i<4; i++)
                    y[0][6+i] = f_corona[i];
                
                
                if (fabs(u[0])>1) {
                    dh/=fabs(u[0]);
                    
                }

                
            }
            
            timeout_corona++;
            
            gamma_p->timeout_corona++;
        };
        
        /*________________________________________________________________________________________________________________*/
        
        
        t_prev=y[0][0];
        r_prev=y[0][1];
        theta_prev=y[0][2];
        phi_prev=y[0][3];
        y04_prev=y[0][4];
        y05_prev=y[0][5];
        f_prev[0]=y[0][6];
        f_prev[1]=y[0][7];
        f_prev[2]=y[0][8];
        f_prev[3]=y[0][9];
        adjustment=1;


        if (gamma_p->timeout==1000001)
            timeout=1000001;
        
        // check if we are hitting the horizon.
        double mu    = cos(y[0][2]);
        int    index = DISK_HORIZON_MU-linBin(-1.,1.,DISK_HORIZON_MU,mu)-1.;
        rh           = disk_p->horizon[index];
        
#ifdef PRINT 
        // pipe results to file
        if((flag>=0)&&((timeout<10000)||(timeout%20==0)))
        {
           os <<flag<<" "<<y[0][0]<<" "<<y[0][1]<<" "<<y[0][2]<<" "<<y[0][3]<<" "<<
            u[0]<<" "<<u[1]<<" "<<u[2]<<" "<<u[3]<<" "<<endl;
        };
#endif
        
    } while ( (fabs(dist2)<.5) &&  (y[0][1]>1.02*rh) && (y[0][1]<disk_p->rmax) && (timeout++ < 100000));
        
    gamma_p->xBL[0]=y[0][0]; // t
    gamma_p->xBL[1]=y[0][1]; // r
    gamma_p->xBL[2]=y[0][2]; // theta
    gamma_p->xBL[3]=y[0][3]; // phi
    
    for (int i=0;i<4;i++) u[i] *= gamma_p->E;
    
    gamma_p->uBL[0] = u[0]; // d_theta/d_lambda
    gamma_p->uBL[1] = u[1]; // d_r/d_lambda
    gamma_p->uBL[2] = u[2]; // d_theta/d_lambda
    gamma_p->uBL[3] = u[3]; // d_phi/d_lambda
    
    gamma_p->fBL[0]=y[0][6];
    gamma_p->fBL[1]=y[0][7];
    gamma_p->fBL[2]=y[0][8];
    gamma_p->fBL[3]=y[0][9];
    
    if(y[0][1]>=disk_p->rmax)
    {
        // transform into coordinate-stationary frame.
        double emui[4][4];
        double eimu[4][4];
        
        CS_emui(disk_p->M,disk_p->a,disk_p->epsilon3,y[0],emui);
        migs(emui,4,eimu);
        
        double uCS[4];
        double fCS[4];
        
        for (int i=0;i<4;i++){
            uCS[i]=0.;
            fCS[i]=0;
            for (int mu=0; mu<4; mu++) {
                uCS[i] += eimu[i][mu]*u[mu];   
                fCS[i] += eimu[i][mu]*y[0][6+mu];   
            }
            gamma_p->uCS[i] = uCS[i];
            gamma_p->fCS[i] = fCS[i];
        }
        
#ifdef DEBUG_SCAT
        if (gamma_p->nScatter>0)
            cout <<"uCS "<<uCS[0]<<" "<<uCS[1]<<" "<<uCS[2]<<" "<<uCS[3]<<endl<<endl;
#endif        
        
        g00 = metric(0,0,disk_p->M,disk_p->a,disk_p->epsilon3,y[0]);
        g03 = metric(0,3,disk_p->M,disk_p->a,disk_p->epsilon3,y[0]);
        g11 = metric(1,1,disk_p->M,disk_p->a,disk_p->epsilon3,y[0]);
        g22 = metric(2,2,disk_p->M,disk_p->a,disk_p->epsilon3,y[0]);
        g33 = metric(3,3,disk_p->M,disk_p->a,disk_p->epsilon3,y[0]);

        dist2 = g00*sq(u[0]) + g33*sq(u[3]) + 2.*g03*u[0]*u[3] + g11*sq(u[1]) + g22*sq(u[2]); // velocities
        
        gamma_p->NuCS=0.; 
        for (int i=0; i<4; i++)
            gamma_p->NuCS += eta[i]*sq(uCS[i]);
        
        gamma_p->NfCS=0.;
        for (int i=0; i<4; i++)
            gamma_p->NfCS += eta[i]*sq(fCS[i]);
        
        // get Stokes parameters in this frame.
        
        double kVV[3];
        normalize(&uCS[1],kVV);
        
        double ethetaN[3]={0.,-1.,0.}; 
        double etheta_perp[3];
        perp(ethetaN,kVV,etheta_perp);
        
        double f[3];
        normalize(&fCS[1],f);
        
        double f_perp[3];
        perp(f,kVV,f_perp);
        
        double help[3];
        cross(ethetaN,kVV,help);
        double sign=dot(help,f_perp);
        double costest=dot(etheta_perp,f_perp);
        double chi= atan2(sign,costest);
        if(chi<0)
            chi+=PI;
        
        gamma_p->stokes[0] = gamma_p->weight;
        gamma_p->stokes[1] = gamma_p->weight * gamma_p->pol * cos(2.*chi);
        gamma_p->stokes[2] = gamma_p->weight * gamma_p->pol * sin(2.*chi);
        
        double testchi=getChi(gamma_p->stokes[1],gamma_p->stokes[2]);
        if (abs(testchi-chi)>0.01) {
            cout<<" chi me "<<chi<<" chitest "<<testchi<<endl;

        }
        
    }
    else {
        for (int i=0; i<4; i++){
            gamma_p->uCS[i] = 0.;
            gamma_p->fCS[i] = 0.;
        }
        gamma_p->NuCS=dist2; 
        gamma_p->NfCS=dist5;
        gamma_p->stokes[0] = gamma_p->weight;
        gamma_p->stokes[1] = 0.;
        gamma_p->stokes[2] = 0.;
    }

    gamma_p->timeout=timeout;        
}
/* Physics functions */

void tab24(double mu, double &Is, double &pol)
/* from Chandrasekhar "Radiative Transfer", p. 248 */
{
    double res;
#define Ntab24 21
    
    static double t_mu [Ntab24] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};
    static double Il_v [Ntab24] = {0.18294,0.21613,0.24247,0.26702,0.29057,0.31350,0.33599,0.35817,0.38010,0.40184,
        0.42343,0.44489,0.46624,0.48750,0.50869,0.52981,0.55087,0.57189,0.59286,0.61379,0.63469};
    static double Ir_v [Ntab24] = {0.23147,0.25877,0.28150,0.30299,0.32381,0.34420,0.36429,0.38417,0.40388,0.42346,
        0.44294,0.46233,0.48165,0.50092,0.52013,0.53930,0.55844,0.57754,0.59661,0.61566,0.63469};
    
    if (mu<0.) mu=0.;
    if (mu>1.) mu=1.;
    
    double Ir  = get_value_lin(t_mu,Ir_v,Ntab24,mu);
    double Il  = get_value_lin(t_mu,Il_v,Ntab24,mu);
    Is  = mu*(Ir+Il);
    pol = (Ir-Il)/(Ir+Il);
}

void tab25(double *F, double mu0, double phi0, double *I, double mu, double phi)
/* from Chandrasekhar "Radiative Transfer", p. 259 */
// be careful:
// the number of photons is propotional to the energy flux (in PF)
// we want to multiply a photon's weight with a factor which averages out to 1 when we average over all scattering processes (actually we could check that we get this right)
//
// important:
// let's call incoming flux "F"; Chandrasekhar's "F" is F/PI;
// if we want to have unit energy flux per unit area of the disk coming in, we have to feed the program with F/mu.
// outgoing is I, we have to multiply with mu to get the outgoing <flux per area>. 
//
// STILL TO DO: CHECK THET U ACQUIRES THE RIGHT SIGN.
//
{
#define Ntab25 21
    static double t_mu [Ntab25] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};
    
    static double Tpsi  [Ntab25] = {0.00000,0.04075,0.09144,0.15105,0.21916,0.29557,0.38014,0.47276,0.57338,0.68195,0.79843,0.92277,1.05497,1.19501,1.34286,1.49852,1.66198,1.83321,2.01223,2.19902,2.39357};
    static double Tphi  [Ntab25] = {1.00000,1.12988,1.20976,1.26850,1.31108,1.33973,1.35569,1.35971,1.35228,1.33374,1.30437,1.26432,1.21375,1.15279,1.08153,1.00003,0.90836,0.80655,0.69468,0.57276,0.44083};
    static double Tchi  [Ntab25] = {1.00000,1.10352,1.18638,1.26329,1.33687,1.40821,1.47801,1.54664,1.61435,1.68132,1.74772,1.81362,1.87911,1.94425,2.00907,2.07365,2.13799,2.20213,2.26609,2.32990,2.39356};
    static double Tceta [Ntab25] = {0.00000,0.01824,0.03764,0.05780,0.07852,0.09969,0.12121,0.14303,0.16510,0.18738,0.20984,0.23247,0.25523,0.27812,0.30112,0.32421,0.34739,0.37065,0.39398,0.41738,0.44083};
    static double Th1mu [Ntab25] = {1.00000,1.07167,1.11602,1.14837,1.17155,1.18685,1.19487,1.19599,1.19030,1.17774,1.15816,1.13118,1.09624,1.05256,0.99899,0.93381,0.85435,0.75611,0.63033,0.45471,0.00000};
    static double Th2   [Ntab25] = {1.00000,1.04967,1.08621,1.11762,1.14552,1.17075,1.19383,1.21508,1.23476,1.25308,1.27019,1.28624,1.30132,1.31554,1.32895,1.34166,1.35371,1.36515,1.37601,1.38638,1.39625};
    
    mu= fabs(mu);
    if (mu<0.) mu=0.;
    if (mu>1.) mu=1.;
    
    double res;
    double S[3][3];
    double m1[3][3],m2[3][3],m3[3][3],m4[3][3];
    double Q[3][3]={1.,0.,0.,0.,1.,0.,0.,0.,2.};
    
    double psiMu   = get_value_lin(t_mu,Tpsi,Ntab25,mu);
    double psiMu0  = get_value_lin(t_mu,Tpsi,Ntab25,mu0);
    
    double phiMu   = get_value_lin(t_mu,Tphi,Ntab25,mu);
    double phiMu0  = get_value_lin(t_mu,Tphi,Ntab25,mu0);
    
    double chiMu   = get_value_lin(t_mu,Tchi,Ntab25,mu);
    double chiMu0  = get_value_lin(t_mu,Tchi,Ntab25,mu0);
    
    double cetaMu  = get_value_lin(t_mu,Tceta,Ntab25,mu);
    double cetaMu0 = get_value_lin(t_mu,Tceta,Ntab25,mu0);
    
    double h1Mu    = get_value_lin(t_mu,Th1mu,Ntab25,mu);
    double h1Mu0   = get_value_lin(t_mu,Th1mu,Ntab25,mu0);
    
    double h2Mu    = get_value_lin(t_mu,Th2,Ntab25,mu);
    double h2Mu0   = get_value_lin(t_mu,Th2,Ntab25,mu0);
    
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            S[i][j]=m1[i][j]=m2[i][j]=m3[i][j]=m4[i][j]=0.;
    
    double st=sqrt(2.);
    
    m1[0][0]=0.75 * psiMu;
    m1[0][1]=0.75 * st * phiMu; 
    m1[1][0]=0.75 * chiMu;
    m1[1][1]=0.75 * st * cetaMu;
    
    m2[0][0]=     psiMu0;
    m2[0][1]=     chiMu0; 
    m2[1][0]=st * phiMu0;
    m2[1][1]=st * cetaMu0;
    
     
    Matrix_Mult(m1,m2,m1);
    
    
    m3[0][0]= -4.* mu * mu0 * cos(phi0-phi); 
    m3[0][2]=  2.* mu       * sin(phi0-phi);
    m3[2][0]=  2.* mu0      * sin(phi0-phi); 
    m3[2][2]=                 cos(phi0-phi);
    
    double fac = 0.75 * h1Mu * h1Mu0;
    
    Matrix_SMult(fac,m3,m3);
    
    Matrix_Add(m1,m3,S);
    
    m4[0][0]= sq(mu)*sq(mu0)*cos(2*(phi0-phi));
    m4[0][1]=-sq(mu)*cos(2*(phi0-phi));
    m4[0][2]=-sq(mu)*mu0*sin(2*(phi0-phi));
    
    m4[1][0]=-sq(mu0)*cos(2*(phi0-phi));
    m4[1][1]= cos(2*(phi0-phi));
    m4[1][2]= mu0*sin(2*(phi0-phi));
    
    m4[2][0]=-mu*sq(mu0)*sin(2*(phi0-phi));
    m4[2][1]= mu*sin(2*(phi0-phi));
    m4[2][2]=-mu*mu0*cos(2*(phi0-phi));
    
    fac = 0.75 * h2Mu * h2Mu0;
    Matrix_SMult(fac,m4,m4);
    
    
    Matrix_Add(S,m4,S);
    
    fac = 1./(1./mu0+1./mu);
    Matrix_SMult(fac,S,S);
    
    Matrix_Mult(Q,S,S);
    Matrix_SMult(2.*PI*mu/(4.*mu*mu0*PI),S,S); 
    // we multiply with 2 PI because: Il+Ir integrated over upper hemisphere gives 1. 
    // if it were constant, then we had Il+Ir = 1/ 2 PI everywhere; we would have to multiply with 2 PI to get a weight normalized to 1.
    // We multiply with mu to get the energy flux per area of the accretion disk. we divide by PI because we use F (and not PI x F) as the incoming flux. 
    
    
    MatrixVector_Mult(S,F,I);
    
}

double ChristoffelKerr1(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return -(M*((p.y1_2 + p.a_2)*p.pow2_y2_2_p_a_2cos_y2_2*(-2*p.y1_6 + p.y1_3*hair*(3*p.y1 - 8*M)*p.M_2 + 2*p.y1*p.a_2*(-p.y1_3 + hair*p.M_2*(p.y1 + 2*M))*p.cos_y2_2 + p.a_4*(2*p.y1_2 - hair*p.M_2)*p.cos_y2_4 + 2*p.a_6*p.cos_y2_6) + p.y1*p.a_2*hair*p.M_3*(2*p.y1*p.a_2*hair*p.M_2*(p.y1 + M)*p.cos_y2_2 + p.a_4*(4*p.y1_2 - hair*p.M_2)*p.cos_y2_4 + p.y1_3*(4*p.y1*(p.y1_2 + p.a_2) + hair*(3*p.y1 - 2*M)*p.M_2 + 4*p.y1*p.a_2*p.cos_2y2))*p.sin_y2_2))/(2.*p.pow2_y2_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4)*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr2(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (p.y1*p.a_2*M*p.cos_y2*p.sin_y2*(2*p.a_6*(-4*p.y1_2*(p.a_2 + p.y1*(p.y1 - 2*M)) + (p.y1_2 + p.a_2)*hair*p.M_2)*p.cos_y2_6 - 2*p.a_8*(p.a_2 + p.y1*(p.y1 - 2*M))*p.cos_y2_8 - p.y1*p.a_4*p.cos_y2_4*(12*p.y1_3*(p.y1_2 + p.a_2 - 2*p.y1*M) + hair*p.M_2*(-6*p.y1_3 - 6*p.y1*p.a_2 + 6*p.y1_2*M + 5*p.a_2*M - 8*p.y1*p.M_2) + p.a_2*hair*p.M_3*p.cos_2y2) + 2*p.y1*p.a_2*p.cos_y2_2*(-4*p.y1_5*(p.a_2 + p.y1*(p.y1 - 2*M)) + p.y1_2*hair*p.M_2*(3*p.a_2*(p.y1 - 2*M) + p.y1*(3*p.y1_2 - 6*p.y1*M + 8*p.M_2)) + p.a_2*hair*p.M_3*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2) + 2*p.y1_2*(-(p.y1_6*(p.a_2 + p.y1*(p.y1 - 2*M))) + 2*p.y1*p.hair_2*p.M_7 + p.y1_3*hair*p.M_2*(p.a_2*(p.y1 - 3*M) + p.y1*(p.y1_2 - 3*p.y1*M + 4*p.M_2)) + p.a_2*hair*p.M_3*(p.y1_3 + hair*(p.y1 - M)*p.M_2)*p.sin_y2_2)))/(p.pow2_y1_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4)*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr4(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return -(M*((p.y1_2 + p.a_2)*p.pow2_y2_2_p_a_2cos_y2_2*(-2*p.y1_6 + p.y1_3*hair*(3*p.y1 - 8*M)*p.M_2 + 2*p.y1*p.a_2*(-p.y1_3 + hair*p.M_2*(p.y1 + 2*M))*p.cos_y2_2 + p.a_4*(2*p.y1_2 - hair*p.M_2)*p.cos_y2_4 + 2*p.a_6*p.cos_y2_6) + p.y1*p.a_2*hair*p.M_3*(2*p.y1*p.a_2*hair*p.M_2*(p.y1 + M)*p.cos_y2_2 + p.a_4*(4*p.y1_2 - hair*p.M_2)*p.cos_y2_4 + p.y1_3*(4*p.y1*(p.y1_2 + p.a_2) + hair*(3*p.y1 - 2*M)*p.M_2 + 4*p.y1*p.a_2*p.cos_2y2))*p.sin_y2_2))/(2.*p.pow2_y2_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4)*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr7(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (a*M*p.sin_y2_2*(-(p.y1_7*(3*p.y1_2*(p.y1_3 + 2*hair*p.M_3) + p.a_2*(p.y1_3 + 4*hair*p.M_3))) + p.a_2*(2*p.y1*p.a_4*(-9*p.y1_5 + p.a_2*(p.y1_3 + hair*p.M_3))*p.cos_y2_6 + p.y1_2*p.a_6*(-7*p.y1_2 + 3*p.a_2)*p.cos_y2_8 + p.a_8*(-p.y1 + a)*(p.y1 + a)*p.cos_y2_10 - 2*p.y1_3*p.a_2*p.cos_y2_4*(p.y1_2*(11*p.y1_3 + p.y1*p.a_2 + 3*hair*p.M_3) - p.a_2*hair*p.M_3*p.sin_y2_2) - p.y1_2*p.cos_y2_2*(13*p.y1_8 + 3*p.y1_6*p.a_2 + 6*p.y1_3*(2*p.y1_2 + p.a_2)*hair*p.M_3 - p.a_2*p.hair_2*p.M_6*p.sin_y2_2) + p.y1_4*hair*p.M_3*((2*p.y1_3 - hair*p.M_3)*p.sin_y2_2 + p.y1*p.a_2*p.sin_2y2_2))))/(p.pow2_y1_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4)*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr8(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (p.y1*p.a_2*M*p.cos_y2*p.sin_y2*(2*p.a_6*(-4*p.y1_2*(p.a_2 + p.y1*(p.y1 - 2*M)) + (p.y1_2 + p.a_2)*hair*p.M_2)*p.cos_y2_6 - 2*p.a_8*(p.a_2 + p.y1*(p.y1 - 2*M))*p.cos_y2_8 - p.y1*p.a_4*p.cos_y2_4*(12*p.y1_3*(p.y1_2 + p.a_2 - 2*p.y1*M) + hair*p.M_2*(-6*p.y1_3 - 6*p.y1*p.a_2 + 6*p.y1_2*M + 5*p.a_2*M - 8*p.y1*p.M_2) + p.a_2*hair*p.M_3*p.cos_2y2) + 2*p.y1*p.a_2*p.cos_y2_2*(-4*p.y1_5*(p.a_2 + p.y1*(p.y1 - 2*M)) + p.y1_2*hair*p.M_2*(3*p.a_2*(p.y1 - 2*M) + p.y1*(3*p.y1_2 - 6*p.y1*M + 8*p.M_2)) + p.a_2*hair*p.M_3*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2) + 2*p.y1_2*(-(p.y1_6*(p.a_2 + p.y1*(p.y1 - 2*M))) + 2*p.y1*p.hair_2*p.M_7 + p.y1_3*hair*p.M_2*(p.a_2*(p.y1 - 3*M) + p.y1*(p.y1_2 - 3*p.y1*M + 4*p.M_2)) + p.a_2*hair*p.M_3*(p.y1_3 + hair*(p.y1 - M)*p.M_2)*p.sin_y2_2)))/(p.pow2_y2_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4)*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr11(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (p.y1*p.a_3*M*p.cos_y2*p.sin_y2_3*(2*p.y1*p.a_6*(4*p.y1*(p.a_2 + p.y1*(p.y1 - 2*M)) - hair*p.M_3)*p.cos_y2_6 + 2*p.a_8*(p.a_2 + p.y1*(p.y1 - 2*M))*p.cos_y2_8 + p.y1*p.a_4*p.cos_y2_4*(hair*p.M_3*(5*p.a_2 - 8*p.y1*M) + 12*p.y1_3*(p.y1_2 + p.a_2 - 2*p.y1*M) + p.a_2*hair*p.M_3*p.cos_2y2) + 2*p.y1_2*p.a_2*p.cos_y2_2*(-((p.y1_3 + hair*p.M_3)*(-4*p.y1_3 + 8*p.y1_2*M + hair*p.M_3)) + p.a_2*(4*p.y1_4 + 5*p.y1*hair*p.M_3) + p.y1*p.a_2*hair*p.M_3*p.cos_2y2) + 2*p.y1_2*(p.y1_3*hair*(3*p.a_2 + 2*p.y1*(p.y1 - 2*M))*p.M_3 - p.y1*p.hair_2*p.M_6*(p.y1 + 2*M) + p.y1_6*(p.y1_2 + p.a_2 - 2*p.y1*M) + p.a_2*hair*p.M_3*(-p.y1_3 + hair*p.M_3)*p.sin_y2_2)))/(p.pow2_y1_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4)*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr13(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (a*M*p.sin_y2_2*(-(p.y1_7*(3*p.y1_2*(p.y1_3 + 2*hair*p.M_3) + p.a_2*(p.y1_3 + 4*hair*p.M_3))) + p.a_2*(2*p.y1*p.a_4*(-9*p.y1_5 + p.a_2*(p.y1_3 + hair*p.M_3))*p.cos_y2_6 + p.y1_2*p.a_6*(-7*p.y1_2 + 3*p.a_2)*p.cos_y2_8 + p.a_8*(-p.y1 + a)*(p.y1 + a)*p.cos_y2_10 - 2*p.y1_3*p.a_2*p.cos_y2_4*(p.y1_2*(11*p.y1_3 + p.y1*p.a_2 + 3*hair*p.M_3) - p.a_2*hair*p.M_3*p.sin_y2_2) - p.y1_2*p.cos_y2_2*(13*p.y1_8 + 3*p.y1_6*p.a_2 + 6*p.y1_3*(2*p.y1_2 + p.a_2)*hair*p.M_3 - p.a_2*p.hair_2*p.M_6*p.sin_y2_2) + p.y1_4*hair*p.M_3*((2*p.y1_3 - hair*p.M_3)*p.sin_y2_2 + p.y1*p.a_2*p.sin_2y2_2))))/(p.pow2_y1_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4)*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr14(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (p.y1*p.a_3*M*p.cos_y2*p.sin_y2_3*(2*p.y1*p.a_6*(4*p.y1*(p.a_2 + p.y1*(p.y1 - 2*M)) - hair*p.M_3)*p.cos_y2_6 + 2*p.a_8*(p.a_2 + p.y1*(p.y1 - 2*M))*p.cos_y2_8 + p.y1*p.a_4*p.cos_y2_4*(hair*p.M_3*(5*p.a_2 - 8*p.y1*M) + 12*p.y1_3*(p.y1_2 + p.a_2 - 2*p.y1*M) + p.a_2*hair*p.M_3*p.cos_2y2) + 2*p.y1_2*p.a_2*p.cos_y2_2*(-((p.y1_3 + hair*p.M_3)*(-4*p.y1_3 + 8*p.y1_2*M + hair*p.M_3)) + p.a_2*(4*p.y1_4 + 5*p.y1*hair*p.M_3) + p.y1*p.a_2*hair*p.M_3*p.cos_2y2) + 2*p.y1_2*(p.y1_3*hair*(3*p.a_2 + 2*p.y1*(p.y1 - 2*M))*p.M_3 - p.y1*p.hair_2*p.M_6*(p.y1 + 2*M) + p.y1_6*(p.y1_2 + p.a_2 - 2*p.y1*M) + p.a_2*hair*p.M_3*(-p.y1_3 + hair*p.M_3)*p.sin_y2_2)))/(p.pow2_y1_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4)*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr16(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (M*(p.y1_3*(2*p.y1_3 - 3*p.y1*hair*p.M_2 + 8*hair*p.M_3) + 2*p.y1*p.a_2*(p.y1_3 - hair*p.M_2*(p.y1 + 2*M))*p.cos_y2_2 + p.a_4*(-2*p.y1_2 + hair*p.M_2)*p.cos_y2_4 - 2*p.a_6*p.cos_y2_6)*((p.a_2 + p.y1*(p.y1 - 2*M))*p.pow2_y1_2_p_a_2cos_y2_2 + p.y1*p.a_2*hair*p.M_3*p.sin_y2_2))/(2.*p.pow5_y1_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4));
}


double ChristoffelKerr19(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (a*M*(-(p.y1_3*(p.y1_3 + 4*hair*p.M_3)) - p.y1*p.a_2*(p.y1_3 - 2*hair*p.M_3)*p.cos_y2_2 + p.y1_2*p.a_4*p.cos_y2_4 + p.a_6*p.cos_y2_6)*p.sin_y2_2*((p.a_2 + p.y1*(p.y1 - 2*M))*p.pow2_y1_2_p_a_2cos_y2_2 + p.y1*p.a_2*hair*p.M_3*p.sin_y2_2))/(p.pow5_y1_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4));
}

double ChristoffelKerr21(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (5*p.y1_2 - p.a_2*p.cos_y2_2 + (p.y1*hair*p.M_3*(-3*p.y1_2 + p.a_2*p.cos_y2_2))/(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4) + (p.pow2_y2_2_p_a_2cos_y2_2*(p.y1_2*(-3*p.a_2 + p.y1*(-5*p.y1 + 8*M)) + p.a_2*(-p.y1 + a)*(p.y1 + a)*p.cos_y2_2))/((p.a_2 + p.y1*(p.y1 - 2*M))*p.pow2_y2_2_p_a_2cos_y2_2 + p.y1*p.a_2*hair*p.M_3*p.sin_y2_2))/(2.*(p.y1_3 + p.y1*p.a_2*p.cos_y2_2));
}


double ChristoffelKerr22(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return p.a_2*p.cos_y2*p.sin_y2*(-(1/(p.y1_2 + p.a_2*p.cos_y2_2)) - (2*(p.y1_2 + p.a_2*p.cos_y2_2))/(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4) + (p.y1*(2*p.y1*(p.a_2 + p.y1*(p.y1 - 2*M)) - hair*p.M_3) + 2*p.a_2*(p.a_2 + p.y1*(p.y1 - 2*M))*p.cos_y2_2)/((p.a_2 + p.y1*(p.y1 - 2*M))*p.pow2_y1_2_p_a_2cos_y2_2 + p.y1*p.a_2*hair*p.M_3*p.sin_y2_2));
}


double ChristoffelKerr25(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return p.a_2*p.cos_y2*p.sin_y2*(-(1/(p.y1_2 + p.a_2*p.cos_y2_2)) - (2*(p.y1_2 + p.a_2*p.cos_y2_2))/(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4) + (p.y1*(2*p.y1*(p.a_2 + p.y1*(p.y1 - 2*M)) - hair*p.M_3) + 2*p.a_2*(p.a_2 + p.y1*(p.y1 - 2*M))*p.cos_y2_2)/((p.a_2 + p.y1*(p.y1 - 2*M))*p.pow2_y1_2_p_a_2cos_y2_2 + p.y1*p.a_2*hair*p.M_3*p.sin_y2_2));
}


double ChristoffelKerr26(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return -((p.y1*(p.a_2 + p.y1*(p.y1 - 2*M) + (p.y1*p.a_2*hair*p.M_3*p.sin_y2_2)/p.pow2_y1_2_p_a_2cos_y2_2))/(p.a_2*p.cos_y2_2 + p.y1*(p.y1 + (hair*p.M_3)/(p.y1_2 + p.a_2*p.cos_y2_2))));
}


double ChristoffelKerr28(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (a*M*(-(p.y1_3*(p.y1_3 + 4*hair*p.M_3)) - p.y1*p.a_2*(p.y1_3 - 2*hair*p.M_3)*p.cos_y2_2 + p.y1_2*p.a_4*p.cos_y2_4 + p.a_6*p.cos_y2_6)*p.sin_y2_2*((p.a_2 + p.y1*(p.y1 - 2*M))*p.pow2_y1_2_p_a_2cos_y2_2 + p.y1*p.a_2*hair*p.M_3*p.sin_y2_2))/(p.pow5_y1_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4));
}


double ChristoffelKerr31(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return -(p.sin_y2_2*((p.a_2 + p.y1*(p.y1 - 2*M))*p.pow2_y1_2_p_a_2cos_y2_2 + p.y1*p.a_2*hair*p.M_3*p.sin_y2_2)*(2*p.y1_9 + p.a_2*(2*p.y1*p.a_6*p.cos_y2_8 - p.y1_3*M*(2*p.y1_3 + hair*p.M_2*(3*p.y1 + 8*M))*p.sin_y2_2 + 2*p.cos_y2_6*(4*p.y1_3*p.a_4 + p.a_6*M*p.sin_y2_2) + 2*p.y1_2*p.cos_y2_2*(4*p.y1_5 - p.a_2*M*(p.y1_2 + hair*p.M_2)*p.sin_y2_2) + p.a_2*p.cos_y2_4*(12*p.y1_5 + p.a_2*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2) + p.y1*p.a_2*hair*p.M_4*p.sin_2y2_2)))/(2.*p.pow5_y1_2_p_a_2cos_y2_2*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4));
}


double ChristoffelKerr32(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (-2*p.y1*p.a_2*M*p.cos_y2*(p.y1_4 + p.y1*hair*p.M_2*(-p.y1 + 3*M) + p.a_2*(2*p.y1_2 - hair*p.M_2)*p.cos_y2_2 + p.a_4*p.cos_y2_4)*p.sin_y2)/p.pow5_y1_2_p_a_2cos_y2_2;
}


double ChristoffelKerr35(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (p.y1*a*M*(8*p.y1_6 + 16*p.y1_4*p.a_2 + 11*p.y1_2*p.a_4 + 3*p.a_6 + 8*p.y1*(p.y1_2 + 2*p.a_2)*hair*p.M_3 + 4*p.a_2*(2*p.y1_4 + 3*p.y1_2*p.a_2 + p.a_4 - 2*p.y1*hair*p.M_3)*p.cos_2y2 + p.a_4*(p.y1_2 + p.a_2)*p.cos_4y2)*p.sin_2y2)/(8.*p.pow5_y1_2_p_a_2cos_y2_2);
}


double ChristoffelKerr37(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (p.a_2*p.cos_y2*(p.y1_2 + p.a_2 - 2*p.y1*M + (16*p.y1_2*(p.y1_2 + p.a_2)*p.hair_2*p.M_6)/pow(2*p.y1_2 + p.a_2 + p.a_2*p.cos_2y2,4) + (8*p.y1*hair*p.M_3*(p.a_2 + p.y1*(p.y1 + M)))/pow(2*p.y1_2 + p.a_2 + p.a_2*p.cos_2y2,2) - (4*p.y1*hair*p.M_3)/(2*p.y1_2 + p.a_2 + p.a_2*p.cos_2y2))*p.sin_y2)/((p.y1_2 + p.a_2*p.cos_y2_2)*sqr(p.a_2 + p.y1*(p.y1 - 2*M) + (p.y1*p.a_2*hair*p.M_3*p.sin_y2_2)/p.pow2_y1_2_p_a_2cos_y2_2));
}


double ChristoffelKerr38(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return p.y1/(p.y1_2 + p.a_2*p.cos_y2_2);
}


double ChristoffelKerr41(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return p.y1/(p.y1_2 + p.a_2*p.cos_y2_2);
}


double ChristoffelKerr42(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return -((p.a_2*p.cos_y2*p.sin_y2)/(p.y1_2 + p.a_2*p.cos_y2_2));
}


double ChristoffelKerr44(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (p.y1*a*M*(8*p.y1_6 + 16*p.y1_4*p.a_2 + 11*p.y1_2*p.a_4 + 3*p.a_6 + 8*p.y1*(p.y1_2 + 2*p.a_2)*hair*p.M_3 + 4*p.a_2*(2*p.y1_4 + 3*p.y1_2*p.a_2 + p.a_4 - 2*p.y1*hair*p.M_3)*p.cos_2y2 + p.a_4*(p.y1_2 + p.a_2)*p.cos_4y2)*p.sin_2y2)/(8.*p.pow5_y1_2_p_a_2cos_y2_2);
}


double ChristoffelKerr47(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return -((p.cos_y2*p.sin_y2*((p.y1_2 + p.a_2)*p.pow4_y1_2_p_a_2cos_y2_2 + p.y1*p.a_2*M*(4*p.y1_6 + 7*p.y1_4*p.a_2 + p.y1*hair*p.M_2*(2*p.y1_2*(p.y1 + 2*M) + p.a_2*(3*p.y1 + 5*M)) + 2*p.a_4*(6*p.y1_2 + hair*p.M_2)*p.cos_y2_4 + 4*p.a_6*p.cos_y2_6 + p.y1*p.a_2*(5*p.y1_3 + hair*(p.y1 - M)*p.M_2)*p.cos_2y2)*p.sin_y2_2 + 2*p.y1*p.a_6*M*p.cos_y2_2*(2*p.y1_2 + hair*p.M_2 + p.a_2*p.cos_y2_2)*p.sin_y2_4))/p.pow5_y1_2_p_a_2cos_y2_2);
}


double ChristoffelKerr49(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (a*M*(p.y1_2 - p.a_2*p.cos_y2_2)*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4))/(p.pow2_y1_2_p_a_2cos_y2_2*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr50(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return -(p.y1*a*(p.a_2 + p.y1*(p.y1 - 2*M))*M*(2*p.y1_2 + p.a_2 + p.a_2*p.cos_2y2)*(8*p.y1_4 + 8*p.y1_2*p.a_2 + 3*p.a_4 + 8*p.y1*hair*p.M_3 + 4*p.a_2*(2*p.y1_2 + p.a_2)*p.cos_2y2 + p.a_4*p.cos_4y2)*cot(p.y2))/(8.*p.pow3_y1_2_p_a_2cos_y2_2*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr52(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (a*M*(p.y1_2 - p.a_2*p.cos_y2_2)*(p.y1_4 + p.y1*hair*p.M_3 + 2*p.y1_2*p.a_2*p.cos_y2_2 + p.a_4*p.cos_y2_4))/(p.pow2_y1_2_p_a_2cos_y2_2*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr55(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (2*p.y1*p.pow3_y1_2_p_a_2cos_y2_2*(p.y1*(p.y1 - 2*M) + p.a_2*p.cos_y2_2) + p.a_2*M*(-2*p.y1_6 - p.y1_3*hair*p.M_2*(3*p.y1 + 2*M) - 2*p.y1*p.a_2*(p.y1_3 + hair*(p.y1 - M)*p.M_2)*p.cos_y2_2 + p.a_4*(2*p.y1_2 + hair*p.M_2)*p.cos_y2_4 + 2*p.a_6*p.cos_y2_6)*p.sin_y2_2)/(2.*p.pow2_y1_2_p_a_2cos_y2_2*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr56(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return -(p.y1*a*(p.a_2 + p.y1*(p.y1 - 2*M))*M*(2*p.y1_2 + p.a_2 + p.a_2*p.cos_2y2)*(8*p.y1_4 + 8*p.y1_2*p.a_2 + 3*p.a_4 + 8*p.y1*hair*p.M_3 + 4*p.a_2*(2*p.y1_2 + p.a_2)*p.cos_2y2 + p.a_4*p.cos_4y2)*cot(p.y2))/(8.*p.pow3_y1_2_p_a_2cos_y2_2*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr59(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (((p.y1_2 + p.a_2)*pow(2*p.y1_2 + p.a_2 + p.a_2*p.cos_2y2,3)*(p.a_2 + 2*p.y1*(p.y1 - 2*M) + p.a_2*p.cos_2y2)*cot(p.y2))/16. + p.y1*p.a_2*M*((p.a_2*p.cos_y2_3*(7*p.a_4 + 16*p.y1_2*(3*p.y1_2 - 2*p.y1*M + hair*p.M_2) + 8*p.a_2*(4*p.y1_2 - p.y1*M + hair*p.M_2) + 8*p.a_2*(p.a_2 + p.y1*(2*p.y1 - M))*p.cos_2y2 + p.a_4*p.cos_4y2)*p.sin_y2)/4. + p.y1*(2*p.y1_5 - 2*p.y1_4*M + p.y1_3*hair*p.M_2 - 2*p.y1*hair*p.M_4 + p.a_2*(p.y1_3 + hair*p.M_2*(p.y1 + M))*p.sin_y2_2)*p.sin_2y2))/(p.pow2_y1_2_p_a_2cos_y2_2*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr61(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (2*p.y1*p.pow3_y1_2_p_a_2cos_y2_2*(p.y1*(p.y1 - 2*M) + p.a_2*p.cos_y2_2) + p.a_2*M*(-2*p.y1_6 - p.y1_3*hair*p.M_2*(3*p.y1 + 2*M) - 2*p.y1*p.a_2*(p.y1_3 + hair*(p.y1 - M)*p.M_2)*p.cos_y2_2 + p.a_4*(2*p.y1_2 + hair*p.M_2)*p.cos_y2_4 + 2*p.a_6*p.cos_y2_6)*p.sin_y2_2)/(2.*p.pow2_y1_2_p_a_2cos_y2_2*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr62(double M, double a, double hair, const double *y, const ChristoffelParams &p)
{
    return (((p.y1_2 + p.a_2)*pow(2*p.y1_2 + p.a_2 + p.a_2*p.cos_2y2,3)*(p.a_2 + 2*p.y1*(p.y1 - 2*M) + p.a_2*p.cos_2y2)*cot(p.y2))/16. + p.y1*p.a_2*M*((p.a_2*p.cos_y2_3*(7*p.a_4 + 16*p.y1_2*(3*p.y1_2 - 2*p.y1*M + hair*p.M_2) + 8*p.a_2*(4*p.y1_2 - p.y1*M + hair*p.M_2) + 8*p.a_2*(p.a_2 + p.y1*(2*p.y1 - M))*p.cos_2y2 + p.a_4*p.cos_4y2)*p.sin_y2)/4. + p.y1*(2*p.y1_5 - 2*p.y1_4*M + p.y1_3*hair*p.M_2 - 2*p.y1*hair*p.M_4 + p.a_2*(p.y1_3 + hair*p.M_2*(p.y1 + M))*p.sin_y2_2)*p.sin_2y2))/(p.pow2_y1_2_p_a_2cos_y2_2*(-(p.y1_3*(p.y1_2 + p.a_2)*(-p.y1 + 2*M)) + p.a_2*(p.a_2*(p.y1_2 + p.a_2)*p.cos_y2_4 - p.y1*p.cos_y2_2*(p.a_2*(-2*p.y1 + M) + 2*p.y1_2*(-p.y1 + M) + p.a_2*M*p.cos_2y2) + p.y1*M*(2*p.y1_2 + hair*p.M_2)*p.sin_y2_2)));
}


double ChristoffelKerr(int rho, int mu, int nu, double M, double a, double hair,const double *y, const ChristoffelParams &p) /* double t,double r,double theta,double phi */
{
    static const std::vector< double (*)(double, double, double, const double*, const ChristoffelParams&) > chr = {
        NULL,
        &ChristoffelKerr1,
        &ChristoffelKerr2,
        NULL,
        &ChristoffelKerr4,
        NULL,
        NULL,
        &ChristoffelKerr7,
        &ChristoffelKerr8,
        NULL,
        NULL,
        &ChristoffelKerr11,
        NULL,
        &ChristoffelKerr13,
        &ChristoffelKerr14,
        NULL,
        &ChristoffelKerr16,
        NULL,
        NULL,
        &ChristoffelKerr19,
        NULL,
        &ChristoffelKerr21,
        &ChristoffelKerr22,
        NULL,
        NULL,
        &ChristoffelKerr25,
        &ChristoffelKerr26,
        NULL,
        &ChristoffelKerr28,
        NULL,
        NULL,
        &ChristoffelKerr31,
        &ChristoffelKerr32,
        NULL,
        NULL,
        &ChristoffelKerr35,
        NULL,
        &ChristoffelKerr37,
        &ChristoffelKerr38,
        NULL,
        NULL,
        &ChristoffelKerr41,
        &ChristoffelKerr42,
        NULL,
        &ChristoffelKerr44,
        NULL,
        NULL,
        &ChristoffelKerr47,
        NULL,
        &ChristoffelKerr49,
        &ChristoffelKerr50,
        NULL,
        &ChristoffelKerr52,
        NULL,
        NULL,
        &ChristoffelKerr55,
        &ChristoffelKerr56,
        NULL,
        NULL,
        &ChristoffelKerr59,
        NULL,
        &ChristoffelKerr61,
        &ChristoffelKerr62,
        NULL
    };
    
    
    int x = make_christoffel_index(rho, mu, nu);
    
    if (x >= chr.size()) {
        cout << "Undefined value in Christoffel:" << x << endl;
        return 0;
    }
    
    auto c = chr[x];
    if (c) return c(M, a, hair, y, p);
    
    return 0;
}

double metric(int mu, int nu, double M, double a, double epsilon3,double *y)
{
    
    double h=(pow(M,3)*y[1]*epsilon3)/pow(pow(y[1],2) + pow(a,2)*pow(cos(y[2]),2),2);
    
    int    x = nu+4.*mu;
    double res=0.;
    
    switch (x) {
        case 0: res = -1 - h + (2*(1 + h)*M*y[1])/(pow(y[1],2) + pow(a,2)*pow(cos(y[2]),2)); break;
            
        case 1: res = 0; break;
            
        case 2: res = 0; break;
            
        case 3: res = (-2*a*(1 + h)*M*y[1]*pow(sin(y[2]),2))/(pow(y[1],2) + pow(a,2)*pow(cos(y[2]),2)); break;
            
        case 4: res = 0; break;
            
        case 5: res = ((1 + h)*(pow(y[1],2) + pow(a,2)*pow(cos(y[2]),2)))/(pow(a,2) + y[1]*(-2*M + y[1]) + pow(a,2)*h*pow(sin(y[2]),2)); break;
            
        case 6: res = 0; break;
            
        case 7: res = 0; break;
            
        case 8: res = 0; break;
            
        case 9: res = 0; break;
            
        case 10: res = pow(y[1],2) + pow(a,2)*pow(cos(y[2]),2); break;
            
        case 11: res = 0; break;
            
        case 12: res = (-2*a*(1 + h)*M*y[1]*pow(sin(y[2]),2))/(pow(y[1],2) + pow(a,2)*pow(cos(y[2]),2)); break;
            
        case 13: res = 0; break;
            
        case 14: res = 0; break;
            
        case 15: res = pow(sin(y[2]),2)*(pow(a,2) + pow(y[1],2) + (pow(a,2)*(y[1]*(2*(1 + h)*M + h*y[1]) + pow(a,2)*h*pow(cos(y[2]),2))*pow(sin(y[2]),2))/(pow(y[1],2) + pow(a,2)*pow(cos(y[2]),2))); break;                        
        default:
            cout << "Undefined value in metric:"<<x;
    };
    return res;
}

void PF_eimu(double M, double a, double epsilon3, double *y ,double eimu[4][4], double emui[4][4]) // transforms BL to PF and vsv
{
    
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            eimu[i][j]=0.;
            emui[i][j]=0.;
        }
    }
    
    double r=y[1]*sin(y[2]);
    
    //define momentum coeficeints to find the energy and angular momentum
    
    double P4= pow((pow(a,2)* M *pow((pow(r,3) + pow(M,3)*epsilon3),6)* (4* pow(r,7) +16* pow(M,3)*pow( r,4)*epsilon3 - 6* pow(M,2)*pow(r,5)*epsilon3 +9 *pow(a,2)* pow(M,5)*pow(epsilon3,2))* pow((pow(r,4)* (-2* M + r) +pow(a,2) *(pow(r,3) +pow(M,3)*epsilon3)),2)),0.5);
    double P1=pow(a,2)* M*pow(r,4)* pow((pow(r,3) +pow(M,3)*epsilon3),2)* (12*pow(a,2)*pow(M,3)*epsilon3*pow((pow(r,3)+pow(M,3)*epsilon3),2)-pow(r,4)*(4*pow(r,6)* (-5 *M + 3 *r) +2* pow(M,2)* pow(r,3)* (-8* pow(M,2) + 3* pow(r,2))*epsilon3 + pow(M,5)* (40 *pow(M,2) - 48 *M* r + 15*pow(r,2))*pow(epsilon3,2)));
    
    double P2=2 *(2*pow(r,4)* (-M* P4 + pow(r,20)) +M *pow(r,12)* (2* pow(r,9)* (-12 *pow(M,2)+ 16* M* r - 7 *pow(r,2)) +pow(M,2)* pow((-2* M + r),2)*epsilon3* (-3* (8 *M - 3 *r)*pow(r,6) -6* pow(M,3)* (5* M - 2* r) *pow(r,3)*epsilon3+pow(M,6) *(-12* M + 5 *r)*pow(epsilon3,2))));
    
    double P3= -8* pow(a,2)* M* pow((pow(r,3) + pow(M,3)*epsilon3),2)* (2*pow(r,3) +5 *pow(M,3)*epsilon3) +pow(r,4)* pow((6* M *pow(r,3) - 2*pow(r,4) + 12* pow(M,4)*epsilon3 - 5* pow(M,3)* r *epsilon3),2);
    
    double P5= (epsilon3* pow(M,3) +pow(r,3))* (12*epsilon3*pow(a,6)*pow(M,3)* pow((epsilon3* pow(M,3) -2 *pow(r,3)),2)* pow((epsilon3 *pow(M,3) + pow(r,3)),4) +pow(a,4)*pow(r,4)*pow((epsilon3*pow(M,3) + pow(r,3)),2)* (-40*pow(epsilon3,4)* pow(M,13) +40*pow(epsilon3,4)*pow(M,12)* r - 15*pow(epsilon3,4)* pow(M,11)* pow(r,2) +128*pow(epsilon3,3)*pow(M,10)*pow(r,3) - 296*pow(epsilon3,3)*pow(M,9)*pow(r,4) +54.*pow(epsilon3,3)*pow(M,8)*pow(r,5) - 924*pow(epsilon3,2)*pow(M,7)*pow(r,6) +276*pow(epsilon3,2)*pow(M,6)*pow(r,7) - 36*pow(epsilon3,2)*pow(M,5)*pow(r,8) -880*epsilon3* pow(M,4)*pow(r,9) + 304*epsilon3*pow(M,3)*pow(r,10) -24*epsilon3*pow(M,2)*pow(r,11) - 112* M *pow(r,12) + 16* pow(r,13)) -2 *pow(a,2)*pow(r,8)* (48*pow(epsilon3,5)*pow(M,17) - 12*pow(epsilon3,5)*pow(M,16)* r -52*pow(epsilon3,5)*pow(M,15)*pow(r,2) +3*pow(epsilon3,4)*pow(M,14)*pow(r,3)*(5*epsilon3 + 88) -720*pow(epsilon3,4)*pow(M,13)*pow(r,4) + 298*pow(epsilon3,4)*pow(M,12)*pow(r,5) -3*pow(epsilon3,3)*pow(M,11)*pow(r,6)* (13*epsilon3 + 480) +516*pow(epsilon3,3)*pow(M,10)*pow(r,7) + 2 *pow(epsilon3,3)*pow(M,9)*pow(r,8) -6*pow(epsilon3,2)*pow(M,8)*pow(r,9)* (3*epsilon3 + 508) +2292*pow(epsilon3,2)*pow(M,7)* pow(r,10) - 628*pow(epsilon3,2)*pow(M,6)*pow(r,11) +12 *epsilon3*pow(M,5)*pow(r,12)* (5*epsilon3 - 134) +1188*epsilon3*pow(M,4)*pow(r,13) - 296*epsilon3*pow(M,3)*pow(r,14) + 24*pow(M,2)*pow(r,15)*(epsilon3 - 9) + 120* M*pow(r,16) - 16 *pow(r,17)) -pow(r,14)* pow((epsilon3*pow(M,3) + 6* M*pow(r,2) -2* pow(r,3)),2)* (96*pow(epsilon3,2)*pow(M,7) - 76*pow(epsilon3,2)*pow(M,6)*r +15*pow(epsilon3,2)*pow(M,5)*pow( r,2) + 72*epsilon3*pow(M,4)*pow(r,3) -44*epsilon3*pow(M,3)*pow(r,4) + 6*epsilon3*pow(M,2)*pow(r,5) + 12* M* pow(r,6) -4*pow(r,7))) -4*P4* (pow(a,2)* pow((epsilon3*pow(M,3) - 2* pow(r,3)),2)* (epsilon3* pow(M,3) + pow(r,3)) +6*epsilon3* pow(M,3)*pow(r,5)* (epsilon3* pow(M,3) + 6*M* pow(r,2) - 2*pow(r,3)));
    
    double P6= -6* M *pow(r,5) + 2*pow(r,6) -6*pow(M,4)*pow(r,2)*epsilon3 + pow(M,3)*pow(r,3)*epsilon3 - pow(M,6)*pow(epsilon3,2);
    
    double Energy=(1/pow(r,6))* pow(((P1 + P2)/P3),0.5);
    
    
    // in Kerr metric the root paprameter in the mathematica code is -1 for a<0.99 for a=0.99 root is 1.6126619934663224. To change it to the JP metric we need to include this root from mathematica code
    
    
    double root=-1;
    double sign=1;
    if (a==0.99) {
        root=1.6126619934663224;
    }
    if (a==0.9) {
        root=1.9340243143470757;
    }
    if (a==0.95) {
        root=1.781789048103268;
    }
    if (a==0.98) {
        root=1.6619345900612479;
    }
    if (a==0.998) {
        root=1.5673976204938782;
    }
    
    if (r<root) {
        sign=-1;
    }
    double Lz=(1/(pow(r,4)* P6))*(sign*pow(M*(pow(r,3) +epsilon3*pow(M,3))*(P5/P3),0.5) -6*a*M*(pow(r,3) +epsilon3*pow(M,3))*pow(((P1 +P2)/P3),0.5));
    
    // for equatorial plane theta=90
    
    double gg00=(1 - (2*M)/r)*(-1 - (pow(M,3)*epsilon3/pow(r,3)));
    double gg03= -((2*a*M*(pow(r,3)+pow(M,3)*epsilon3))/pow(r,4));
    double gg33=(pow(r,6)+sq(a)*(2*M+r)*(pow(r,3)+pow(M,3)*epsilon3))/pow(r,4);
    double gg11=(pow(r,5)+pow(M,3)*sq(r)*epsilon3)/(pow(r,4)*(-2*M+r)+sq(a)*(pow(r,3)+pow(M,3)*epsilon3));
    double gg22=sq(r);
    
    double g00i_equ=-1*((sq(r)*(pow(r,6)+sq(a)*(2*M+r)*(pow(r,3)+pow(M,3)*epsilon3)))/((pow(r,3)+pow(M,3)*epsilon3)*(pow(r,4)*(-2*M+r)+sq(a)*(pow(r,3)+pow(M,3)*epsilon3))));
    double g03i_equ=-1*((2* M* sq(r)* a)/(pow(r,4)*(-2*M+r)+ sq(a)*(pow(r,3)+pow(M,3)*epsilon3)));
    double g30i_equ=g03i_equ;
    double g33i_equ=(sq(r)*(-2*M+r))/(pow(r,4)*(-2*M+r)+pow(a,2)*(pow(r,3)+pow(M,3)*epsilon3));
    
    
    
    double uu0 = g00i_equ*(-1*Energy) + g03i_equ*Lz;
    double uu3 = g30i_equ*(-1*Energy) + g33i_equ*Lz;
        
    
    double AA1 = -(((gg03*uu0 + gg33*uu3)*sqrt(-(pow(gg00*uu0 + gg03*uu3,2)/((pow(gg03,2) - gg00*gg33)*(gg00*pow(uu0,2) + uu3*(2*gg03*uu0 + gg33*uu3))))))/(gg00*uu0 + gg03*uu3));
    
    double BB1 = sqrt(pow(gg00*uu0 + gg03*uu3,2)/((-pow(gg03,2) + gg00*gg33)*(gg00*pow(uu0,2) + uu3*(2*gg03*uu0 + gg33*uu3))));
    
    double BB=BB1;
    double AA=AA1;
    
    if (BB1<=0) {
        AA=-1*AA1;
        BB=-1*BB1;
    }
    
    
    
    emui[0][0]=uu0;
    emui[0][3]=AA;
    emui[1][1]=1/(pow(gg11,0.5));
    emui[2][2]=1/(pow(gg22,0.5));
    emui[3][0]=uu3;
    emui[3][3]=BB;
    
    
    eimu[0][0]=(emui[1][1]*emui[2][2]*emui[3][3])/(-emui[0][3]*emui[1][1]*emui[2][2]*emui[3][0] + emui[0][0]*emui[1][1]*emui[2][2]*emui[3][3]);
    eimu[0][3]=-1*((emui[0][3]*emui[1][1]*emui[2][2])/(-emui[0][3]*emui[1][1]*emui[2][2]*emui[3][0] + emui[0][0]*emui[1][1]*emui[2][2]*emui[3][3]));
    eimu[1][1]=(-emui[0][3]*emui[2][2]*emui[3][0] + emui[0][0]*emui[2][2]*emui[3][3])/(-emui[0][3]*emui[1][1]*emui[2][2]*emui[3][0] + emui[0][0]*emui[1][1]*emui[2][2]*emui[3][3]);
    eimu[2][2]=(-emui[0][3]*emui[1][1]*emui[3][0] + emui[0][0]*emui[1][1]*emui[3][3])/(-emui[0][3]*emui[1][1]*emui[2][2]*emui[3][0] + emui[0][0]*emui[1][1]*emui[2][2]*emui[3][3]);
    eimu[3][0]=-1*((emui[3][0]*emui[1][1]*emui[2][2])/(-emui[0][3]*emui[1][1]*emui[2][2]*emui[3][0] + emui[0][0]*emui[1][1]*emui[2][2]*emui[3][3]));
    eimu[3][3]=(emui[0][0]*emui[1][1]*emui[2][2])/(-emui[0][3]*emui[1][1]*emui[2][2]*emui[3][0] + emui[0][0]*emui[1][1]*emui[2][2]*emui[3][3]);
    
    
}



void CS_emui(double M, double a, double epsilon3,double *y,double emui[4][4])
{
    for(int i=0;i<4;i++)
        for(int j=0;j<4;j++)
            emui[i][j]=0.;
    
    emui[1][1]=1./sqrt(metric(1,1,M,a,epsilon3,y));
    emui[2][2]=1./sqrt(metric(2,2,M,a,epsilon3,y));
    
    double gg00=metric(0,0,M,a,epsilon3,y);
    double gg03=metric(0,3,M,a,epsilon3,y);
    double gg33=metric(3,3,M,a,epsilon3,y);
    
    double AA=-(gg03/sqrt(gg00*(-sq(gg03) + gg00*gg33)));
    double BB=sqrt(gg00/(-sq(gg03) + gg00*gg33));
    emui[0][0]=1./sqrt(fabs(gg00));
    emui[0][3]=AA;
    emui[3][3]=BB;
}

void Scatter(Photon *gamma_p, double u[4], double f_corona[4], double Tcorona, double emunu[4][4], double enumu[4][4], double fac){
    
    
    //transform u ZAMO
    double uZAMO[4];
    double fZAMO[4];
    
    for (int i=0;i<4;i++)
    {
        uZAMO[i]=0.;
        fZAMO[i]=0;
        for (int mu=0; mu<4; mu++)
        {
            uZAMO[i] += enumu[i][mu]*u[mu];
            fZAMO[i] += enumu[i][mu]*f_corona[mu];
        }
    }
    
    
    if (uZAMO[0]<0) {
        if(gamma_p->E>0){
            cout<<"ZAMO photons number "<< gamma_p->totalN<<endl;
        }
        else{
            uZAMO[0]*=-1;
            uZAMO[1]*=-1;
            uZAMO[2]*=-1;
            uZAMO[3]*=-1;
        }
        
    }
    
    double mag1= -1*uZAMO[0]*fZAMO[0]+uZAMO[1]*fZAMO[1]+uZAMO[2]*fZAMO[2]+uZAMO[3]*fZAMO[3];
    double mag2= -1*uZAMO[0]*uZAMO[0]+uZAMO[1]*uZAMO[1]+uZAMO[2]*uZAMO[2]+uZAMO[3]*uZAMO[3];
    if(abs(mag1)>0.1 || abs(mag2)>0.1)
        cout<<"mag1 in corona "<<mag1<<" mag2 "<<mag2<<" number "<<gamma_p->totalN<<endl;
    

    
    //give electron random direction
    double cos_theta_e=1-(2*ran1(&seed));
   // double theta_e = PI*ran1(&seed);
    double phi_e = 2.*PI*ran1(&seed);
    
    double n[4]={0.,0.,0.,0.};
    n[1] = sqrt(1-sq(cos_theta_e))*cos(phi_e);
    n[2] = -cos_theta_e;
    n[3] = sqrt(1-sq(cos_theta_e))*sin(phi_e);

    
    
    //Find Gamma (Lorentz Factor of Electron) Using distribution (Maxwell-Boltzman) (using Schnittman&Krolik 2013, Appendix B)
    
    double gamma_e= Gamma_MB(Tcorona);
    double beta_e = sqrt(1 - 1/pow(gamma_e,2));
    
    
    // find angle between elecron and photon before scattering in ZAMO frame
    
    double Cos_angle_oframe = 0;
    double ph_mag_oframe = sqrt(pow(uZAMO[1],2)+pow(uZAMO[2],2)+pow(uZAMO[3],2));
    double n_mag_oframe = sqrt(pow(n[1],2)+pow(n[2],2)+pow(n[3],2));
    
    for (int i=1; i<4; i++)
    Cos_angle_oframe += uZAMO[i]*n[i];
    
    Cos_angle_oframe/=(ph_mag_oframe*n_mag_oframe);
    
    
    double  weight_e_p=1-beta_e*Cos_angle_oframe;  
    
    //construct Lorentz transformation and inverse (Misner,Thorne & Wheeler 1973, pg.69)
    double LAMBDA[4][4];
    double LAMBDAi[4][4];  //inverse of LAMBDA
    
    Lorentz(gamma_e, beta_e, n, LAMBDA, LAMBDAi); // Lorentz Transformation Metrix
    
    
    //transform u and f into electron Rest Frame
    double uERF[4];
    double fERF[4];
    
    for (int i=0; i<4; i++)
    {
        uERF[i]=0.;
        fERF[i]=0.;

        for (int j=0; j<4; j++)
        {
            uERF[i] += LAMBDA[i][j]*uZAMO[j];
            fERF[i] += LAMBDA[i][j]*fZAMO[j];

        }
    }
    
    for (int i=1; i<4; i++) {
        fERF[i]=fERF[i]-(fERF[0]/uERF[0])*uERF[i];
    }
    fERF[0]=0;
    
    //give photon a random scattering direction in ERF frame
    double cos_theta_ERF=1-(2*ran1(&seed));
    //double theta_ERF = PI*ran1(&seed);
    double phi_ERF = 2.*PI*ran1(&seed);
    
    double m[4]={1.,0.,0.,0.};
    m[1] = sqrt(1-sq(cos_theta_ERF))*cos(phi_ERF);
    m[2] = -cos_theta_ERF;
    m[3] = sqrt(1-sq(cos_theta_ERF))*sin(phi_ERF);
    
    
    //FIND ANGLE BETWEEN uERF AND M[]
    double um_Cos_angle = 0;
    double uERF_mag = sqrt(pow(uERF[1],2)+pow(uERF[2],2)+pow(uERF[3],2));
    double m_mag = 1;
    
    for (int i=1; i<4; i++)
        um_Cos_angle += (uERF[i]*m[i]);
    
    um_Cos_angle /= (uERF_mag*m_mag);
    
    double theta_um=acos(um_Cos_angle);
    
    // define basis and finding the f vecotr in these basis
    
    double k1n[3];
    normalize(&uERF[1],k1n);
    
    double k2n[3];
    normalize(&m[1],k2n);
    
    double e_perp[3];
    cross(k2n,k1n,e_perp);
    normalize(e_perp,e_perp);
    
    
    double e_para[3];
    cross(k1n,e_perp,e_para);
    normalize(e_para,e_para);
    
    
    double f1n[3];
    normalize(&fERF[1],f1n);
    
    double Cos_polAngl=dot(f1n,e_para);
    double Sin_polAngl=dot(f1n,e_perp);
    
    double jam=sq(Cos_polAngl)+sq(Sin_polAngl);
    if((jam<0.95)||(jam>1.05)){
        cout<<" jam "<<sq(Cos_polAngl)+sq(Sin_polAngl)<<" num "<<gamma_p->totalN<<endl;
        gamma_p->timeout=1000001;
    }
    
    
    double fERF_test[3];
    for (int i=0; i<3; i++) {
        fERF_test[i]=Cos_polAngl*e_para[i]+Sin_polAngl*e_perp[i];
    }
    if (fabs(dot(fERF_test,k1n))>0.01) {
        cout<<" There is a problem in finding the normal vector for fERF in first step "<<endl;
        exit(1);
    }
    
    double PolAngl=atan2(Sin_polAngl,Cos_polAngl);
    if(PolAngl<0)
        PolAngl+=PI;
    
    double stokes_noscat[3]={1, (gamma_p->pol *cos(2*PolAngl)) , (gamma_p->pol *sin(2*PolAngl))};
    
    double testchi=getChi(stokes_noscat[1],stokes_noscat[2]);
    
    if (abs(testchi-PolAngl)>0.01) {
        cout<<"In corona chi me "<<PolAngl<<" chitest "<<testchi<<endl;
        
    }
    
    // finding new energy
    double Escale=(uERF[0]*abs(gamma_p->E))/511.;
    
    double f=1/(1.+Escale*(1-um_Cos_angle));
    double uERF0save=uERF[0]*abs(gamma_p->E);
    
    for (int i=0; i<4; i++)
        uERF[i]=f*uERF0save*m[i];
    
    
    //defin FANO matrix
    
    double FANO[4][4];
    
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            FANO[i][j]=0;
        }
    }
    FANO[0][0]=pow(um_Cos_angle,2);
    FANO[1][1]=1;
    FANO[2][2]=um_Cos_angle;
    FANO[3][3]=um_Cos_angle;
    
    
    // New stokes after scattering
    double I_noscat[4]={0.5*(stokes_noscat[0]+stokes_noscat[1]),0.5*(stokes_noscat[0]-stokes_noscat[1]), stokes_noscat[2],0};
    double I_scat[4];
    
    
    
    for (int i=0; i<4; i++) {
        I_scat[i]=0;
        for (int j=0; j<4; j++) {
            I_scat[i] += FANO[i][j]*I_noscat[j];
        }
    }
    
    double stokes_scat[3]={I_scat[0]+I_scat[1],I_scat[0]-I_scat[1],I_scat[2]};
    // Use rejection method to deal with Klein Nishina regime:
    double normFac=stokes_scat[0];
    if (fac<normFac) {
        // In the Thomson regime, this will reject 1/3 of all events.
        // This is good because it gives the right cos(theta) distribution,
        // however, we have to make up for the lost events by
        // weighing the events that make it with the weight 1.5
        
        
        
        // define the new basis with scatt
        double e_perp_prime[3];
        for (int i=0; i<3; i++) {
            e_perp_prime[i]=e_perp[i];
        }
        
        double e_para_prime[3];
        cross(k2n,e_perp_prime,e_para_prime);
        normalize(e_para_prime,e_para_prime);
        
        PolAngl= getChi(stokes_scat[1],stokes_scat[2]);
        
        double fERF_scat[4];
        fERF_scat[0]=0;
        for (int i=1; i<4; i++) {
            fERF_test[i-1]=cos(PolAngl)*e_para_prime[i-1]+sin(PolAngl)*e_perp_prime[i-1];
            fERF_scat[i]=cos(PolAngl)*e_para_prime[i-1]+sin(PolAngl)*e_perp_prime[i-1];
        }
        if (fabs(dot(fERF_test,k2n))>0.01) {
            cout<<" There is a problem in finding the normal vector for fERF_scat in second step "<<endl;
            exit(1);
        }
        
        gamma_p->pol = sqrt(pow((stokes_scat[1]),2)+pow(stokes_scat[2],2))/stokes_scat[0];
        
        //normalization weight is:
        double normFac=(((2*uERF[0]*(2+uERF[0]*(1+uERF[0])*(8+uERF[0])))/(pow(1+2*uERF[0],2)))+(-2+(-2+uERF[0])*uERF[0])*log(1+2*uERF[0]))/(4*pow(uERF[0],3));
        
        
        
        
        if (stokes_scat[0]<0) {
            
            cout<<" stokes_scat[0] <0 "<<stokes_scat[0]<<" ph number "<<gamma_p->totalN<<endl;
            exit(1);
            
        }
        
        //transform u and f back into ZAMO
        for (int i=0; i<4; i++)
        {
            uZAMO[i]=0.;
            fZAMO[i]=0;
            for (int j=0; j<4; j++)
            {
                uZAMO[i] += LAMBDAi[i][j]*uERF[j];
                fZAMO[i] += LAMBDAi[i][j]*fERF_scat[j];
                
            }
            
        }
        
        
        
        //transform u and f back into BL coordinates
        for (int mu=0;mu<4;mu++)
        {
            u[mu]=0.;
            f_corona[mu]=0;
            for (int i=0; i<4; i++)
            {
                u[mu] += emunu[mu][i] *uZAMO[i];
                f_corona[mu] += emunu[mu][i] *fZAMO[i];
            }
        }
        
        
        gamma_p->weight *= weight_e_p;
        
        gamma_p->nScatter_Corona=1+gamma_p->nScatter_Corona;
        
    }
}


void BZ_emunu(double M, double a, double epsilon3, double *y ,double emunu[4][4], double enumu[4][4]) // emunu transforms BL to ZAMO, enumu, inverse metrix, transfroms ZAMO to BL
{
    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            emunu[i][j]=0;
            enumu[i][j]=0;
        }
    }
    
    double Sigma_metric= sq(y[1])+sq(a*cos(y[2]));
    double Delta_metric= sq(y[1])-2*M*y[1]+sq(a);
    double A_metric=sq(sq(y[1])+sq(a))-sq(a)*Delta_metric*sq(sin(y[2]));
    
    double rho_sq=sq(y[1])+sq(a*cos(y[2]));
    double delta=sq(y[1])-2*M*y[1]+sq(a);
    double alpha=sqrt(rho_sq*delta/(rho_sq*delta+2*M*y[1]*(sq(a)+sq(y[1]))));
    double omega=2*M*y[1]*a/(rho_sq*delta+2*M*y[1]*(sq(a)+sq(y[1])));
    double omegabar_sq=(rho_sq*delta+2*M*y[1]*(sq(a)+sq(y[1])))*sq(sin(y[2]))/rho_sq;
    
    emunu[0][0]=1/alpha;
    emunu[1][1]=sqrt(delta/rho_sq);
    emunu[2][2]=sqrt(1/rho_sq);
    emunu[3][3]=sqrt(1/omegabar_sq);
    emunu[3][0]=omega/alpha;
    
    enumu[0][0]=sqrt((Sigma_metric*Delta_metric)/A_metric);
    enumu[1][1]=sqrt(Sigma_metric/Delta_metric);
    enumu[2][2]=sqrt(Sigma_metric);
    enumu[3][3]=sqrt(A_metric/Sigma_metric)*sin(y[2]);
    enumu[3][0]=-2*M*a*y[1]*sin(y[2])/sqrt(Sigma_metric*A_metric);
    
}

void Lorentz(double gamma_e, double beta_e, double *n ,double LAMBDA[4][4], double LAMBDAi[4][4]) // Lorentz Transformation Metrix
{
    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            LAMBDAi[i][j]=0;
            LAMBDA[i][j]=0;
        }
    }
    
    for (int j=1; j<4; j++)
    {
        LAMBDA[0][j] = -beta_e*gamma_e*n[j];
        LAMBDA[j][0] = LAMBDA[0][j];
        LAMBDA[j][j] = (gamma_e - 1)*n[j]*n[j] + 1;
        
        LAMBDAi[0][j] = beta_e*gamma_e*n[j]; //beta -> -beta
        LAMBDAi[j][0] = LAMBDAi[0][j];
        LAMBDAi[j][j] = LAMBDA[j][j];
    }
    LAMBDA[0][0] = gamma_e;
    LAMBDA[1][2] = (gamma_e - 1)*n[1]*n[2];
    LAMBDA[1][3] = (gamma_e - 1)*n[1]*n[3];
    LAMBDA[2][1] = LAMBDA[1][2];
    LAMBDA[2][3] = (gamma_e - 1)*n[2]*n[3];
    LAMBDA[3][1] = LAMBDA[1][3];
    LAMBDA[3][2] = LAMBDA[2][3];
    
    LAMBDAi[0][0] = LAMBDA[0][0];
    LAMBDAi[1][2] = LAMBDA[1][2];
    LAMBDAi[1][3] = LAMBDA[1][3];
    LAMBDAi[2][1] = LAMBDAi[1][2];
    LAMBDAi[2][3] = LAMBDA[2][3];
    LAMBDAi[3][1] = LAMBDAi[1][3];
    LAMBDAi[3][2] = LAMBDAi[2][3];
    
}

double Gamma_MB(double Tcorona)  // find gamma using Mawell- Boltzman distribution
{
    int N=1000;
    double THETA = Tcorona/511; //511 comes from converting mec^2 in jouls to kev because THETA is kT_e/m_ec^2
    double gamma_e=0;
    for (int i=0; i<N; i++) {
        
        double lambda0 = ran1(&seed);
        double x0= RootFinder(THETA,lambda0);
        double beta_e = sqrt(1 - 1/pow(x0,2));
        
        double lambda1 =ran1(&seed);
        
        if (beta_e>lambda1) {
            gamma_e=x0;
            i=N;
        }
        
    }
    
    if (gamma_e==0) {
        cout<<" Cannot find any Gamma for N="<<N<<endl;
    }
    else{
        return gamma_e;
    }
    
}


// === general purpose functions ===

void load_table(double *t_x,int &numbins, char *fname1,int checknum)
/* load table from file */
{        
    cout <<"Opening "<<fname1<<"... ";
    
    ifstream fp_in;
    
    fp_in.open(fname1,ios::in); 
    
    numbins=0;
    double val;
    
    while (fp_in >>val)
		t_x[numbins++]=val;
    
    cout <<"Read "<<numbins<<" values: "<<"Val[0]= "<<t_x[0]<<" ... Val["<<numbins-1<<"]="<<t_x[numbins-1]<<" . "<<endl;
    
	fp_in.close();
    
    if (numbins!=checknum) {cout << "The table does not fit:"<<numbins<<" read but "<<checknum<<" expected "<<endl; exit(1);}    
}

double get_value_log(double *t_x,double *t_v,int numbins,double x)
{
	int i;
	double rest,res;
	double delta;
    
	if ((x<t_x[0])||(x>t_x[numbins-1]))
    {
        printf("Error in get_table : %f %f %f\n",x,t_x[0],t_x[numbins-1]);
        exit (-1);
    }
	
	for (i=0;i<numbins-1;i++)
		if (x<t_x[i]) break;
    
	delta =  log(t_x[i])-log(t_x[i-1]);
	rest  = (log(x     )-log(t_x[i-1]))/delta;
	res   = (1.-rest)*t_v[i-1] + rest*t_v[i];
	return res;
}

double get_value_lin(double *t_x,double *t_v,int numbins,double x)
{
	int i;
	double rest,res;
	double delta;
    
	if ((x<t_x[0])||(x>t_x[numbins-1]))
    {
        printf("Error in get_table : %f %f %f\n",x,t_x[0],t_x[numbins-1]);
        exit (-1);
    }
	
	for (i=0;i<numbins-1;i++)
		if (x<t_x[i]) break;
    
	delta =  t_x[i]-t_x[i-1];
	rest  = (x-t_x[i-1])/delta;
	res   = (1.-rest)*t_v[i-1] + rest*t_v[i];
	return res;
}

int linBin(double x1,double x2,int numbins,double x)
{
    double delta=(x2-x1)/(double)numbins;
    
 	if ((x<x1)||(x>x2))
    {
        printf("Error in linBin : %f %f %f\n",x1,x2,x);
        exit (-1);
    }
    
    return (int)((x-x1)/delta);
}

int logBin(double x1,double x2,int numbins,double x)
{
    double delta=(log(x2)-log(x1))/(double)numbins;
    
 	if ((x<x1)||(x>x2))
    {
        printf("Error in logBin : %f %f %f\n",x1,x2,x);
        exit (-1);
    }
    
    return (int)((log(x)-log(x1))/delta);
}

double binLog(double x1,double x2,int numbins,int nbin)
{
    double delta=(log(x2)-log(x1))/(double)numbins;
    double lval =exp(log(x1)+delta*((double)nbin+0.5));
    
    return lval;
}


double RootFinder(double THETA,double lambda0) //newton method
{
    double x0= 1.+THETA;     //giving the first estimate for root //I checked before the root is around this
    double xi=x0;
    double xi_1;
    double gamma;
    double E=0.001;   //error for the root, So this means that the answer is 0 for this root by the error of 0.001
    double root=0.;
    int N=2000;  //The number of loop for newton method
    
    for (int i=0; i<N; i++) {
        
        gamma=xi;
        
        double f_xi=(1-(exp(-gamma/THETA)*(2*sq(THETA)+2*THETA*gamma+sq(gamma))/(exp(-1/THETA)*(2*sq(THETA)+2*THETA+1))))-lambda0;
        
        double df_xi=(1/(exp(-1/THETA)*(2*sq(THETA)+2*THETA+1)))*exp(-gamma/THETA)*(sq(gamma)/THETA);
        
        xi_1= xi-(f_xi/df_xi);
        
        gamma=xi_1;
        
        double f_xi_1=(1-(exp(-gamma/THETA)*(2*sq(THETA)+2*THETA*gamma+sq(gamma))/(exp(-1/THETA)*(2*sq(THETA)+2*THETA+1))))-lambda0;
        
        if (abs(f_xi_1)<E) {
            root=gamma;
            i=N;
        }
        
        else{
            xi=xi_1;
        }
        
    }
    
    if (root==0) {
        cout<<" root not found "<<"lambda0 "<<lambda0<<" Theta "<<THETA<<endl;
    }
    else{
        return root;
    }
}

double cot(double x) 
{ return (1./tan(x)); }

void Matrix_SMult(double c,double a1[][3], double a2[][3])
{
    for(int i = 0; i < 3; i++) 
        for(int j = 0; j < 3; j++) 
            a2[i][j] =  c*a1[i][j];
}

void Matrix_Add(double a1[][3], double a2[][3], double a3[][3])
{
    for(int i = 0; i < 3; i++) 
        for(int j = 0; j < 3; j++) 
            a3[i][j] =  a1[i][j]+a2[i][j];
}

void Matrix_Mult(double a1[][3], double a2[][3], double a3[][3])
{
    double d[3][3];
    
    for(int i = 0; i < 3; i++) 
        for(int j = 0; j < 3; j++) 
            d[i][j] =  dot3(a1, a2, i, j);
    
    for(int i = 0; i < 3; i++) 
        for(int j = 0; j < 3; j++) 
            a3[i][j]=d[i][j];
}

double dot3(double a1[][3], double a2[][3], int aRow, int bCol)
{
    double sum = 0;
	for(int k = 0; k < 3; k++)
	    sum += a1[aRow][k] * a2[k][bCol];
    return sum;
}

void MatrixVector_Mult(double a1[][3], double v[3], double r[3])
{
    double d[3];
    for(int i = 0; i < 3; i++) {
        d[i] = 0;
        for(int k = 0; k < 3; k++)
            d[i] += a1[i][k] * v[k];
    }
    for(int i = 0; i < 3; i++) r[i]=d[i];
}

void Matrix_Print(double *F,int n, int l)
{
    for (int i=0;i<n;i++,cout<<endl)
        for (int j=0;j<l;j++)
            cout<<F[i*l+j]<<" ";
}

double dot(double *a, double *b)
{
    double sum=0.;
	for(int i = 0; i < 3; i++)
	    sum+=a[i]*b[i];
    return sum;
}

void normalize(double *v,double *n)
{
    double sum=0.;
    
    for (int i=0;i<3;i++) sum+=sq(v[i]);
    sum=sqrt(sum);
    
    for (int i=0;i<3;i++) n[i] = v[i]/sum;
}

void perp(double *v,double *k,double *p)
{
    double knorm[3];
    normalize(k,knorm);
    
    double s=dot(v,knorm);
    
    for (int i=0;i<3;i++) p[i]=v[i]-s*knorm[i];
    normalize(p,p);
}

void cross(double *a, double *b, double *c)
{
    double d[3];
    d[0]=a[1]*b[2]-a[2]*b[1];
    d[1]=a[2]*b[0]-a[0]*b[2];
    d[2]=a[0]*b[1]-a[1]*b[0];
    for (int i=0;i<3;i++) c[i]=d[i];
}

double angle(double *aNorm,double *bNorm)
{
    double s=dot(aNorm,bNorm);
    return acos(s);
}

double getChi(double q, double u)
{
    double res = atan(u/q)/2.;
    if(q < 0.) res += PI/2.;
    else if (u < 0.) res += PI;
    return res;
}


// to simulate planck distribution I ussed the following functions
// These functions are from the website: https://people.sc.fsu.edu/~jburkardt/cpp_src/prob/prob.cpp
//

double planck_sample ( double a, double b, int &seed )
{
    double a2;
    double b2;
    double c2;
    double g;
    double x;
    int z;
    //
    a2 = 0.0;
    b2 = 1.0;
    c2 = b + 1.0;
    
    g = gamma_sample ( a2, b2, c2, seed );
    
    z = zipf_sample ( c2, seed );
    
    x = g / ( a * ( double ) ( z ) );
    
    return x;
}

int zipf_sample ( double a, int &seed )
{
    double b;
    double t;
    double u;
    double v;
    double w;
    int x;
    
    b = pow ( 2.0, ( a - 1.0 ) );
    
    for ( ; ; )
    {
        u = r8_uniform_01 ( seed );
        v = r8_uniform_01 ( seed );
        w = ( int ) ( 1.0 / pow ( u, 1.0 / ( a - 1.0 ) ) );
        
        t = pow ( ( w + 1.0 ) / w, a - 1.0 );
        
        if ( v * w * ( t - 1.0 ) * b <= t * ( b - 1.0 ) )
        {
            break;
        }
        
    }
    
    x = ( int ) w;
    
    return x;
}

double gamma_sample ( double a, double b, double c, int &seed )
{
    double a1 =   0.3333333;
    double a2 = - 0.2500030;
    double a3 =   0.2000062;
    double a4 = - 0.1662921;
    double a5 =   0.1423657;
    double a6 = - 0.1367177;
    double a7 =   0.1233795;
    double bcoef;
    double co;
    double d;
    double e;
    double e1 = 1.0;
    double e2 = 0.4999897;
    double e3 = 0.1668290;
    double e4 = 0.0407753;
    double e5 = 0.0102930;
    double euler = 2.71828182845904;
    double p;
    double q;
    double q0;
    double q1 =  0.04166669;
    double q2 =  0.02083148;
    double q3 =  0.00801191;
    double q4 =  0.00144121;
    double q5 = -0.00007388;
    double q6 =  0.00024511;
    double q7 =  0.00024240;
    double r;
    double s;
    double si;
    double s2;
    double t;
    double u;
    double v;
    double w;
    double x;
    //
    //  Allow C = 0.
    //
    if ( c == 0.0 )
    {
        x = a;
        return x;
    }
    //
    //  C < 1.
    //
    if ( c < 1.0 )
    {
        for ( ; ; )
        {
            u = r8_uniform_01 ( seed );
            t = 1.0 + c / euler;
            p = u * t;
            
            s = exponential_01_sample ( seed );
            
            if ( p < 1.0 )
            {
                x = exp ( log ( p ) / c );
                if ( x <= s )
                {
                    break;
                }
            }
            else
            {
                x = - log ( ( t - p ) / c );
                if ( ( 1.0 - c ) * log ( x ) <= s )
                {
                    break;
                }
            }
        }
        
        x = a + b * x;
        return x;
    }
    //
    //  1 <= C.
    //
    else
    {
        s2 = c - 0.5;
        s = sqrt ( c - 0.5 );
        d = sqrt ( 32.0 ) - 12.0 * sqrt ( c - 0.5 );
        
        t = normal_01_sample ( seed );
        x = pow ( ( sqrt ( c - 0.5 ) + 0.5 * t ), 2 );
        
        if ( 0.0 <= t )
        {
            x = a + b * x;
            return x;
        }
        
        u = r8_uniform_01 ( seed );
        
        if ( d * u <= t * t * t )
        {
            x = a + b * x;
            return x;
        }
        
        r = 1.0 / c;
        
        q0 = ( ( ( ( ( (q7*r+ q6 )*r+ q5 )*r+ q4 )*r+ q3 )*r+ q2 )*r+ q1) * r;
        
        if ( c <= 3.686 )
        {
            bcoef = 0.463 + s - 0.178 * s2;
            si = 1.235;
            co = 0.195 / s - 0.079 + 0.016 * s;
        }
        else if ( c <= 13.022 )
        {
            bcoef = 1.654 + 0.0076 * s2;
            si = 1.68 / s + 0.275;
            co = 0.062 / s + 0.024;
        }
        else
        {
            bcoef = 1.77;
            si = 0.75;
            co = 0.1515 / s;
        }
        
        if ( 0.0 < sqrt ( c - 0.5 ) + 0.5 * t )
        {
            v = 0.5 * t / s;
            
            if ( 0.25 < r8_abs ( v ) )
            {
                q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v );
            }
            else
            {
                q = q0 + 0.5 * t * t * ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
            }
            
            if ( log ( 1.0 - u ) <= q )
            {
                x = a + b * x;
                return x;
            }
        }
        
        for ( ; ; )
        {
            e = exponential_01_sample ( seed );
            
            u = r8_uniform_01 ( seed );
            
            u = 2.0 * u - 1.0;
            t = bcoef + r8_abs ( si * e ) * r8_sign ( u );
            
            if ( -0.7187449 <= t )
            {
                v = 0.5 * t / s;
                
                if ( 0.25 < r8_abs ( v ) )
                {
                    q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v );
                }
                else
                {
                    q = q0 + 0.5 * t * t * ( ( ( ( ( (a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1) * v;
                }
                
                if ( 0.0 < q )
                {
                    if ( 0.5 < q )
                    {
                        w = exp ( q ) - 1.0;
                    }
                    else
                    {
                        w = ( ( ( (e5   * q+ e4 ) * q+ e3 ) * q+ e2 ) * q+ e1 ) * q;
                    }
                    
                    if ( co * r8_abs ( u ) <= w * exp ( e - 0.5 * t * t ) )
                    {
                        x = a + b * pow ( s + 0.5 * t, 2 );
                        return x;
                    }
                }
            }
        }
    }
}

double r8_uniform_01 ( int &seed )
{
    int k;
    double r;
    
    k = seed / 127773;
    
    seed = 16807 * ( seed - k * 127773 ) - k * 2836;
    
    if ( seed < 0 )
    {
        seed = seed + 2147483647;
    }
    
    r = ( double ) ( seed ) * 4.656612875E-10;
    
    return r;
}

double exponential_01_sample ( int &seed )
{
    double cdf;
    double x;
    
    cdf = r8_uniform_01 ( seed );
    
    x = - log ( 1.0 - cdf );
    
    return x;
}

double normal_01_sample ( int &seed )
{
    
    const double pi = 3.14159265358979323;
    double r1;
    double r2;
    double x;
    
    r1 = r8_uniform_01 ( seed );
    r2 = r8_uniform_01 ( seed );
    x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );
    
    return x;
}

double r8_abs ( double x )
{
    double value;
    
    if ( 0.0 <= x )
    {
        value = x;
    }
    else
    {
        value = -x;
    }
    return value;
}

double r8_sign ( double x )

{
    if ( x < 0.0 )
    {
        return ( -1.0 );
    }
    else
    {
        return ( 1.0 );
    }
}


/* Random number generator ran1 from Computers in Physics */
/* Volume 6 No. 5, 1992, 522-524, Press and Teukolsky */
/* To generate real random numbers 0.0-1.0 */
/* Should be seeded with a negative integer */
#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836
#define NTAB 32
#define EPS (1.2E-07)
#define MAX(a,b) (a>b)?a:b
#define MIN(a,b) (a<b)?a:b

double ran1(long *idum)
{
    int j,k;
	static int iv[NTAB],iy=0;
	void nrerror();
    static double NDIV = 1.0/(1.0+(IM-1.0)/NTAB);
    static double RNMX = (1.0-EPS);
    static double AM = (1.0/IM);
    
	if ((*idum <= 0) || (iy == 0)) {
		*idum = MAX(-*idum,*idum);
        for(j=NTAB+7;j>=0;j--) {
			k = *idum/IQ;
			*idum = IA*(*idum-k*IQ)-IR*k;
			if(*idum < 0) *idum += IM;
            if(j < NTAB) iv[j] = *idum;
        }
		iy = iv[0];
	}
	k = *idum/IQ;
	*idum = IA*(*idum-k*IQ)-IR*k;
	if(*idum<0) *idum += IM;
    j = iy*NDIV;
    iy = iv[j];
    iv[j] = *idum;
    return MIN(AM*iy,RNMX);
}
#undef IA 
#undef IM 
#undef IQ
#undef IR
#undef NTAB
#undef EPS 
#undef MAX
#undef MIN

/************************************************** **********************/
/* Please Note: */
/* */
/* (1) This computer program is written by Tao Pang in conjunction with */
/* his book, "An Introduction to Computational Physics," published */
/* by Cambridge University Press in 1997. */
/* */
/* (2) No warranties, express or implied, are made for this program. */
/* */
/************************************************** **********************/

void migs (double a[][NMAX],int n,double x[][NMAX])
{
    int i,j,k;
    double b[NMAX][NMAX];
    int indx[NMAX];
    
    if (n > NMAX)
    {
        printf("The matrix dimension is too large.\n");
        exit(1);
    }
    
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            b[i][j] = 0;
        }
    }
    for (i = 0; i < n; ++i)
    {
        b[i][i] = 1;
    }
    
    elgs (a,n,indx);
    
    for (i = 0; i < n-1; ++i)
    {
        for (j = i+1; j < n; ++j)
        {
            for (k = 0; k < n; ++k)
            {
                b[indx[j]][k] = b[indx[j]][k]-a[indx[j]][i]*b[indx[i]][k];
            }
        }
    }
    
    for (i = 0; i < n; ++i)
    {
        x[n-1][i] = b[indx[n-1]][i]/a[indx[n-1]][n-1];
        for (j = n-2; j >= 0; j = j-1)
        {
            x[j][i] = b[indx[j]][i];
            for (k = j+1; k < n; ++k)
            {
                x[j][i] = x[j][i]-a[indx[j]][k]*x[k][i];
            }
            x[j][i] = x[j][i]/a[indx[j]][j];
        }
    }
}

void elgs (double a[][NMAX],int n,int indx[NMAX])

/* Function to perform the partial-pivoting Gaussian elimination.
 a[][] is the original matrix in the input and transformed
 matrix plus the pivoting element ratios below the diagonal
 in the output. indx[] records the pivoting order.
 Copyright (c) Tao Pang 2001. */
{
    int i, j, k, itmp;
    double c1, pi, pi1, pj;
    double c[NMAX];
    
    if (n > NMAX)
    {
        printf("The matrix dimension is too large.\n");
        exit(1);
    }
    
    /* Initialize the index */
    
    for (i = 0; i < n; ++i)
    {
        indx[i] = i;
    }
    
    /* Find the rescaling factors, one from each row */
    
    for (i = 0; i < n; ++i)
    {
        c1 = 0;
        for (j = 0; j < n; ++j)
        {
            if (fabs(a[i][j]) > c1) c1 = fabs(a[i][j]);
        }
        c[i] = c1;
    }
    
    /* Search the pivoting (largest) element from each column */ 
    
    for (j = 0; j < n-1; ++j)
    {
        pi1 = 0;
        for (i = j; i < n; ++i)
        {
            pi = fabs(a[indx[i]][j])/c[indx[i]];
            if (pi > pi1)
            {
                pi1 = pi;
                k = i;
            }
        }
        
        /* Interchange the rows via indx[] to record pivoting order */
        
        itmp = indx[j];
        indx[j] = indx[k];
        indx[k] = itmp;
        for (i = j+1; i < n; ++i)
        {
            pj = a[indx[i]][j]/a[indx[j]][j];
            
            /* Record pivoting ratios below the diagonal */
            
            a[indx[i]][j] = pj;
            
            /* Modify other elements accordingly */
            
            for (k = j+1; k < n; ++k)
            {
                a[indx[i]][k] = a[indx[i]][k]-pj*a[indx[j]][k];
            }
        }
    }
}
