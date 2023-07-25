// for data files with  Emis out put
//Finding polarization with mean of 5 subsets

// C++
#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <iomanip>   // format manipulation
#include <string>
#include <cmath>

// C
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>    

// Root
#include <TROOT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2.h>
#include <TLine.h>
#include <TFile.h>
#include <TChain.h>
#include "TStyle.h"
#include "TROOT.h"
#include "TColor.h"
#include "TFile.h"


using namespace std;

#include "jpsch.h"
#include "macro.h"

static long seed = -99;
double fitf(double *v, double *par);


// remove the "NO_" to fill a tree with the information of individual tracks
// you can look at the information with tree->Draw() ...
#define No_TreeFill

void wait()
{
	char in[50];
	cout << "Please press enter..."<<endl;
	cin >> in;
	cout << "You are very nice."<< in <<endl;
}

double getChi(double q, double u)
{
    double res = atan(u/q)/2.;
    if(q < 0.) res += PI/2.;
    else if (u < 0.) res += PI;
    return res;
}

double simpsint(int N,double *dg, double *di, double *dx,float a,float b);
double overplotPolarization(TH2F *I,TH2F *Q,TH2F *U);
void smooth(TH2F *I);
void smooth2(TH2F *J);


void plot(int nConfig,int nBin,double pmax)
{
// input parameters:
// nConfig: gives the spacetime we simulated
// nBin: gives the # of events we want to use
// pmax: gives the scale of the image in gravitational radii. pmax=20 gives good results.
    
// information about the file names of the simulated event samples
    
    char outNames[90][90][90]=
    {

        
         {"/Volumes/Work Drive/alltransfer/pol_agn_basic_sphere/5/out02_exspkn5_1.dat","/Volumes/Work Drive/alltransfer/pol_agn_basic_sphere/5/out02_exspkn5_2.dat"}, //9
        
        
        //  {"/Volumes/Work Drive/alltransfer/pol_agn_basic_sphere/out02_spT50t3.dat","/Volumes/Work Drive/alltransfer/pol_agn_basic_sphere/out02_spT50t3_2.dat"}, //4
        
        
        
    };

// how many data files should we read in (depends on how many events we simulated)

    cout <<"nConfig "<<nConfig<<endl;
    int outNum[20]={2,2};
    for (int i=0;i<10;i++)
        cout <<outNames[i][0]<<endl;
        
    gROOT->Reset();
    
// define a funky palette so that the images look nice
    Int_t MyPalette[100];
    Double_t r[]    = {0., 0.0, 1.0, 1.0, 1.0};
    Double_t g[]    = {0., 0.0, 0.0, 1.0, 1.0};
    Double_t b[]    = {0., 1.0, 0.0, 0.0, 1.0};
    Double_t stop[] = {0., .25, .50, .75, 1.0};
    Int_t FI = TColor::CreateGradientColorTable(5, stop, r, g, b, 100);
    for (int i=0;i<100;i++) MyPalette[i] = FI+i;
    
    gStyle->SetPalette(100, MyPalette);
    
// load the information about the accretion disk
    Disk *disk_p=new Disk;
    
    loadConfig(nConfig,disk_p);

// now we start reading in the information
    cout <<endl<<"Reading "<<nBin<<" events."<<endl;

    PhotonRed *gammaRed_p = new PhotonRed;
    int size = sizeof(PhotonRed);
    
    cout <<"Reading data in chunks of "<<size<<" units."<<endl;
    int total=0;
    
// we want to do the analysis for thetaBins different inclinations:
#define thetaBins 3
    double thetaCd[thetaBins]={45,60,75};
    double kmax = pmax/disk_p->rmax;
    
    cout <<"kmax "<<kmax<<endl;
    
// define dtheta window for the arrival position
#define dTheta 4.
#define nset 11
// save the results into a root file.
    char fileName[50];
    sprintf(fileName ,"out%d.root",(int)nConfig);

    TFile out(fileName,"recreate");
	out.cd();
    
    
    disk_p->rCorona_ends=6; //in r_g unit //BB
    disk_p->rCorona_start=disk_p->ISCO;

    // define 2-D histograms for holding the intensity maps
    TCanvas *cImage[thetaBins];
    TCanvas *cImageN[thetaBins];
    TH2F  *map[thetaBins],*mapN[thetaBins];
    TH2F  *mapI [thetaBins],*mapQ [thetaBins],*mapU [thetaBins];
    TH2F  *mapIN[thetaBins],*mapQN[thetaBins],*mapUN[thetaBins];
    
    
    // define 1-D histograms for other interesting results
    TCanvas *ChistosAll[thetaBins], *CHistosemisivityweight[thetaBins], *cHistosFlux[thetaBins], *Chistospol[thetaBins], *ChistospolAll[thetaBins], *Chistos[thetaBins];
    
    TH1F  *Flux[thetaBins],*CoronaFlux[thetaBins], *DiskFlux[thetaBins], *hI[nset][thetaBins], *hQ[nset][thetaBins], *hU[nset][thetaBins], *hPi[thetaBins], *hPDir[thetaBins], *hICorona[thetaBins], *hIDisk[thetaBins], *hQCorona[thetaBins], *hQDisk[thetaBins], *hUCorona[thetaBins], *hUDisk[thetaBins], *hPDirCorona[thetaBins], *hPDirDisk[thetaBins], *hPiCorona[thetaBins], *hPiDisk[thetaBins];
    TH1F *hIc[thetaBins],*hQc[thetaBins],*hUc[thetaBins];
    
    
    char name[100];
#define SEED  42
    
    // define binning for energy spectra histograms
#define NB 6
    double e1=0.1;
    double e2=100;
    double x1=log10(e1);
    double x2=log10(e2);
    
#define NBr 50
    
    
#define nMap 128
#define nPol 10
    
    
 
    
    for (int i=0; i<thetaBins;i++)
    {
        //Canveses
        
        sprintf(name ,"CImage%2d",(int)thetaCd[i]);
        cImage[i]  = new TCanvas(name,name,1+40*i,1+40*i,600,600);
        
        sprintf(name ,"CImageN%2d",(int)thetaCd[i]);
        cImageN[i] = new TCanvas(name,name,20+40*i,20+40*i,600,600);
        
        sprintf(name ,"CHistosAll%2d",(int)thetaCd[i]);
        ChistosAll[i] = new TCanvas(name,name,30+40*i,30+40*i,800,1000);
        ChistosAll[i]->SetLeftMargin(0.15); ChistosAll[i]->SetRightMargin(0.15);
    
        
        sprintf(name ,"cHistosFlux%2d",(int)thetaCd[i]);
        cHistosFlux[i] = new TCanvas(name,name,30+40*i,30+40*i,800,1000);
        cHistosFlux[i]->SetLeftMargin(0.15); cHistosFlux[i]->SetRightMargin(0.15);
        
        
        sprintf(name ,"Chistospol%2d",(int)thetaCd[i]);
        Chistospol[i] = new TCanvas(name,name,30+40*i,30+40*i,800,1000);
        Chistospol[i]->SetLeftMargin(0.15); Chistospol[i]->SetRightMargin(0.15);
        Chistospol[i]->Divide(1,2);
        
        sprintf(name ,"ChistospolAll%2d",(int)thetaCd[i]);
        ChistospolAll[i] = new TCanvas(name,name,30+40*i,30+40*i,800,1000);
        ChistospolAll[i]->SetLeftMargin(0.15); ChistospolAll[i]->SetRightMargin(0.15);
        ChistospolAll[i]->Divide(1,2);
        
        
        //Histograms
        for (int j=0; j<nset; j++) {
            sprintf(name ,"hI%d%d",(int)j,(int)thetaCd[i]);
            hI[j][i]   = new TH1F(name,";log(Energy); I",NB,x1,x2);
            sprintf(name ,"hQ%d%d",(int)j,(int)thetaCd[i]);
            hQ[j][i]   = new TH1F(name,";log(Energy); Q",NB,x1,x2);
            sprintf(name ,"hU%d%d",(int)j,(int)thetaCd[i]);
            hU[j][i]   = new TH1F(name,";log(Energy); U",NB,x1,x2);
          
            
        }
        
        sprintf(name ,"Flux%d",(int)thetaCd[i]);
        Flux[i]   = new TH1F(name,";log(Energy); E^2* dE/dN",NB,x1,x2);
        
        sprintf(name ,"hIc%d",(int)thetaCd[i]);
        hIc[i]   = new TH1F(name,";log(Energy); E^2* dE/dN",NB,x1,x2);
        
        sprintf(name ,"hQc%d",(int)thetaCd[i]);
        hQc[i]   = new TH1F(name,";log(Energy); E^2* dE/dN",NB,x1,x2);
        
        sprintf(name ,"hUc%d",(int)thetaCd[i]);
        hUc[i]   = new TH1F(name,";log(Energy); E^2* dE/dN",NB,x1,x2);
        
        sprintf(name ,"CoronaFlux%d",(int)thetaCd[i]);
        CoronaFlux[i]   = new TH1F(name,";log(Energy); E^2* dE/dN",NB,x1,x2);
        
        sprintf(name ,"DiskFlux%d",(int)thetaCd[i]);
        DiskFlux[i]   = new TH1F(name,";log(Energy); E^2* dE/dN",NB,x1,x2);
        
        
        sprintf(name ,"hICorona%d",(int)thetaCd[i]);
        hICorona[i]   = new TH1F(name,";log(Energy); I",NB,x1,x2);
        sprintf(name ,"hIDisk%d",(int)thetaCd[i]);
        hIDisk[i]   = new TH1F(name,";log(Energy); I",NB,x1,x2);

        
        
        sprintf(name ,"hQCorona%d",(int)thetaCd[i]);
        hQCorona[i]   = new TH1F(name,";log(Energy); Q",NB,x1,x2);
        sprintf(name ,"hQDisk%d",(int)thetaCd[i]);
        hQDisk[i]   = new TH1F(name,";log(Energy); Q",NB,x1,x2);
        
        
        sprintf(name ,"hUCorona%d",(int)thetaCd[i]);
        hUCorona[i]   = new TH1F(name,";log(Energy); U",NB,x1,x2);
        sprintf(name ,"hUDisk%d",(int)thetaCd[i]);
        hUDisk[i]   = new TH1F(name,";log(Energy); U",NB,x1,x2);
        
        sprintf(name ,"hPDir%d",(int)thetaCd[i]);
        hPDir[i]   = new TH1F(name,";log(Energy); Pol. Deg.",NB,x1,x2);
        sprintf(name ,"hPDirCorona%d",(int)thetaCd[i]);
        hPDirCorona[i]   = new TH1F(name,";log(Energy); Pol. Deg.",NB,x1,x2);
        sprintf(name ,"hPDirDisk%d",(int)thetaCd[i]);
        hPDirDisk[i]   = new TH1F(name,";log(Energy); Pol. Deg.",NB,x1,x2);
        
        sprintf(name ,"hPi%d",(int)thetaCd[i]);
        hPi[i]   = new TH1F(name,";log(Energy); Pol. Frac.",NB,x1,x2);
        sprintf(name ,"hPiCorona%d",(int)thetaCd[i]);
        hPiCorona[i]   = new TH1F(name,";log(Energy); Pol. Frac.",NB,x1,x2);
        sprintf(name ,"hPiDisk%d",(int)thetaCd[i]);
        hPiDisk[i]   = new TH1F(name,";log(Energy); Pol. Frac.",NB,x1,x2);
        
      
        
        // 2-D Maps
        
        // intensities
        
        sprintf(name ,"Map%2d",(int)thetaCd[i]);
        map[i]  = new TH2F(name,";-k_{r}^{#hat{#phi}}/k_{r}^{#hat{t}} [r_{g}/r_{r}];k_{r}^{#hat{#theta}}/k_{r}^{#hat{t}} [r_{g}/r_{r}]",nMap,-pmax,pmax,nMap,-pmax,pmax);
        
        sprintf(name ,"MapN%2d",(int)thetaCd[i]);
        mapN[i] = new TH2F(name,";-k_{r}^{#hat{#phi}}/k_{r}^{#hat{t}} [r_{g}/r_{r}];k_{r}^{#hat{#theta}}/k_{r}^{#hat{t}} [r_{g}/r_{r}]",nMap,-pmax,pmax,nMap,-pmax,pmax);
        map[i]->SetStats(kFALSE);
        mapN[i]->SetStats(kFALSE);
        
        // Stokes Parameters
        
        sprintf(name ,"IMap%2d",(int)thetaCd[i]);
        mapI[i]  = new TH2F(name,";-k_{r}^{#hat{#phi}}/k_{r}^{#hat{t}} [r_{g}/r_{r}];k_{r}^{#hat{#theta}}/k_{r}^{#hat{t}} [r_{g}/r_{r}]",nPol,-pmax,pmax,nPol,-pmax,pmax);
        
        sprintf(name ,"QMap%2d",(int)thetaCd[i]);
        mapQ[i]  = new TH2F(name,";-k_{r}^{#hat{#phi}}/k_{r}^{#hat{t}} [r_{g}/r_{r}];k_{r}^{#hat{#theta}}/k_{r}^{#hat{t}} [r_{g}/r_{r}]",nPol,-pmax,pmax,nPol,-pmax,pmax);
        
        sprintf(name ,"UMap%2d",(int)thetaCd[i]);
        mapU[i]  = new TH2F(name,";-k_{r}^{#hat{#phi}}/k_{r}^{#hat{t}} [r_{g}/r_{r}];k_{r}^{#hat{#theta}}/k_{r}^{#hat{t}} [r_{g}/r_{r}]",nPol,-pmax,pmax,nPol,-pmax,pmax);
        
        sprintf(name ,"IMapN%2d",(int)thetaCd[i]);
        mapIN[i]  = new TH2F(name,";-k_{r}^{#hat{#phi}}/k_{r}^{#hat{t}} [r_{g}/r_{r}];k_{r}^{#hat{#theta}}/k_{r}^{#hat{t}} [r_{g}/r_{r}]",nPol,-pmax,pmax,nPol,-pmax,pmax);
        
        sprintf(name ,"QMapN%2d",(int)thetaCd[i]);
        mapQN[i]  = new TH2F(name,";-k_{r}^{#hat{#phi}}/k_{r}^{#hat{t}} [r_{g}/r_{r}];k_{r}^{#hat{#theta}}/k_{r}^{#hat{t}} [r_{g}/r_{r}]",nPol,-pmax,pmax,nPol,-pmax,pmax);
        
        sprintf(name ,"UMapN%2d",(int)thetaCd[i]);
        mapUN[i]  = new TH2F(name,";-k_{r}^{#hat{#phi}}/k_{r}^{#hat{t}} [r_{g}/r_{r}];k_{r}^{#hat{#theta}}/k_{r}^{#hat{t}} [r_{g}/r_{r}]",nPol,-pmax,pmax,nPol,-pmax,pmax);
        
        
       
        
    }

// we also store the events in a big root tree,
// so that we can look at them with tree->Draw()
// using different cuts.
    
    double x,y,z,kr,kt,kp,kx,ky,pol;
    
    TTree *tree = new TTree("T","data"); 
	//tree->Branch("r0",&(gammaRed_p->r0),"r0/F");
	tree->Branch("T",&(gammaRed_p->T),"T/F");
	tree->Branch("u0PF",&(gammaRed_p->u0PF[0]),"u0PF[4]/F");
	tree->Branch("E",&(gammaRed_p->E),"E/F");
	tree->Branch("nS",&(gammaRed_p->nScatter),"nS/I");
    
	tree->Branch("x",&x,"x/F");
	//tree->Branch("y",&y,"y/F");
	tree->Branch("z",&z,"z/F");
	tree->Branch("kr",&kr,"kr/F");
	tree->Branch("kt",&kt,"kt/F");
	//tree->Branch("kp",&kp,"kp/F");
	tree->Branch("kx",&kt,"kx/F");
	tree->Branch("ky",&kp,"ky/F");
    
    
//  tree->Branch("xSC",&(gammaRed_p->xSC[0]),"xSC[4]/D");
//	tree->Branch("uSC",&(gammaRed_p->uSC[0]),"uSC[4]/D");
    
	
    tree->Branch("xBL",&(gammaRed_p->xBL[0]),"xBL[4]/F");
	tree->Branch("uCS",&(gammaRed_p->uCS[0]),"uCS[4]/F");
	//tree->Branch("fCS",&(gammaRed_p->fCS[0]),"fCS[4]/F");
	//tree->Branch("NuCS",&(gammaRed_p->NuCS),"NuCS/F");
	//tree->Branch("NfCS",&(gammaRed_p->NfCS),"NfCS/F");
	tree->Branch("stokes",&(gammaRed_p->stokes[0]),"stokes[4]/F");
	tree->Branch("pol",&pol,"pol/F");
    
    
    int numb=0;

    int seti;
//
//  ** Now we do the real work **
//  we loop over all files, and read in the information, and fill the histograms and maps
//
    
    for (int fnum=0;fnum<outNum[nConfig];fnum++)
    {
        ifstream i_file(outNames[nConfig][fnum],ios::in | ios::binary);   //open file for writing
        if ( !i_file.is_open()) {cout << "couldn't open output file "<< outNames[nConfig][fnum] << endl; return;}
        
        cout <<endl<<endl<<">>> Opened "<<outNames[nConfig][fnum]<<endl<<endl;
        
        while ((total++<nBin)&&(i_file.read((char*)gammaRed_p,size)))
        {
            
            //print the detailed event information for the first 5 events
            if ((total<5)||(total%500000==0)) cout <<total<<" "<<" r0 "<<gammaRed_p->r0<<" T "<<gammaRed_p->T<<" E "<<gammaRed_p->E <<" nS "<<gammaRed_p->nScatter<<" xBL "<<gammaRed_p->xBL[0]<<" "<<gammaRed_p->xBL[1]<<" "<<gammaRed_p->xBL[2]<<" "<<gammaRed_p->xBL[3]<<" "<<" uCS "<<gammaRed_p->uCS[0]<<" "<<gammaRed_p->uCS[1]<<" "<<gammaRed_p->uCS[2]<<" "<<gammaRed_p->uCS[3]<<" "<<endl;
            
            seti=total%nset;
            
            // we only simulated the emission into the upper hemisphere;
            // we should have simulated the emission into the lower hemisphere also
            // however, because of the symmetry of the accretion disk, we can
            // simply use the photons arriving below the accretion disk and
            // "flip them up". In this way, we get all the events that we would have
            // gotten, had we simulated the disk emission into the lower hemisphere.
            if (gammaRed_p->xBL[2]>PI/2.) {
                gammaRed_p->xBL[2] = PI/2. - (gammaRed_p->xBL[2]-PI/2.);
                gammaRed_p->uCS[2] = -gammaRed_p->uCS[2];
                gammaRed_p->stokes[2]=-gammaRed_p->stokes[2];
            }
            
            // transform from spherical coordinates of the final photon position into x-y-z coordinates
            x=gammaRed_p->xBL[1]*sin(gammaRed_p->xBL[2])*cos(gammaRed_p->xBL[3]);
            y=gammaRed_p->xBL[1]*sin(gammaRed_p->xBL[2])*sin(gammaRed_p->xBL[3]);
            z=gammaRed_p->xBL[1]*cos(gammaRed_p->xBL[2]);
            
            // get the wave vector of the photon in the coordinate stationary frame
            kr=gammaRed_p->uCS[1];
            kt=gammaRed_p->uCS[2];
            kp=gammaRed_p->uCS[3];
            
            // get the normalized wave vector components in x and y direction.
            kx=-gammaRed_p->uCS[3]/fabs(gammaRed_p->uCS[0]);
            ky= gammaRed_p->uCS[2]/fabs(gammaRed_p->uCS[0]);

            // get the polarization fraction if the photon was successfully tracked
            if (gammaRed_p->stokes[0]>0.)
                pol=sqrt( sq(gammaRed_p->stokes[1])+sq(gammaRed_p->stokes[2]) )/ gammaRed_p->stokes[0];
            else
                pol=0.;
            
#ifdef TreeFill 
            if ((numb<1000000)||((numb++)%10==0)) tree->Fill(); // only fill the tree with 1 in 10 events
#endif            
           
           
        
            if ((gammaRed_p->xBL[1]>=disk_p->rmax))
            {
                for (int i=0;i<thetaBins;i++)
                    if ((gammaRed_p->xBL[2]>(thetaCd[i]-dTheta)/RtD)&&
                        (gammaRed_p->xBL[2]<(thetaCd[i]+dTheta)/RtD))
                    {
                        
                       
                        // if photon scattered in the corona

                        if ((gammaRed_p->nScatter_Corona!=0)) {
                           
                            CoronaFlux[i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[0]);
                            
                            hICorona[i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[0]);
                            hQCorona[i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[1]);
                            hUCorona[i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[2]);
                            
                        }
                        
                        if (gammaRed_p->uCS[0]>e1 && gammaRed_p->uCS[0]<e2) {
                            

                            Flux[i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[0]);
                            hI[seti][i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[0]);
                            hQ[seti][i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[1]);
                            hU[seti][i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[2]);

                        }
                        

                      // photons did not scatter in the corona
                           
                        if ((gammaRed_p->nScatter_Corona==0)) {
                            
                            
                            DiskFlux[i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[0]);
                           
                            hIDisk[i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[0]);
                            hQDisk[i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[1]);
                            hUDisk[i]->Fill(log10(gammaRed_p->uCS[0]),gammaRed_p->uCS[0]*gammaRed_p->stokes[2]);
                         
                            
                        }
                        
                        
                        
                     
                        //map for
                        
                       
                        if ((gammaRed_p->nScatter_Corona==1)&&(gammaRed_p->nScatter==0)&&(gammaRed_p->uCS[0]>e1 && gammaRed_p->uCS[0]<e2)) {
                           
                                map[i]->Fill(kx *pmax/kmax,ky *pmax/kmax,gammaRed_p->stokes[0]/100);
                                
                                mapI[i]->Fill(kx*pmax/kmax,ky*pmax/kmax,gammaRed_p->u0PF[0]*gammaRed_p->stokes[0]);
                                mapQ[i]->Fill(kx*pmax/kmax,ky*pmax/kmax,gammaRed_p->u0PF[0]*gammaRed_p->stokes[1]);
                                mapU[i]->Fill(kx*pmax/kmax,ky*pmax/kmax,gammaRed_p->u0PF[0]*gammaRed_p->stokes[2]); //
                           
                        }
                        
                        if ((gammaRed_p->uCS[0]>e1 && gammaRed_p->uCS[0]<e2)) {

                            mapN[i]->Fill(kx *pmax/kmax,ky *pmax/kmax,gammaRed_p->stokes[0]/100.);
                            
                            mapIN[i]->Fill(kx*pmax/kmax,ky*pmax/kmax,gammaRed_p->u0PF[0]*gammaRed_p->stokes[0]);
                            mapQN[i]->Fill(kx*pmax/kmax,ky*pmax/kmax,gammaRed_p->u0PF[0]*gammaRed_p->stokes[1]);
                            mapUN[i]->Fill(kx*pmax/kmax,ky*pmax/kmax,gammaRed_p->u0PF[0]*gammaRed_p->stokes[2]);
                        }
                        
                        
                        
                    }
            }
        }
        
        i_file.close();  
    }

//
// ** Done with the real work **
//

    
    

for (int i=0;i<thetaBins;i++){
    
    double sumI[NB];
    double sumQ[NB];
    double sumU[NB];
    
    for (int iE=0; iE<NB; iE++) {
        sumI[iE]=0.;
        sumQ[iE]=0.;
        sumU[iE]=0.;
    }

    double sigmafrac[NB],sigmaang[NB];
    double sigmaI[NB],sigmaQ[NB],sigmaU[NB];
    for (int j=0; j<nset; j++) {
        
        //smothing

        for (int iE=1; iE<(NB-1); iE++) {
            
            hI[j][i]->SetBinContent(iE+1,(hI[j][i]->GetBinContent(iE)+hI[j][i]->GetBinContent(iE+1)+hI[j][i]->GetBinContent(iE+2))/3);
            hQ[j][i]->SetBinContent(iE+1,(hQ[j][i]->GetBinContent(iE)+hQ[j][i]->GetBinContent(iE+1)+hQ[j][i]->GetBinContent(iE+2))/3);
            hU[j][i]->SetBinContent(iE+1,(hU[j][i]->GetBinContent(iE)+hU[j][i]->GetBinContent(iE+1)+hU[j][i]->GetBinContent(iE+2))/3);
            
        }
        hI[j][i]->SetBinContent(1,(hI[j][i]->GetBinContent(1)+hI[j][i]->GetBinContent(2))/2);
        hI[j][i]->SetBinContent(NB,(hI[j][i]->GetBinContent(NB)+hI[j][i]->GetBinContent(NB-1))/2);
        hQ[j][i]->SetBinContent(1,(hQ[j][i]->GetBinContent(1)+hQ[j][i]->GetBinContent(2))/2);
        hQ[j][i]->SetBinContent(NB,(hQ[j][i]->GetBinContent(NB)+hQ[j][i]->GetBinContent(NB-1))/2);
        hU[j][i]->SetBinContent(1,(hU[j][i]->GetBinContent(1)+hU[j][i]->GetBinContent(2))/2);
        hU[j][i]->SetBinContent(NB,(hU[j][i]->GetBinContent(NB)+hU[j][i]->GetBinContent(NB-1))/2);
        
        
        for (int iE=0; iE<NB; iE++) {
            sumI[iE]+=hI[j][i]->GetBinContent(iE+1);
            sumQ[iE]+=hQ[j][i]->GetBinContent(iE+1);
            sumU[iE]+=hU[j][i]->GetBinContent(iE+1);
            
        }
        
       

    }

    for (int iE=0; iE<NB; iE++) {
        hIc[i]->SetBinContent(iE+1,sumI[iE]/(double)nset);
        hQc[i]->SetBinContent(iE+1,sumQ[iE]/(double)nset);
        hUc[i]->SetBinContent(iE+1,sumU[iE]/(double)nset);
        
    }
    
    for (int iE=0;iE<NB;iE++)
    {
        
        sigmafrac[iE]=0;
        sigmaang[iE]=0;
        
        double res=0.;
        if (hIc[i]->GetBinContent(iE+1)>0.)
            res= sqrt( sq(hQc[i]->GetBinContent(iE+1))+sq(hUc[i]->GetBinContent(iE+1)) )/ hIc[i]->GetBinContent(iE+1);
        
        hPi[i]->SetBinContent(iE+1,res*100);
        
        if(hQc[i]->GetBinContent(iE+1)!=0)
            hPDir[i]->SetBinContent(iE+1,getChi(hQc[i]->GetBinContent(iE+1),hUc[i]->GetBinContent(iE+1))*180./PI);
        else
            hPDir[i]->SetBinContent(iE+1,0);
    }
    
    
    
    
    for (int iE=0; iE<NB; iE++) {
        double sum1=0;
        double sum2=0;
        double sum3=0;

        for (int j=0; j<nset; j++) {
            sum1+=sq(hI[j][i]->GetBinContent(iE+1)-hIc[i]->GetBinContent(iE+1));
            sum2+=sq(hQ[j][i]->GetBinContent(iE+1)-hQc[i]->GetBinContent(iE+1));
            sum3+=sq(hU[j][i]->GetBinContent(iE+1)-hUc[i]->GetBinContent(iE+1));

        }
        sigmaI[iE]=sqrt(sum1/(nset*(nset-1)));
        sigmaQ[iE]=sqrt(sum2/(nset*(nset-1)));
        sigmaU[iE]=sqrt(sum3/(nset*(nset-1)));
        
        double megfac= (sq(hQc[i]->GetBinContent(iE+1)*sigmaQ[iE])+sq(hUc[i]->GetBinContent(iE+1)*sigmaU[iE]))/(sq(hIc[i]->GetBinContent(iE+1))*(sq(hQc[i]->GetBinContent(iE+1))+sq(hUc[i]->GetBinContent(iE+1)))) + ((sq(hQc[i]->GetBinContent(iE+1))+sq(hUc[i]->GetBinContent(iE+1))))*sq(sigmaI[iE])/pow(hIc[i]->GetBinContent(iE+1),4);
        sigmafrac[iE]=sqrt(megfac)*100;
        
        double magang=sq(sigmaU[iE]*hQc[i]->GetBinContent(iE+1)/(sq(hQc[i]->GetBinContent(iE+1))+sq(hUc[i]->GetBinContent(iE+1))))+sq(sigmaQ[iE]*hUc[i]->GetBinContent(iE+1)/(sq(hQc[i]->GetBinContent(iE+1))+sq(hUc[i]->GetBinContent(iE+1))));
        sigmaang[iE]=0.5*sqrt(magang)*180./PI;
        
        cout<<"sigma for bin "<<iE<<" is "<<sigmafrac[iE]<<" "<<sigmaang[iE]<<endl;


    }
    
    
    
    

    hPi[i]->SetStats(kFALSE);
    hPDir[i]->SetStats(kFALSE);
    
    Chistospol[i]->cd(1);  hPi[i]->DrawCopy();
    Chistospol[i]->cd(2);  hPDir[i]->DrawCopy();
    Chistospol[i]->SaveAs(".png");
    
    for (int iE=0;iE<NB;iE++)
    {
        double res=0.;
        if (hICorona[i]->GetBinContent(iE+1)>0.)
            res= sqrt( sq(hQCorona[i]->GetBinContent(iE+1))+sq(hUCorona[i]->GetBinContent(iE+1)) )/ hICorona[i]->GetBinContent(iE+1);
        
        hPiCorona[i]->SetBinContent(iE+1,res*100);
        if(hQCorona[i]->GetBinContent(iE+1)!=0)
            hPDirCorona[i]->SetBinContent(iE+1,getChi(hQCorona[i]->GetBinContent(iE+1),hUCorona[i]->GetBinContent(iE+1))*180./PI);
        else
            hPDirCorona[i]->SetBinContent(iE+1,0);
    }
    for (int iE=0;iE<NB;iE++)
    {
        double res=0.;
        if (hIDisk[i]->GetBinContent(iE+1)>0.)
            res= sqrt( sq(hQDisk[i]->GetBinContent(iE+1))+sq(hUDisk[i]->GetBinContent(iE+1)) )/ hIDisk[i]->GetBinContent(iE+1);
        
        hPiDisk[i]->SetBinContent(iE+1,res*100);
        if(hQDisk[i]->GetBinContent(iE+1)!=0)
            hPDirDisk[i]->SetBinContent(iE+1,getChi(hQDisk[i]->GetBinContent(iE+1),hUDisk[i]->GetBinContent(iE+1))*180./PI);
        else
            hPDirDisk[i]->SetBinContent(iE+1,0);
    }
    
    hPiDisk[i]->SetStats(kFALSE);
    hPDirDisk[i]->SetStats(kFALSE);
    hPiCorona[i]->SetStats(kFALSE);
    hPDirCorona[i]->SetStats(kFALSE);
    
    ChistospolAll[i]->cd(1);
    hPi[i]->SetLineStyle(1); hPi[i]->GetYaxis()->SetRangeUser(0,7); hPi[i]->DrawCopy();
    hPiDisk[i]->SetLineStyle(2);  hPiDisk[i]->DrawCopy("SAME");
    hPiCorona[i]->SetLineStyle(3); hPiCorona[i]->DrawCopy("SAME");
    
    ChistospolAll[i]->cd(2);
    hPDir[i]->SetLineStyle(1); hPDir[i]->DrawCopy();
    hPDirDisk[i]->SetLineStyle(2);  hPDirDisk[i]->DrawCopy("SAME");
    hPDirCorona[i]->SetLineStyle(3); hPDirCorona[i]->DrawCopy("SAME");
    
    ChistospolAll[i]->SaveAs(".png");
    ChistospolAll[i]->Close();
    
    
    Flux[i]->SetStats(kFALSE);
    DiskFlux[i]->SetStats(kFALSE);
    CoronaFlux[i]->SetStats(kFALSE);
    
    
    ChistosAll[i]->cd();  gPad->SetLogy();
    Flux[i]->SetLineStyle(1); Flux[i]->DrawCopy();
    DiskFlux[i]->SetLineStyle(2);  DiskFlux[i]->DrawCopy("SAME");
    CoronaFlux[i]->SetLineStyle(3); CoronaFlux[i]->DrawCopy("SAME");
    
    ChistosAll[i]->SaveAs(".png");
    ChistosAll[i]->Close();
    
    TF1 *func = new TF1("fit",fitf,log10(e1),log10(e2),2);
    func->SetParameters(1E+9,-1.);
    func->SetParNames("Constant","Gamma");
    
    Flux[i]->Fit("fit","+Q","same",-0.1,1.);
    
    
    double p1 = func->GetParameter(1);
    cout <<"Fitted  Weighted  E^2 dN/dE index: "<<p1<<endl;
    
    
    cHistosFlux[i]->cd();  gPad->SetLogy();  Flux[i]->DrawCopy("");
    cHistosFlux[i]->SaveAs(".png");
    cHistosFlux[i]->Close();

    
    char fileNameOut2[100];
    ofstream print;
    sprintf(fileNameOut2, "fluxpol%d.txt", (int) i);
    print.open(fileNameOut2, ios::out);
    cout<<"print it"<<fileNameOut2<<endl;
    for (int iE=0; iE<NB; iE++) {
        
        print<<hPi[i]->GetXaxis()->GetBinCenter(iE+1)<<" "<<Flux[i]->GetBinContent(iE+1)<<" "<<hPi[i]->GetBinContent(iE+1)<<" "<<hPDir[i]->GetBinContent(iE+1)<<" "<<sigmafrac[iE]<<" "<<sigmaang[iE]<<endl;
    
    }
    print.close();
    
    
    
}

    //plot maps and histos
    
    for (int i=0;i<thetaBins;i++){
        cout <<"Drawing i "<<i<<" theta " <<thetaCd[i]<<endl;
        
        cImage[i]->cd();
        
        smooth(map[i]);
        map[i]->Smooth(1,"k5b");
        smooth2(map[i]);
        
        map[i]->SetMinimum(0.);
        map[i]->DrawCopy("colz");
        cout <<"Including Scattering, Theta-Bin "<<i<<" Max. Pol. "<<overplotPolarization(mapI[i],mapQ[i],mapU[i])<<endl;
        
        cImage[i]->Update();
        cImage[i]->SaveAs(".png");
        cImage[i]->Close();
        
        cImageN[i]->cd();
        
        smooth(mapN[i]);
        
        mapN[i]->Smooth(1,"k5b");
        smooth2(mapN[i]);
        mapN[i]->SetMinimum(0.);
        mapN[i]->DrawCopy("colz");
        cout <<"Excluding Scattering, Theta-Bin "<<i<<" Max. Pol. "<<overplotPolarization(mapIN[i],mapQN[i],mapUN[i])<<endl;
        
        cImageN[i]->Update();
        cImageN[i]->SaveAs(".png");
        cImageN[i]->Close();
        
    }

#ifdef TreeFill

    char fname[50];
    sprintf(fname,"track.dat");
    ifstream t_file(fname,ios::in);   //open file for writing

    if ( !t_file.is_open()) {cout << "couldn't open output file" << endl; exit(1);}

    int tnumber;
    double bl0,bl1,bl2,bl3;
    double k0,k1,k2,k3;
    
    TTree *track = new TTree("track","track"); 
	track->Branch("n",&tnumber,"n/I");
	track->Branch("bl0",&bl0,"bl0/D");
	track->Branch("bl1",&bl1,"bl1/D");
	track->Branch("bl2",&bl2,"bl2/D");
	track->Branch("bl3",&bl3,"bl3/D");
    
    track->Branch("x",&x,"x/D");
	track->Branch("y",&y,"y/D");
	track->Branch("z",&z,"z/D");
    
	track->Branch("k0",&k0,"k0/D");
	track->Branch("k1",&k1,"k1/D");
	track->Branch("k2",&k2,"k2/D");
	track->Branch("k3",&k3,"k3/D");
    
    while (t_file >>tnumber>>bl0>>bl1>>bl2>>bl3>>k0>>k1>>k2>>k3)
    {
        x=bl1*sin(bl2)*cos(bl3);
        y=bl1*sin(bl2)*sin(bl3);
        z=bl1*cos(bl2);
        track->Fill();
    };
    
    t_file.close();

#endif  
    
    out.Write();
	out.Close();

    cout <<"Done."<<endl;
}

/* custom functions */

void loadConfig(int nConfig, Disk *disk_p)
{
    char root[20][30]={
        "Parm_1_0p0_0",
        "Parm_1_0p5_0",    
        "Parm_1_0p9_0",
        "Parm_1_0p95_0",
        "Parm_1_0p5_0a25",
        "Parm_1_0p5_0a4",
        "Parm_1_0p5-2p5",
        "Parm_1_0p5+2p5",
        "Parm_1_0p5-5",
//        "Parm_1_0p5+5",    
        "Parm_1_0p5+6p32915",    
        "Parm_1_0p99-5", // 10 
        "Parm_1_0p99-2p5",
        "Parm_1_0p5-30p611"};
    
    char emui[6][30]={"lmi00","lmi03","lmi11","lmi22","lmi30","lmi33"};
    char eimu[6][30]={"lim00","lim03","lim11","lim22","lim30","lim33"};
    
    int a[4][4]={  
        0,-1,-1,1,
        -1,2,-1,-1,
        -1,-1,3,-1,
        4,-1,-1,5};
    
    
    char fname[50]; 
    int dummy;
    
    cout<<"Loading ";
    cout<<root[nConfig]<<endl;
    
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

void smooth(TH2F *J)
{
    double av;
    for (int i=1;i<=J->GetNbinsX();i++)
        for (int j=1;j<=J->GetNbinsY();j++)
        {
            av =  100.*J->GetBinContent(i,j);
            if (av>0.) J->SetBinContent(i,j,log10(av));
            else J->SetBinContent(i,j,0.001);
        }
}

void smooth2(TH2F *J)
{
    double av;
    for (int i=1;i<=J->GetNbinsX();i++)
        for (int j=1;j<=J->GetNbinsY();j++)
        {
            av =  J->GetBinContent(i,j);
            if ((av>0.)&&(av<10.)) J->SetBinContent(i,j,av);
            else J->SetBinContent(i,j,0.001);
        }
}

double overplotPolarization(TH2F *I,TH2F *Q,TH2F *U)
{
    double polDG[100][100];
    double polDir[100][100];
    
    double maxPol=-1.;
    
    double kmax=I->GetXaxis()->GetXmax();
    
    cout<<" i and j max "<<I->GetNbinsX()<<" "<<I->GetNbinsY()<<endl;
    
    for (int i=1;i<=I->GetNbinsX();i++)
        for (int j=1;j<=I->GetNbinsY();j++)
        {
            double denum=I->GetBinContent(i,j);
            if (denum>0.)
            {
                polDG[i][j]  = sqrt(Q->GetBinContent(i,j)*Q->GetBinContent(i,j)+U->GetBinContent(i,j)*U->GetBinContent(i,j))/denum;
                if(Q->GetBinContent(i,j)!=0)
                    polDir[i][j] = getChi(Q->GetBinContent(i,j),U->GetBinContent(i,j));
                else
                    polDir[i][j]=0;
                
               // cout<<"i , j "<<i<<" "<<j<<" pol angle "<<polDir[i][j]<<" pol fraction "<<polDG[i][j]<<" I,Q,U "<<I->GetBinContent(i,j)<<" "<<Q->GetBinContent(i,j)<<" "<<U->GetBinContent(i,j)<<endl;
                
                if (polDG[i][j]>maxPol) maxPol=polDG[i][j];
            }
            else
            {
                polDG[i][j]=0;
                polDir[i][j]=0.;
            }
        }
    
    int npol= I->GetNbinsX();
    
    double kx0 = -kmax;
    double ky0 = -kmax;
    double delta  = 2.*kmax/(double)npol;
    double maxL   = 0.98 * 2.*kmax/(double)npol / 2.;
    if (maxPol>0.3) maxPol=0.3;
    maxPol=0.35;
    maxPol=0.08;
    
    for (int i=1;i<=I->GetNbinsX();i++)
        for (int j=1;j<=I->GetNbinsY();j++)
        {
            double kxc = kx0+((double)i-0.5)*delta;
            double kyc = ky0+((double)j-0.5)*delta;
            double length=maxL * polDG[i][j]/maxPol;
            if (length>20.*maxL) length=20.*maxL;
            double dx = length * sin(polDir[i][j]);
            double dy = length * cos(polDir[i][j]);
            
           // cout<<"dx dy "<<dx<<" "<<dy<<" lenght "<<length<<endl;
            
            TLine *l = new TLine(kxc-dx,kyc-dy,kxc+dx,kyc+dy);
            l->SetLineWidth(3);
            l->SetLineColor(18);
            l->SetLineColor(1);
            
            l->Draw();
            
           
            
            //            delete l;
        }
    return maxPol;
}



double fitf(double *v, double *par)
{
    double fitval = par[0]*pow(pow(10.,v[0]),-par[1]);
    return fitval;
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

#define MAXDIMS 10001
double simpsint(int N,double *dg, double *di, double *dx,float a,float b)
/*     integrate "dg" from "a" to "b"
 dg(1)... dg(n) corresponds
 dg(a)... dg(b) equally spaced
 
 N has to be odd
 in di(x) int_a^x is stored
 in dx(i) the abscissa is stored */
{
    /* step width */
    float    delta;
    double   dresult;
    
    /* dummy Array */
    double f[MAXDIMS];
    
    /* Dummy Variables */
    int i,ind;
    
    delta   = ((double)(b-a))/((double)(N-1));
    
    if (N>MAXDIMS)
    {
        printf("MAXDIMS too small\n");
        exit(1);
    }
    if ((N%2)!=1)
    {
        b -= delta;
        N--;
    }
    
    /* Copy array */
    for (i=0;i<N;i++)
        f[i] = dg[i];
    
    dresult = (double) 0.;
    di[0]   = (double) 0.;
    dx[0]   = (double) a;
    
    for (i=0; i<(int)((N-1)/2); i++)
    {
        ind       = 2*i+1;
        dresult  += (delta/(double)3.)*(f[ind-1]+(double)4.*f[ind]+f[ind+1]);
        
        /* Store dresults in "di" */
        
        di[ind]     = (di[ind-1]+dresult)/(double)2.;
        di[ind+1] = dresult;
        
        /* Store dresults in "x" */
        dx[ind  ] = dx[ind-1]+delta;
        dx[ind+1] = dx[ind  ]+delta;
    }
    return dresult;
}




