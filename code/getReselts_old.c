#include <math.h>
#include <stdio.h>  
#include<stdlib.h>
#define NDATA 800
#define NTRAJ 50
#define NBLOCKS 1
#define CHOP 0
main()
{
//double edot = 5.0;
double Gdot = 0.0;
double chi = 1.0;
//double pref = ((4.0*chi) + ((1.0 - chi)*(1.0 - chi)))/(sqrt(chi) + (1.0 - chi));
//printf("prefactor = %f \n", pref);
int Nc = 1;
double vol = 2232.232033;
double A, xin, detEr, xinEr, phiEr;
A = 4.0*M_PI*Nc/vol;
double phi;
int i, j, k, ib;
int NDATA1 = NDATA - CHOP;
int block = NDATA/NBLOCKS; // size of the block
double time, det;
double Re2;
double Rg2;
double X;
double EtaShr;
double EtaExt;
double EtaMix;
double EtaMixNew;
double EtaMixnew;
double FirstNormalStressDiff;
double g11, g12, g13, g21, g22, g23, g31, g32, g33;
double EtaMixNewTemp;
double Re2mean = 0.0; 
double Rg2mean = 0.0;
double Re2traj[300], Rg2traj[300];
double Re2block = 0.0, Rg2block = 0.0; // block avergaing over the trajectory
double Re2bsq = 0.0, Rg2bsq = 0.0; 
double Xmean = 0.0;
double EtaShrmean = 0.0;
double EtaExtmean = 0.0;
double EtaMixmean = 0.0;
double g11mean = 0.0; double g12mean = 0.0; double g13mean = 0.0;
double g21mean = 0.0; double g22mean = 0.0; double g23mean = 0.0;
double g31mean = 0.0; double g32mean = 0.0; double g33mean = 0.0;
double EtaMixmeanNew = 0.0;
double NormalStressDiffmean = 0.0;
/*
double Re2M=0.0;
double Rg2M=0.0;
double EtaShrM=0.0;
double EtaExtM=0.0;
double EtaMixM=0.0;
*/
double Re2er=0.0;
double Rg2er=0.0;
double Re2trajER = 0.0, Rg2trajER = 0.0;
double Re2ber, Rg2ber;
df (i < CHOP) {
                fscanf(fp, " %lf %lf %lf %lf %lf %d \n",&d6,&d2,&d3,&d4,&d5,&d1);
                }
                if (i >= CHOP) {
                fscanf(fp, "%lf %lf %lf %lf %lf %d \n", &time, &Re2, &Rg2, &dummy2, &dummy3, &dummy1);
//                printf("after fscanf %lf \t %lf \n", Re2, Rg2);
//                ouble Xer = 0.0;
double EtaShrer=0.0;
double EtaExter=0.0;
double EtaMixer=0.0;
double NormalStressDiffer = 0.0;
double g11er = 0.0; double g12er = 0.0; double g13er = 0.0;
double g21er = 0.0; double g22er = 0.0; double g23er = 0.0;
double g31er = 0.0; double g32er = 0.0; double g33er = 0.0;
double EtaMixerNew=0.0;
int dummy1 = 0;
double dummy2=0.0;
double dummy3=0.0;
int d1 = 0;
double d2 = 0.0; 
double d3 = 0.0;
double d4 = 0.0;
double d5 = 0.0;
double d6 = 0.0;
double d7 = 0.0;
double d8 = 0.0;
double d9 = 0.0; double d10 = 0.0; double d11 = 0.0; double d12 = 0.0;
double d13 = 0.0; double d14 = 0.0; double d15 = 0.0; double d16 = 0.0;
double d17 = 0.0; double d18 = 0.0; 
const char baseProp[] = "traj_";
const char extProp[] = ".txt";
char filename [ FILENAME_MAX ];
for (k=0; k<NTRAJ; k++){
     Re2traj[k] = 0;
     Rg2traj[k] = 0;
   }
for (k=0;k<NTRAJ;k++) {
	FILE * fp;
		sprintf(filename, "%s%d%s", baseProp, k+101, extProp);
	//printf("filename = %s \n", filename);
	fp = fopen(filename,"r");
	j = 0;
	for (i=0;i<NDATA;i++) {
                
		if (i < CHOP) {
		fscanf(fp, " %lf %lf %lf %lf %lf %d \n",&d6,&d2,&d3,&d4,&d5,&d1);
		} 
		if (i >= CHOP) {
		fscanf(fp, "%lf %lf %lf %lf %lf %d \n", &time, &Re2, &Rg2, &dummy2, &dummy3, &dummy1);
//                printf("after fscanf %lf \t %lf \n", Re2, Rg2);
		//j = j + 1;
		Re2mean = Re2mean + Re2;
                Re2traj[k] += Re2; 
                
		Rg2mean = Rg2mean + Rg2;
                Rg2traj[k] += Rg2; 
               
                while (j < block) {
                  Re2block += Re2;
                  Rg2block += Rg2;
                  j = j+1;
                } 
                if ((i+1)%block == 0) {
                  Re2bsq += (Re2block/block)*(Re2block/block);  // block averages
                  Rg2bsq += (Rg2block/block)*(Rg2block/block);
                 
                  Re2block = 0;
                  Rg2block = 0;  
                  j = 0;
                }
                
		Xmean = Xmean + X;
                EtaMixnew = - ((sqrt(chi)*(EtaExt - EtaMix)+(1.0 - chi)*EtaShr)/(Nc*Gdot*(4.0*(chi)+(1.0 - chi)*(1.0-chi)))); 
		//EtaShrmean = EtaShrmean + EtaShr;
		//EtaExtmean = EtaExtmean + EtaExt;
		EtaMixmean = EtaMixmean + EtaMixnew;
                FirstNormalStressDiff =  -(EtaExt - EtaMix)/(Nc*Gdot);
                NormalStressDiffmean =  NormalStressDiffmean + FirstNormalStressDiff;
                g11mean = g11mean + g11; g12mean = g12mean + g12; g13mean = g13mean + g13;
		g21mean = g21mean + g21; g22mean = g22mean + g22; g23mean = g23mean + g23;
		g31mean = g31mean + g31; g32mean = g32mean + g32; g33mean = g33mean + g33;
		//EtaMixmeanNew += pref*(((edot*edot*EtaExt)+(gdot*gdot*EtaShr))/((4.0*edot*edot)+(gdot*gdot)));
              	}
	}
	fclose(fp);
}
// Calculating mean values
Re2mean = Re2mean/(NTRAJ*NDATA1);
Rg2mean = Rg2mean/(NTRAJ*NDATA1);
Xmean = Xmean/(NTRAJ*NDATA1);
// calculating the average over each trajectory
for (k=0; k<NTRAJ; k++){
     Re2traj[k] = Re2traj[k]/NDATA1;
     Rg2traj[k] = Rg2traj[k]/NDATA1;
   }
//EtaShrmean = EtaShrmean/(NTRAJ*NDATA1);
//EtaExtmean = EtaExtmean/(NTRAJ*NDATA1);
EtaMixmean = EtaMixmean/(NTRAJ*NDATA1);
NormalStressDiffmean =  NormalStressDiffmean /(NTRAJ*NDATA1);
g11mean = g11mean/(NTRAJ*NDATA1); g12mean = g12mean/(NTRAJ*NDATA1); g13mean = g13mean/(NTRAJ*NDATA1);
g21mean = g21mean/(NTRAJ*NDATA1); g22mean = g22mean/(NTRAJ*NDATA1); g23mean = g23mean/(NTRAJ*NDATA1);
g31mean = g31mean/(NTRAJ*NDATA1); g32mean = g32mean/(NTRAJ*NDATA1); g33mean = g33mean/(NTRAJ*NDATA1);
//EtaMixmeanNew = EtaMixmeanNew/(NTRAJ*NDATA1);
// Calculating standard errors
for (k=0;k<NTRAJ;k++) {
        FILE * fp;
        sprintf(filename, "%s%d%s", baseProp, k+101, extProp);
        //printf("filename = %s \n", filename);
        //printf("Rg2erchk = %lf \n", Rg2er);
        fp = fopen(filename,"r");
        j = 0;
        for (i=0;i<NDATA;i++) {
                if (i < CHOP) {
		 fscanf(fp, " %lf %lf %lf %lf %lf %d \n",&d6,&d2,&d3,&d4,&d5,&d1);
                }
                if (i >= CHOP) {
                fscanf(fp, "%lf %lf %lf %lf %lf %d \n", &time, &Re2, &Rg2, &dummy2, &dummy3, &dummy1);


                j = j + 1;
                Re2er = Re2er + ((Re2 - Re2mean)*(Re2 - Re2mean));
		Rg2er = Rg2er + ((Rg2 - Rg2mean)*(Rg2 - Rg2mean));
		Xer = Xer + ((X - Xmean)*(X - Xmean));
                EtaMixnew = - ((sqrt(chi)*(EtaExt - EtaMix)+(1.0 - chi)*EtaShr)/(Nc*Gdot*(4.0*(chi)+(1.0 - chi)*(1.0-chi))));
                FirstNormalStressDiff =  -(EtaExt - EtaMix)/(Nc*Gdot);
		//EtaShrer = EtaShrer + ((EtaShr - EtaShrmean)*(EtaShr - EtaShrmean));
		//EtaExter = EtaExter + ((EtaExt - EtaExtmean)*(EtaExt - EtaExtmean));
		EtaMixer = EtaMixer + ((EtaMixnew - EtaMixmean)*(EtaMixnew - EtaMixmean));
                NormalStressDiffer = NormalStressDiffer + ((FirstNormalStressDiff - NormalStressDiffmean)*(FirstNormalStressDiff - NormalStressDiffmean));
		g11er = g11er + ((g11 - g11mean)*(g11 - g11mean));
		g12er = g12er + ((g12 - g12mean)*(g12 - g12mean));
		g13er = g13er + ((g13 - g13mean)*(g13 - g13mean));
                g21er = g21er + ((g21 - g21mean)*(g21 - g21mean));
                g22er = g22er + ((g22 - g22mean)*(g22 - g22mean));
                g23er = g23er + ((g23 - g23mean)*(g23 - g23mean));
                g31er = g31er + ((g31 - g31mean)*(g31 - g31mean));     
                g32er = g32er + ((g32 - g32mean)*(g32 - g32mean));
                g33er = g33er + ((g33 - g33mean)*(g33 - g33mean));
		//EtaMixNewTemp = pref*(((edot*edot*EtaExt)+(gdot*gdot*EtaShr))/((4.0*edot*edot)+(gdot*gdot)));
 		//EtaMixerNew = EtaMixerNew + ((EtaMixNewTemp - EtaMixmeanNew)*(EtaMixNewTemp - EtaMixmeanNew));
                }
        }
        fclose(fp);
}
//printf("Rg2erchk = %lf \n", Rg2er);
Re2er = sqrt(Re2er)/(NTRAJ*NDATA1);
Rg2er = sqrt(Rg2er)/(NTRAJ*NDATA1);
Xer = sqrt(Xer)/(NTRAJ*NDATA1);
//calculating std errors for each trajectory
for (k=0; k<NTRAJ; k++){
     Re2trajER += (Re2traj[k] - Re2mean)*(Re2traj[k] - Re2mean);
     Rg2trajER += (Rg2traj[k] - Rg2mean)*(Rg2traj[k] - Rg2mean);
   }
Re2trajER = sqrt(Re2trajER)/NTRAJ;
Rg2trajER = sqrt(Rg2trajER)/NTRAJ;

Re2ber = sqrt((Re2bsq/(NBLOCKS*NTRAJ) - Re2mean*Re2mean)/(NBLOCKS*NTRAJ));
Rg2ber = sqrt((Rg2bsq/(NBLOCKS*NTRAJ) - Rg2mean*Rg2mean)/(NBLOCKS*NTRAJ));
//EtaShrer = sqrt(EtaShrer)/(NTRAJ*NDATA1);
//EtaExter = sqrt(EtaExter)/(NTRAJ*NDATA1);
EtaMixer = sqrt(EtaMixer)/(NTRAJ*NDATA1);
NormalStressDiffer = sqrt(NormalStressDiffer)/(NTRAJ*NDATA1);
g11er = sqrt(g11er)/(NTRAJ*NDATA1); g12er = sqrt(g12er)/(NTRAJ*NDATA1); g13er = sqrt(g13er)/(NTRAJ*NDATA1);
g21er = sqrt(g21er)/(NTRAJ*NDATA1); g22er = sqrt(g22er)/(NTRAJ*NDATA1); g23er = sqrt(g23er)/(NTRAJ*NDATA1);
g31er = sqrt(g31er)/(NTRAJ*NDATA1); g32er = sqrt(g32er)/(NTRAJ*NDATA1); g33er = sqrt(g33er)/(NTRAJ*NDATA1);
//EtaMixerNew = sqrt(EtaMixerNew)/(NTRAJ*NDATA1);
// Calculation of volume fraction and error on it
// ////////////////////////////////////////////////
det = (g11mean*((g22mean*g33mean) - (g32mean*g23mean))) - (g12mean*((g21mean*g33mean) - (g23mean*g31mean))) + (g13mean*((g21mean*g32mean) - (g22mean*g31mean)));
phi = 4.0*M_PI*Nc*sqrt(3.0*det)/vol;
xin = 3.0*det;
detEr = det*(sqrt(((g11er/g11mean)*(g11er/g11mean)) + ((g22er/g22mean)*(g22er/g22mean)) + ((g33er/g33mean)*(g33er/g33mean))));
xinEr = 3.0*detEr;
phiEr = A*xinEr/(2.0*sqrt(xin));
//////////////////////////////////////////////////////
//PRINTING OUT RESULTS
printf("Results for NDATA1 = %d and NTRAJ = %d are: \n", NDATA1, NTRAJ);
printf("Re2 = %lf +/- %lf \n", Re2mean, Re2er);
printf("Rg2 = %lf +/- %lf \n", Rg2mean, Rg2er);
printf("X = %lf +/- %lf \n", Xmean, Xer);
//printf("EtaShr = %lf +/- %lf \n", EtaShrmean, EtaShrer);
//printf("EtaExt = %lf +/- %lf \n", EtaExtmean, EtaExter);
printf("EtaMixNew = %lf +/- %lf \n", EtaMixmean, EtaMixer);
printf("FirstNormalStreeDiff = %lf +/- %lf \n",NormalStressDiffmean, NormalStressDiffer);
printf("phi = %lf +/- %lf \n", phi, phiEr);
printf("Re2ber = %lf\t Rg2ber = %lf\n", Re2ber, Rg2ber);
//Printing the output Re2 and Rg2 corresponding to each trajectory
/*for (k=0; k<NTRAJ; k++){
     printf("Re2[%d] = %lf\t Rg2[%d] = %lf\n",k, Re2traj[k], k, Rg2traj[k]);
   }
printf("Re2trajER = %lf\t Rg2trajER = %lf\n", Re2trajER, Rg2trajER);
printf("Block averages\n");
for (k=0; k<NBLOCKS; k++){
    printf("%lf\n", Re2block[k]);
   }
for (k=0; k<NBLOCKS; k++){
    printf("%lf\n", Rg2block[k]);
   }
*/
//printf("g11 = %lf +/- %lf \n", g11mean, g11er);
//printf("g12 = %lf +/- %lf \n", g12mean, g12er);
//printf("g13 = %lf +/- %lf \n", g13mean, g13er);
//printf("g21 = %lf +/- %lf \n", g21mean, g21er);
//printf("g22 = %lf +/- %lf \n", g22mean, g22er);
//printf("g23 = %lf +/- %lf \n", g23mean, g23er);
//printf("g31 = %lf +/- %lf \n", g31mean, g31er);
//printf("g32 = %lf +/- %lf \n", g32mean, g32er);
//printf("g33 = %lf +/- %lf \n", g33mean, g33er);
//printf("det = %lf \t phi = %lf \t phiEr = %lf \n", det, phi, phiEr);
printf("For copying to results file \n");
printf("%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf  \t %lf \t %lf \n", Re2mean, Re2er, Rg2mean, Rg2er, Xmean, Xer, EtaMixmean, EtaMixer, NormalStressDiffmean, NormalStressDiffer, phi, phiEr);
}
