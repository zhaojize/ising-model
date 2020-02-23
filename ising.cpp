#include <iostream>
#include <math.h>
#include <random>
using namespace std;

#define  LX   32

#define  LY   32

#define NSITE    (LX*LY)

#define TEMP   (1.0/0.441) 

#define  J   (-1.0)

#define  MCS  100000000

#define  NB   4

int spinconf[LX*LY];

int nbr[LX*LY][NB];

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(0.0, 1.0);

inline double mydrand48()
{
    return dis(gen);
}

void initializeConf(int nbrr[LX*LY][NB], int* sconf, const int &size, double &confenergy, int &confmag)
{
    confmag = 0;

    for (int i = 0; i < size; i++)
    {
	if (mydrand48() >= 0.5) sconf[i] = 1;
        else sconf[i] = -1;	

	confmag += sconf[i];
    }

    confenergy = 0.0;	

    for (int x = 0; x < LX; x++)
    {
	int xm = (x-1+LX)%LX;
        int xp = (x+1)%LX;

	for (int y = 0; y < LY; y++)
	{
	    int ym = (y-1+LY)%LY;
            int yp = (y+1)%LY;	   

	    int xy = x*LY+y; 
	    nbrr[xy][0] = xm*LY+y;
	    nbrr[xy][1] = xp*LY+y;
	    nbrr[xy][2] = x*LY+ym;
	    nbrr[xy][3] = x*LY+yp;

	    confenergy += J*sconf[xy]*(sconf[nbrr[xy][0]]+sconf[nbrr[xy][2]]);    
	}	
    }
}

void updateConf(const int &xy, int* sconf, double &confenergy, int &confmag)
{
    double eo = J*sconf[xy]*(sconf[nbr[xy][0]]+sconf[nbr[xy][1]]+sconf[nbr[xy][2]]+sconf[nbr[xy][3]]); 
    double en = -eo;  

    if (en <= eo || exp(-(en-eo)/TEMP) > mydrand48()) 
    {
	sconf[xy] = -sconf[xy];    
        confenergy += en-eo;

	confmag += sconf[xy]+sconf[xy];
    }
}

int main(void)
{
    double confenergy = 0.0;

    int confmag = 0;

    initializeConf(nbr, spinconf, LX*LY, confenergy, confmag);

    for (int wm =0; wm < 20000; wm++)
    {
	for (int xy = 0; xy < LX*LY; xy++) updateConf(xy, spinconf, confenergy, confmag);	
    }

    double energysum = 0.0;

    double energy2sum = 0.0;

    double magsum = 0.0;

    double mag2sum = 0.0;

    for (long ll = 1; ll <= MCS; ll++)
    {
	for (int i = 0; i < LX*LY; i++)
	{
	    updateConf(i, spinconf, confenergy, confmag);	
        }

        energysum += confenergy;

	energy2sum += confenergy * confenergy;

	magsum += confmag/(NSITE*1.0);

	mag2sum += (confmag/(1.0*NSITE)) * (confmag/(1.0*NSITE));

	if (ll%1000000 == 0)
	{
            cout << "step/10^6 = " << ll/1000000 << endl;	
            cout << "      energy is "  << energysum/ll/2.0/NSITE << endl;	    
	    cout << "      energy/T is : " << energysum/ll/TEMP/2.0/NSITE << endl;
	    cout << "      Cv is  : "  << (energy2sum/ll-(energysum/ll)*(energysum/ll))/(TEMP*TEMP)/NSITE << endl; 
            cout << "      Average magnetization is " << magsum/ll << endl;
	    cout << "      Magnetization squared is " << mag2sum/ll << endl;
	}
    }

    cout<<"energy/T is : "<<energysum/MCS/(LX*LY)/TEMP/2.0<<endl;

    cout<<"Average magnetization is " <<magsum/MCS/(LX*LY)<<endl;

    return 1;
}
