#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "ran1.h"

struct lat_type
{
    int x;    // lattice position x
    int y;    // lattice position y
    int iat;  // atom index in unit-cell
};

const int size = 25;                  // lattice size
const int nat=2;                      // number of atoms per unit cell
const int nsp=nat*size*size;          // number of spin points on lattice
long unsigned int const nmcs=40000;         // number of Monte Carlo steps
int const ntrans=80000;                     // number of transient steps
double const norm=(1.0/double(nsp));  // normalization for averaging
double const norm2=norm*norm;
double const norm4=norm2*norm2;

double T = 8.0;                       // starting point for temperature
const double minT = 5.5;              // minimum temperature
double change = 0.5;                  // size of steps for temperature loop
int lat[size+1][size+1][nat+1];       // 2d lattice for spins
long int seed=436675;                 // seed for random number 
// Ising model:
//   H = -J \sum_{ij} S_i * S_j
double J=5.0;   // 10.0;                      // exchange coupling, in meV
double kb=1.0;  // 8.6173303e-2;              // meV/K

void initialize(int lat[size+1][size+1][nat+1])
{
    for(int iat=1;iat<=nat;iat++)
    {
        for(int y=size;y>=1;y--)
        {
            for(int x=1;x<=size;x++)
            {
                if(ran1(&seed)>=0.5)
                    lat[x][y][iat]= 1;
                else
                    lat[x][y][iat]=-1;
            }
        }
    }
}

void choose_random_pos_lat(lat_type &pos)
{
    pos.x=(int)ceil(ran1(&seed)*(size));
    pos.y=(int)ceil(ran1(&seed)*(size));
    pos.iat=(int)ceil(ran1(&seed)*(nat));
    // std::cout<<pos.x<<"  "<<pos.y<<"  "<<pos.iat<<std::endl;
    if(pos.x>size || pos.y>size || pos.iat>nat)
    {
        std::cout<<"error in array size allocation for random position on lattice";
        std::cout<<pos.x<<", "<<pos.y<<", "<<pos.iat<<std::endl;
        exit;
    }
}

double energy_pos(lat_type &pos)
{
    //periodic boundary conditions
    int up, down, left, right;
    double energy;
    if(pos.y==size)  up=1;
    else  up=pos.y+1;
    
    if(pos.y==1)  down=size;
    else  down=pos.y-1;
    
    if(pos.x==1)  left=size;
    else  left=pos.x-1;
    
    if(pos.x==size)  right=1;
    else  right=pos.x+1;
    
    // energy for specfic position
    if (pos.iat==1)
    {
      int neib=2;
      energy = -J*lat[pos.x][pos.y][pos.iat] * ( lat[left][pos.y][neib] + \
        lat[pos.x][pos.y][neib] + lat[pos.x][down][neib] );
    }
    else
    {
      int neib=1;
      energy = -J*lat[pos.x][pos.y][pos.iat] * ( lat[right][pos.y][neib] + \
        lat[pos.x][pos.y][neib]+lat[pos.x][up][neib] );
    }
    // std::cout << e << std::endl;
    return energy;
}

//function for testing the validity of flipping a spin at a selected position
bool test_flip(lat_type pos, double &de)
{
  // If the spin in site i flips, the change in total energy is: deltaE = -2*J*\sum_{j \in neib} S_i*S_j
  de=-2.0*energy_pos(pos);        // change in energy for specific spin
  if(de<0.0) 
    return true;                  // flip due to lower energy
  else if (ran1(&seed)<exp(-de/kb/T)) 
    return true;                  // flip due to heat bath
  else
    return false;                 // no flip
}

//flip spin at given position
void flip(lat_type pos)
{
  lat[pos.x][pos.y][pos.iat]=-lat[pos.x][pos.y][pos.iat];
}

//function for disregarding transient results
// Discard the transient results
void transient_results()
{
  lat_type pos;
  double de=0;
  for(int imc=1;imc<ntrans;imc++)      //Monte Carlo steps
  {
    for(int isp=1;isp<=nsp;isp++)      //Metropolis steps
    {
      choose_random_pos_lat(pos);
      if(test_flip(pos,de))
      {
        flip(pos);
      }
    }
  }
}

//function for calculating total magnetization of lattice 
double total_magnetization()
{
  double m=0;
  for(int iat=1;iat<=nat;iat++)
  {
     for(int y=size;y>=1;y--)
     {
       for(int x=1;x<=size;x++)
       {
         m+=lat[x][y][iat];
       }
     }
  }
  return m;
}

// Note: This function calculate 2 times the real total energy
//     because the coupling between each spin-pairs are counted twice
double total_energy()
{
  lat_type pos;
  double e=0;
  for(int iat=1;iat<=nat;iat++)
  {
    pos.iat=iat;
    for(int y=size;y>=1;y--)
    {
      pos.y=y;
      for(int x=1;x<=size;x++)
      {
        pos.x=x;
        e+=energy_pos(pos);
      }
    }
  }
  return e;
}

int main(int argc, char **argv)
{
   //declaring variables to be used in calculating the observables
   double E=0, E_avg=0, E2_avg=0, etot=0, e2tot=0;
   double M=0, M_avg=0, M2_avg=0, mtot=0, m2tot=0;
   double Mabs=0, Mabs_avg=0, M4_avg=0, mabstot=0, m4tot=0;
   double de=0;

   std::ofstream DATA("2D-ising.dat", std::ofstream::out);
   lat_type pos;

   //initiliaze lattice to random configuration
   initialize(lat);

   //Tempearature loop
   for(;T>=minT;T=T-change)
   {
     std::cout << "Working on temperature: " << T << std::endl;

     //transient function
     transient_results();

     //observables adopt equilibrated lattice configurations values
     M=total_magnetization();
     Mabs=abs(total_magnetization());
     E=total_energy();

     //initialize summation variables at each temperature step
     etot=0;
     e2tot=0;
     mtot=0;
     m4tot=0;
     mabstot=0;

     //Monte Carlo loop
     for(int imc=1;imc<nmcs;imc++)
     {
       //Metropolis loop
       for(int isp=1;isp<=nsp;isp++)
       {
         choose_random_pos_lat(pos);
         if(test_flip(pos,de))
         {
           flip(pos);
           //adjust observables
           E+=2.0*de; // This is actually 2*real total energy
           M+=2.0*lat[pos.x][pos.y][pos.iat]; // If the spin at site i flips, the total magnetization changes by 2*Si
         }
       }
       //keep summation of observables
       etot+=E/2.0*norm;
       e2tot+=E/2.0*E/2.0*norm2;
       mtot+=M*norm;
       m2tot+=M*M*norm2;
       m4tot+=M*M*M*M*norm4;
       mabstot+=(sqrt(M*M))*norm;

       // DATA<<"  "<<imc<<"  "<<M<<std::endl; 
     }

     //average observables
     E_avg=etot/nmcs;
     E2_avg=e2tot/nmcs;
     M_avg=mtot/nmcs;
     M2_avg=m2tot/nmcs;
     Mabs_avg=mabstot/nmcs;
     M4_avg=m4tot/nmcs;

     //output data to file
     DATA<<"#"<<T<<                                      //temperature
     "\t"<<M_avg<<"\t"<<Mabs_avg<<"\t"<<M2_avg<<        //<M>; <|M|>; <M^2> per spin 
     "\t"<<(M2_avg-(M_avg*M_avg*nsp))/T<<               //susceptibility per spin (X) = dM/dT
     "\t"<<(M2_avg-(Mabs_avg*Mabs_avg*nsp))/T<<         //susceptibility per spin (X') = d|M|/dT
     "\t"<<E_avg<<"\t"<<E2_avg<<                        //<E>; <E^2> per spin
     "\t"<<(E2_avg-(E_avg*E_avg*nsp))/(T*T)<<           //heat capacity (C) per spin
     "\t"<<1-((M4_avg)/(3*M2_avg))<<std::endl;          //cumulant (U_L)
   }
   
   DATA.close();

   return 0;
}
