// Adapted from Jacques Kotze's code
// Ref: https://arxiv.org/abs/0803.0217v1/
// 2D ising model on hexagonal lattice
// Compile: g++ -O3 -std=c++11 2d-ising-2.0.cpp ran1.c -o ising2d-2.exe
// The program automatically reads input.2d-ising and perform calculations.
// Please see input.2d-ising file for explanation

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include "ran1.h"

double const kb=8.6173303e-2;                                 // kb = 8.6173303e-2 meV/K
long int seed=436675;                                         // seed for random number 

struct lat_type
{
    int x;    // lattice position x
    int y;    // lattice position y
    int iat;  // atom index in unit-cell
};

// Ising model for 2-d hexagonal lattice
//   H = -J \sum_{ij} S_i * S_j
//
class ising_2d_hexagonal
{
    private:
      double J1;                           // nearest neighbor exchange coupling, in meV
      double J2;                           // next nearest neighbor 
      double J3;                           // next next nearest neighbor
      double maxT;                         // starting point for temperature, Maximum Temperature
      double minT;                         // minimum temperature
      double Tstep;                        // size of steps for temperature loop
      long unsigned int nmcs;              // number of Monte Carlo steps
      int ntrans;                          // number of transient steps
      int size;                            // lattice size
      int nplot;                           // how many spin-config plots taken from mc simulation for each temperature, written in file spinconf*
      int mcplot;                          // if not equal to 1, then output <M(nstep)> in file "mc_T*" for each temperature
      static int const nat=2;              // number of atoms per unit cell
      int nsp;                             // number of spin points on lattice
      double norm;                         // normalization for averaging
      double norm2;                                 
      double norm4;                               
      int ***lat;
      bool param_set = false;
      bool is_param_set() { return param_set; }

      //Function for disregarding transient results
      void transient_results(double T)
      {
        lat_type pos;
        double de=0;
        for(int imc=1;imc<ntrans;imc++)    //Monte Carlo steps
        {
          for(int isp=1;isp<=nsp;isp++)    //Metropolis steps
          {
            choose_random_pos_lat(pos);
            if(test_flip(pos,de,T))
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
          double etot=0;
          for(int iat=1;iat<=nat;iat++)
          {
              pos.iat=iat;
              for(int y=size;y>=1;y--)
              {
                  pos.y=y;
                  for(int x=1;x<=size;x++)
                  {
                      pos.x=x;
                      etot+=energy_pos(pos);
                  }
              }
          }
          return etot;
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
          double energy = 0.0;
          if(pos.y==size)  up=1;
          else  up=pos.y+1;
          
          if(pos.y==1)  down=size;
          else  down=pos.y-1;
          
          if(pos.x==1)  left=size;
          else  left=pos.x-1;
          
          if(pos.x==size)  right=1;
          else  right=pos.x+1;
          
          // energy for a specfic position
          // contirbution from nearest neighbor
          if (pos.iat==1)
          {
            int neib=2;
            energy = -J1*lat[pos.x][pos.y][pos.iat] * ( lat[left][pos.y][neib] + \
              lat[pos.x][pos.y][neib] + \
              lat[pos.x][down][neib] );
          }
          else
          {
            int neib=1;
            energy = -J1*lat[pos.x][pos.y][pos.iat] * ( lat[right][pos.y][neib] + \
              lat[pos.x][pos.y][neib] + \
              lat[pos.x][up][neib] );
          }
          
          // contribution from next-nearest neighbor
          if (abs(J2) > 1e-3) 
          {
            energy += -J2*lat[pos.x][pos.y][pos.iat] * ( lat[right][pos.y][pos.iat] + \
                lat[left][pos.y][pos.iat] + \
                lat[pos.x][up][pos.iat]   + \
                lat[pos.y][down][pos.iat] + \
                lat[left][up][pos.iat]    + \
                lat[right][down][pos.iat] );
          }
          
          // contribution from next-next-nearest neighbor
          if (abs(J3) > 1e-3)
          {
            if (pos.iat==1)
            {
                int neib = 2;
                energy += -J3*lat[pos.x][pos.y][pos.iat] * ( lat[left][up][neib] + \
                  lat[right][down][neib]  + \
                  lat[left][down][neib]);
            }
            else
            {
                int neib = 1;
                energy += -J3*lat[pos.x][pos.y][pos.iat] * ( lat[left][up][neib] + \
                  lat[right][down][neib]  + \
                  lat[right][up][neib]);
            }
          }
          
          // std::cout << e << std::endl;
          return energy;
      }
      
      //function for testing the validity of flipping a spin at a selected position
      bool test_flip(lat_type pos, double &de, double T)
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

      void print_spin(int nstep, double temp)
      {
         FILE *plotfile;
         char filename[60]="spinconf";
         char tmpstr[20];
         sprintf(tmpstr, "_n%d", nstep);
         strcat(filename, tmpstr);
         sprintf(tmpstr, "_T%lf", temp);
         strcat(filename, tmpstr);
         plotfile = fopen(filename,"w");
         for(int x=1; x<=size; x++)
         {
             for(int y=1; y<=size; y++)
             {
                 for(int iat=1; iat<=nat; iat++)
                 {
                     fprintf( plotfile, " %12.5f  %12.5f  %5d \n ", x+y/2.0+(iat-1)/2.0, y*sqrt(3.0)/2.0+(iat-1)*sqrt(3.0)/6.0, lat[x][y][iat] );
                 }
             }
         }
         fclose(plotfile);
      }

    public:
      void set_parameters()
      {
          FILE *pfin;
          char line[256];
          pfin = fopen("input.2d-ising","r");
          if (pfin == NULL)
          {
              printf("2d-ising-hexagonal::set_parameters: Can't open input.2d-ising");
              exit(1);
          }
          else
          {
              fgets(line, sizeof(line), pfin); sscanf(line, "%lf", &J1);
              fgets(line, sizeof(line), pfin); sscanf(line, "%lf", &J2);
              fgets(line, sizeof(line), pfin); sscanf(line, "%lf", &J3);
              fgets(line, sizeof(line), pfin); sscanf(line, "%lf", &maxT);
              fgets(line, sizeof(line), pfin); sscanf(line, "%lf", &minT);
              fgets(line, sizeof(line), pfin); sscanf(line, "%lf", &Tstep);
              fgets(line, sizeof(line), pfin); sscanf(line, "%lu", &nmcs);
              fgets(line, sizeof(line), pfin); sscanf(line, "%d", &ntrans);
              fgets(line, sizeof(line), pfin); sscanf(line, "%d", &size);
              fgets(line, sizeof(line), pfin); sscanf(line, "%d", &nplot);
              fgets(line, sizeof(line), pfin); sscanf(line, "%d", &mcplot);
              nsp = size*size*nat;
              norm = 1.0/double(nsp);
              norm2 = norm*norm;
              norm4 = norm2*norm2;
              param_set = true;
          }
          fclose(pfin);
      }

      void allocate_lat()
      {
          if (is_param_set())
          {
              lat = new int**[size+1];
              for(int i=0; i<size+1; i++)
              {
                  lat[i] = new int*[size+1];
                  for (int j=0; j<size+1; j++)
                      lat[i][j] = new int[nat+1];
              }
          }
          else
          {
              printf("2d-ising-hexagonal::set_parameters: call set_parameters first.");
              exit(1);
          }
      }

      void deallocate_lat()
      {
          for(int i=0; i<size+1; i++)
          {
              for(int j=0;j<size+1;j++) delete [] lat[i][j];
              delete [] lat[i];
          }
          delete [] lat;
      }

      void initialize()
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

      void mc_simulate()
      {
          //declaring variables to be used in calculating the observables
          double E=0, E_avg=0, E2_avg=0, etot=0, e2tot=0;
          double M=0, M_avg=0, M2_avg=0, mtot=0, m2tot=0;
          double Mabs, Mabs_avg=0, M4_avg=0, mabstot=0, m4tot=0;
          FILE * pfile;
          FILE * mcfile;
          double de=0, T;
          int plotstep;
          lat_type pos;
          if (nplot > 0) plotstep = nmcs/nplot;
          pfile = fopen("2D-ising.dat","w");
          fprintf (pfile, "# J1 =%8.4lf\n# J2 =%8.4lf\n# J3 =%8.4lf\n", J1, J2, J3);
          fprintf (pfile, "# MaxT =%8.4lf\n# MinT =%8.4lf\n# Tstep=%8.4lf\n", maxT, minT, Tstep);
          fprintf (pfile, "# nmcs =%d\n# ntrans =%d\n", nmcs, ntrans);
          fprintf (pfile, "# System size L = %d x %d \n", size, size);
          fprintf (pfile, "# nplot = %d\n", nplot);
          fprintf (pfile, "# %7s ", " T" );                                      // temperature
          fprintf (pfile, "  %8s  %8s  %8s", "M_avg", "Mabs_avg", "M2_avg" );    // <M>; <|M|>; <M^2> per spin 
          fprintf (pfile, "  %8s", "dM/dT" );                                    // susceptibility per spin (X) = dM/dT
          fprintf (pfile, "  %8s", "d|M|/dT" );                                  // susceptibility per spin (X') = d|M|/dT
          fprintf (pfile, "  %8s  %8.4f ", "E_avg", "E^2_avg" );                 // <E>; <E^2> per spin
          fprintf (pfile, "  %8s",  "dE/dT" );                                   // heat capacity (C) per spin
          fprintf (pfile, "  %28s\n","U_L=1-((M4_avg)/(3*M2_avg))" );            // cumulant (U_L)

          //Tempearature loop
          for(T=maxT;T>=minT;T=T-Tstep)
          {
            char filename[40]="mc_T";
            char tmpstr[20];
            sprintf(tmpstr, "%3lf", T);
            strcat(filename, tmpstr);
            if (mcplot==1) mcfile = fopen(filename,"w");
            std::cout << "Working on temperature: " << T << " K." << std::endl;
          
            //transient function
            transient_results(T);
          
            //observables adopt equilibrated lattice configurations values
            M=total_magnetization();
            Mabs=abs(total_magnetization());
            E=total_energy();
          
            //initialize summation variables at each temperature step
            etot=0;
            e2tot=0;
            mtot=0;
            m2tot=0;
            m4tot=0;
            mabstot=0;
            
            //Monte Carlo loop
            for(int imc=1;imc<nmcs;imc++)
            {
              //Metropolis loop
              for(int isp=1;isp<=nsp;isp++)
              {
                choose_random_pos_lat(pos);
                if(test_flip(pos,de,T))
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
            
              if(mcplot==1) fprintf (mcfile, " %12d %12.4f \n", imc, M*norm);

              if(imc%plotstep == 0 && nplot != 0)
              {
                  print_spin(imc, T);
              }
            }
            if(mcplot==1) fclose(mcfile);
            
            //average observables
            E_avg=etot/nmcs;
            E2_avg=e2tot/nmcs;
            M_avg=mtot/nmcs;
            M2_avg=m2tot/nmcs;
            Mabs_avg=mabstot/nmcs;
            M4_avg=m4tot/nmcs;
            
            //output data to file
            fprintf (pfile, "  %7.4f ", T );                                      // temperature
            fprintf (pfile, "  %8.4f  %8.4f  %8.4f", M_avg, Mabs_avg, M2_avg );   // <M>; <|M|>; <M^2> per spin 
            fprintf (pfile, "  %8.4e", (M2_avg-(M_avg*M_avg))/T );                // susceptibility per spin (X) = dM/dT
            fprintf (pfile, "  %8.4e", (M2_avg-(Mabs_avg*Mabs_avg))/T );          // susceptibility per spin (X') = d|M|/dT
            fprintf (pfile, "  %8.4f  %8.4f ", E_avg, E2_avg );                   // <E>; <E^2> per spin
            fprintf (pfile, "  %8.4f",  (E2_avg-(E_avg*E_avg))/(T*T) );           // heat capacity (C) per spin
            fprintf (pfile, "  %8.4f\n",1-((M4_avg)/(3*M2_avg)) );                // cumulant (U_L)
          } // end of Temperature T loop

          fclose(pfile);
      }
};

int main(int argc, char **argv)
{
   ising_2d_hexagonal model;

   //initiliaze lattice to random configuration
   model.set_parameters();
   std::cout << " Read parameters from `input.2d-ising`. Done\n";
   model.allocate_lat();
   std::cout << " Allocate spin lattice. Done\n";
   model.initialize();
   std::cout << " Initialize spin lattice with random spins. Done\n\n";

   model.mc_simulate();
   std::cout << " Finished mc simulations\n";
   model.deallocate_lat();
   std::cout << " Deallcate spin lattice in memory. Done\n Finished!\n";
   return 0;
}

