# 2d-ising-hexagonal

* 2d-ising-1.0.cpp is written with functions. It is the first version. To use this version, you need to set up the parameters in the code, recompile it and run the calculation.

* 2d-ising-2.0.cpp is written with class. It reads input.2d-ising and performs calculations. 


Notes:

How to choose ntrans and nmcs?

At low temperatures, it takes a long time to reach equilibrium state.
The system may trap in local minimum states for a long time.
For example, if I set: J1=5.0, J2=J3=0, kb=1.0, size=25, and T=1.0, it takes about 70000 steps to reach global minimum state
Maybe implement cluster-flip dynamics.

The following are some parameter sets that work well, i.e., the system finds the global minimum with ferromagnetic ground state:

```C
// *** Input parameters ************************************* //
double const J1=5.56;                                         // nearest neighbor exchange coupling, in meV
double const J2=-0.88;                                        // next nearest neighbor 
double const J3= 0.02;                                        // next next nearest neighbor
double const kb=8.6173303e-2;                                 // kb = 8.6173303e-2 meV/K
double T =21.0;                                               // starting point for temperature
double const minT =21.0;                                      // minimum temperature
double const change = 0.25;                                   // size of steps for temperature loop
long unsigned int const nmcs= 80000;                          // number of Monte Carlo steps
int const ntrans= 8000;                                       // number of transient steps
const int size = 25;
// ********************************************************** //


// *** Input parameters ************************************* //
double const J1=5.18;                                         // nearest neighbor exchange coupling, in meV
double const J2=-0.84;                                        // next nearest neighbor 
double const J3=-0.11;                                        // next next nearest neighbor
double const kb=8.6173303e-2;                                 // kb = 8.6173303e-2 meV/K
double T =15.0;                                               // starting point for temperature
double const minT =15.0;                                      // minimum temperature
double const change = 0.25;                                   // size of steps for temperature loop
long unsigned int const nmcs= 800000;                         // number of Monte Carlo steps
int const ntrans= 8000;                                       // number of transient steps
const int size = 25;
// ********************************************************** //


// *** Input parameters ************************************* //
double const J1=4.90;                                         // nearest neighbor exchange coupling, in meV
double const J2=-0.82;                                        // next nearest neighbor 
double const J3=-0.19;                                        // next next nearest neighbor
double const kb=8.6173303e-2;                                 // kb = 8.6173303e-2 meV/K
double T =14.0;                                               // starting point for temperature
double const minT =14.0;                                      // minimum temperature
double const change = 0.25;                                   // size of steps for temperature loop
long unsigned int const nmcs= 800000;                         // number of Monte Carlo steps
int const ntrans= 20000;                                      // number of transient steps
const int size = 25;
// ********************************************************** //
 

// *** Input parameters ************************************* //
double const J1= 4.56;                                        // nearest neighbor exchange coupling, in meV
double const J2=-0.79;                                        // next nearest neighbor 
double const J3=-0.18;                                        // next next nearest neighbor
double const kb=8.6173303e-2;                                 // kb = 8.6173303e-2 meV/K
double T =13.0;                                               // starting point for temperature
double const minT =13.0;                                      // minimum temperature
double const change = 0.25;                                   // size of steps for temperature loop
long unsigned int const nmcs= 800000;                         // number of Monte Carlo steps
int const ntrans= 20000;                                      // number of transient steps
// ********************************************************** //
 

// *** Input parameters ************************************* //
double const J1= 4.34;                                         // nearest neighbor exchange coupling, in meV
double const J2=-0.76;                                        // next nearest neighbor 
double const J3=-0.30;                                        // next next nearest neighbor
double const kb=8.6173303e-2;                                 // kb = 8.6173303e-2 meV/K
double T =14.0;                                               // starting point for temperature
double const minT =14.0;                                      // minimum temperature
double const change = 0.25;                                   // size of steps for temperature loop
long unsigned int const nmcs= 800000;                         // number of Monte Carlo steps
int const ntrans= 20000;                                      // number of transient steps
// ********************************************************** //


// *** Input parameters ************************************* //
double const J1= 4.04;                                         // nearest neighbor exchange coupling, in meV
double const J2=-0.74;                                        // next nearest neighbor 
double const J3=-0.24;                                        // next next nearest neighbor
double const kb=8.6173303e-2;                                 // kb = 8.6173303e-2 meV/K
double T =15.0;                                               // starting point for temperature
double const minT =15.0;                                      // minimum temperature
double const change = 0.25;                                   // size of steps for temperature loop
long unsigned int const nmcs= 800000;                         // number of Monte Carlo steps
int const ntrans= 20000;                                      // number of transient steps
// ********************************************************** //
 

// *** Input parameters ************************************* //
double const J1= 3.92;                                         // nearest neighbor exchange coupling, in meV
double const J2=-0.71;                                        // next nearest neighbor 
double const J3=-0.22;                                        // next next nearest neighbor
double const kb=8.6173303e-2;                                 // kb = 8.6173303e-2 meV/K
double T =15.0;                                               // starting point for temperature
double const minT =15.0;                                      // minimum temperature
double const change = 0.25;                                   // size of steps for temperature loop
long unsigned int const nmcs= 20000;                         // number of Monte Carlo steps
int const ntrans= 1000;                                      // number of transient steps
// ********************************************************** //


// *** Input parameters ************************************* //
double const J1= 3.89;                                         // nearest neighbor exchange coupling, in meV
double const J2=-0.55;                                        // next nearest neighbor 
double const J3= 0.02;                                        // next next nearest neighbor
double const kb=8.6173303e-2;                                 // kb = 8.6173303e-2 meV/K
double T =15.0;                                               // starting point for temperature
double const minT =15.0;                                      // minimum temperature
double const change = 0.25;                                   // size of steps for temperature loop
long unsigned int const nmcs= 20000;                         // number of Monte Carlo steps
int const ntrans= 1000;                                      // number of transient steps
// ********************************************************** //


// *** Input parameters ************************************* //
double const J1= 3.99;                                         // nearest neighbor exchange coupling, in meV
double const J2=-0.44;                                        // next nearest neighbor 
double const J3= 0.22;                                        // next next nearest neighbor
double const kb=8.6173303e-2;                                 // kb = 8.6173303e-2 meV/K
double T =15.0;                                               // starting point for temperature
double const minT =15.0;                                      // minimum temperature
double const change = 0.25;                                   // size of steps for temperature loop
long unsigned int const nmcs= 20000;                         // number of Monte Carlo steps
int const ntrans= 1000;                                      // number of transient steps
// ********************************************************** //


``` 


