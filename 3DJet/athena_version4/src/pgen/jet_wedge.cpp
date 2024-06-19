//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file jet.cpp
//! \brief Sets up a nonrelativistic jet introduced through L-x1 boundary (left edge)
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include <random>
#include <iostream> 
#include <sstream> 
#include <fstream>
#include <algorithm>



// BCs on L-x1 (left edge) of grid with jet inflow conditions
void JetInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh);


// BCs on L-x1 (left edge) of grid with jet inflow conditions
void JetInnerX2_2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
// Make radius of jet and jet variables global so they can be accessed by BC functions
Real d_amb, p_amb, vx_amb, vy_amb, vz_amb, bx_amb, by_amb, bz_amb, den, pres;
Real r_jet, d_jet, p_jet, vx_jet, vy_jet, vz_jet, bx_jet, by_jet, bz_jet, wall_thick;
Real x1_0, x2_0, x3_0, x_0, y_0, z_0, gm1;
Real nuair, nu_h2;
Real Mw_fuel, Mw_air;
Real N; //number of steps before averaging **Turbulence**
int counter;

} // namespace

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(29);
  SetUserOutputVariableName(0, "pressure");
  SetUserOutputVariableName(1, "density");
  SetUserOutputVariableName(2, "vx");
  SetUserOutputVariableName(3, "vy");
  SetUserOutputVariableName(4, "vz");
  SetUserOutputVariableName(5, "temperature");
  SetUserOutputVariableName(6, "mach");
  SetUserOutputVariableName(7, "scalars");
  SetUserOutputVariableName(8, "schlieren");
  SetUserOutputVariableName(9, "Mean_velx");
  SetUserOutputVariableName(10, "Mean_vely");
  SetUserOutputVariableName(11, "Mean_velz");
  SetUserOutputVariableName(12, "Mean_temp");
  SetUserOutputVariableName(13, "Mean_massfraction");
  SetUserOutputVariableName(14, "Mean_dens");
  SetUserOutputVariableName(15, "Mean_mach");
  SetUserOutputVariableName(16, "Mean_Ur");
  SetUserOutputVariableName(17, "Fluct_vx");
  SetUserOutputVariableName(18, "Fluct_vy");
  SetUserOutputVariableName(19, "Fluct_vz");
  SetUserOutputVariableName(20, "Fluct_Ur");
  SetUserOutputVariableName(21, "Re_xx");
  SetUserOutputVariableName(22, "Re_yy");
  SetUserOutputVariableName(23, "Re_zz");
  SetUserOutputVariableName(24, "Re_xy");
  SetUserOutputVariableName(25, "Re_xz");
  SetUserOutputVariableName(26, "Re_zy");
  SetUserOutputVariableName(27, "Re_rr");
  SetUserOutputVariableName(28, "Re_yr");

  return;
  return;
}

int RefinementCondition(MeshBlock *pmb);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // initialize global variables
  d_amb  = pin->GetReal("problem", "d");
  p_amb  = pin->GetReal("problem", "p");
  vx_amb = pin->GetReal("problem", "vx");
  vy_amb = pin->GetReal("problem", "vy");
  vz_amb = pin->GetReal("problem", "vz");
  
  d_jet  = pin->GetReal("problem", "djet");
  p_jet  = pin->GetReal("problem", "pjet");
  vx_jet = pin->GetReal("problem", "vxjet");
  vy_jet = pin->GetReal("problem", "vyjet");
  vz_jet = pin->GetReal("problem", "vzjet");

  r_jet = pin->GetReal("problem", "rjet");
  Mw_fuel = pin->GetReal("problem", "Mwfuel");
  Mw_air = pin->GetReal("problem", "Mwair");
  N = pin->GetReal("problem", "N_steps");

 // wall_thick = pin->GetReal("problem", "wall_thick");

 // nuair = pin->GetReal("problem", "nuair");
 // nu_h2 = pin->GetReal("problem", "nu_h2");

  // get coordinates of center of jet --cylindrical/wedge r, theta and z
  // get coordinates of center of blast, and convert to Cartesian if necessary
  x1_0   = 0.0; // pin->GetOrAddReal("problem", "x1_0", 0.0);
  x2_0   = 0.0; // pin->GetOrAddReal("problem", "x2_0", 0.0);
  x3_0   = 0.0; // pin->GetOrAddReal("problem", "x3_0", 0.0);

  x_0 = x1_0*std::cos(x3_0);
  y_0 = x2_0;
  z_0 = x1_0*std::sin(x3_0);

//counter
int counter;


  // enroll boundary value function pointers
  EnrollUserBoundaryFunction(BoundaryFace::inner_x2, JetInnerX2);


  if(adaptive==true)
      EnrollUserRefinementCondition(RefinementCondition);

  return;
}

//Vectors
const int binNum = 100;
std::vector<std::vector<int>> local_histogram(binNum, std::vector<int>(binNum, 0));
// Create a global histogram to store the results
std::vector<std::vector<int>> global_histogram(binNum, std::vector<int>(binNum, 0));
//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Jet problem
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  gm1 = peos->GetGamma() - 1.0;
  // initialize conserved variables
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) 
      {
            phydro->u(IDN,k,j,i) = d_amb;
            phydro->u(IM1,k,j,i) = d_amb*vx_amb;
            phydro->u(IM2,k,j,i) = d_amb*vy_amb;
            phydro->u(IM3,k,j,i) = d_amb*vz_amb;
            phydro->u(IEN,k,j,i) = p_amb/gm1 + 0.5*d_amb*(SQR(vx_amb)+SQR(vy_amb)+SQR(vz_amb));
            if (NSCALARS > 0){
              pscalars->s(0, k, j, i) = 0.0;
              pscalars->r(0, k, j, i) = 0.0;
            }                                            
      }
    }
  }
  // JetInnerX2_2(this,pcoord,phydro->u,pfield->b,pmy_mesh->time,pmy_mesh->dt,is,ie,js,je,ks,ke,NGHOST);
  return;
}

void MeshBlock::UserWorkInLoop()
{
    Real T_new;
	  Real R_mix_n;
    // Real meanT;
    //  Real meanY;
    // Define the range and discretization parameters
/*const int numBins = 100;
const double minY = 0.0, maxY = 1.0;
const double minT = 0.0, maxT = 1000.0;*/
//std::cout << "Before calculation - vx_jet: " << vx_jet << ", vz_jet: " << vz_jet << std::endl;
	//  MeshBlock *pmb;
	  for (int k = ks; k <= ke; ++k) {
	    for (int j = js; j <= je; ++j) {
	      for (int i = is; i <= ie; ++i) {
	           if(NSCALARS>0){
	           Real u2 = phydro->w(IVX, k, j, i) * phydro->w(IVX, k, j, i) + phydro->w(IVY, k, j, i) * phydro->w(IVY, k, j, i) + phydro->w(IVZ, k, j, i) * phydro->w(IVZ, k, j, i);
	           Real c2 = phydro->w(IPR, k, j, i) / phydro->w(IDN, k, j, i) * peos->GetGamma();
	           Real mach = sqrt(u2 / c2); //calculate mach
	           Real Y = pscalars->r(0, k, j, i); //calculate mass
	           Real M_x = 1.0 / (((1.0-Y) / Mw_air) + (Y / Mw_fuel));
	           Real R_mix = 8.314/M_x;
	           T_new = phydro->w(IPR, k, j, i)/phydro->w(IDN, k, j, i)/R_mix; //calculate Temperature

             Real U_r = std::sqrt(std::pow(phydro->w(IVX, k, j, i), 2) + std::pow(phydro->w(IVZ, k, j, i), 2));

             //std::cout << "U_r=" <<  U_r << std::endl;

             if (std::abs(pmy_mesh->time >= 0.0000224) || 
             std::abs(pmy_mesh->time >= 0.0000544) || 
             std::abs(pmy_mesh->time >= 0.0000864) || 
             std::abs(pmy_mesh->time >= 0.0001184)) {


              if (counter<N) {
              // sum up velocities, Density, mach temp and mass fraction  
              phydro->U_sum(0, k, j, i) += phydro->w(IVX, k, j, i);
			        phydro->U_sum(1, k, j, i) += phydro->w(IVY, k, j, i);
			        phydro->U_sum(2, k, j, i) += phydro->w(IVZ, k, j, i);
			        phydro->U_sum(3, k, j, i) += T_new;
			        phydro->U_sum(4, k, j, i) += Y;
			        phydro->U_sum(5, k, j, i) += phydro->w(IDN, k, j, i);
			        phydro->U_sum(6, k, j, i) += mach;
			        phydro->U_sum(7, k, j, i) += U_r;


           // std::cout << "Accumulated Values: IVX=" <<  phydro->U_sum(0, k, j, i) << std::endl;
              }
              }
           
               
	 	    phydro->U_mean(0, k, j, i) = phydro->U_sum(0, k, j, i)/N;
        phydro->U_mean(1, k, j, i) = phydro->U_sum(1, k, j, i)/N;
        phydro->U_mean(2, k, j, i) = phydro->U_sum(2, k, j, i)/N;
        phydro->U_mean(3, k, j, i) = phydro->U_sum(3, k, j, i)/N;        
        phydro->U_mean(4, k, j, i) = phydro->U_sum(4, k, j, i)/N;
        phydro->U_mean(5, k, j, i) = phydro->U_sum(5, k, j, i)/N;
        phydro->U_mean(6, k, j, i) = phydro->U_sum(6, k, j, i)/N;
        phydro->U_mean(7, k, j, i) = phydro->U_sum(7, k, j, i)/N;

if (std::abs(pmy_mesh->time >= 0.0000224) || 
 std::abs(pmy_mesh->time >= 0.0000544) || 
 std::abs(pmy_mesh->time >= 0.0000864) || 
 std::abs(pmy_mesh->time >= 0.0001184)) {
  if(counter<N){
 //calculate fluctuations

phydro->U_fluc(0, k, j, i) = phydro->w(IVX, k, j, i)-phydro->U_mean(0, k, j, i);
phydro->U_fluc(1, k, j, i) = phydro->w(IVY, k, j, i)-phydro->U_mean(1, k, j, i);
phydro->U_fluc(2, k, j, i) = phydro->w(IVZ, k, j, i)-phydro->U_mean(2, k, j, i);
phydro->U_fluc(3, k, j, i) = U_r-phydro->U_mean(7, k, j, i);


    phydro->fluc_sum(0,k,j,i)+=phydro->U_fluc(0, k, j, i)*phydro->U_fluc(0, k, j, i); //fluctuation multiplication and sum up x dir
    phydro->fluc_sum(1,k,j,i)+=phydro->U_fluc(1, k, j, i)*phydro->U_fluc(1, k, j, i); //fluctuation multiplication and sum up y dir
    phydro->fluc_sum(2,k,j,i)+=phydro->U_fluc(2, k, j, i)*phydro->U_fluc(2, k, j, i); //fluctuation multiplication and sum up z dir
    phydro->fluc_sum(3,k,j,i)+=phydro->U_fluc(0, k, j, i)*phydro->U_fluc(1, k, j, i); //fluctuation multiplication and sum up xy dir
    phydro->fluc_sum(4,k,j,i)+=phydro->U_fluc(0, k, j, i)*phydro->U_fluc(2, k, j, i); //fluctuation multiplication and sum up xz dir
    phydro->fluc_sum(5,k,j,i)+=phydro->U_fluc(1, k, j, i)*phydro->U_fluc(2, k, j, i); //fluctuation multiplication and sum up yz dir
    phydro->fluc_sum(6,k,j,i)+=phydro->U_fluc(3, k, j, i)*phydro->U_fluc(3, k, j, i); //fluctuation multiplication and sum up u_r*_u_r dir
    phydro->fluc_sum(7,k,j,i)+=phydro->U_fluc(3, k, j, i)*phydro->U_fluc(1, k, j, i); //fluctuation multiplication and sum up u_r*_u_y dir

  }
  }
  {
  phydro->Rey_stress(0,k,j,i)=phydro->fluc_sum(0,k,j,i)/N; // Reynolds_stress xx dir
  phydro->Rey_stress(1,k,j,i)=phydro->fluc_sum(1,k,j,i)/N;  // Reynolds_stress yy dir
  phydro->Rey_stress(2,k,j,i)=phydro->fluc_sum(2,k,j,i)/N;  // Reynolds_stress zz dir
  phydro->Rey_stress(3,k,j,i)=phydro->fluc_sum(3,k,j,i)/N;  // Reynolds_stress xy dir
  phydro->Rey_stress(4,k,j,i)=phydro->fluc_sum(4,k,j,i)/N;  // Reynolds_stress xz dir
  phydro->Rey_stress(5,k,j,i)=phydro->fluc_sum(5,k,j,i)/N;  // Reynolds_stress yz dir
  phydro->Rey_stress(6,k,j,i)=phydro->fluc_sum(6,k,j,i)/N;  // Reynolds_stress rr dir
  phydro->Rey_stress(7,k,j,i)=phydro->fluc_sum(7,k,j,i)/N;  // Reynolds_stress yr dir
  
  }
               }

	  if (std::abs(pmy_mesh->time > 3.205e-05 && pmy_mesh->time <0.0000544) || 
       std::abs(pmy_mesh->time > 6.41e-05 && pmy_mesh->time < 0.0000864) || 
       std::abs(pmy_mesh->time > 9.61e-05 && pmy_mesh->time < 0.0001184)) {
  
    phydro->U_sum(0, k, j, i) = 0.0;
    phydro->U_sum(1, k, j, i) = 0.0;
    phydro->U_sum(2, k, j, i) = 0.0;
    phydro->U_sum(3, k, j, i) = 0.0;
    phydro->U_sum(4, k, j, i) = 0.0;
    phydro->U_sum(5, k, j, i) = 0.0;
    phydro->U_sum(6, k, j, i) = 0.0;

phydro->U_mean(0, k, j, i)=0.0;
phydro->U_mean(1, k, j, i)=0.0;
phydro->U_mean(2, k, j, i)=0.0;
phydro->U_mean(3, k, j, i)=0.0;
phydro->U_mean(4, k, j, i)=0.0;
phydro->U_mean(5, k, j, i)=0.0;
phydro->U_mean(6, k, j, i)=0.0;

 phydro->Rey_stress(0,k,j,i)=0;
 phydro->Rey_stress(1,k,j,i)=0;
 phydro->Rey_stress(2,k,j,i)=0;
 phydro->Rey_stress(3,k,j,i)=0;
 phydro->Rey_stress(4,k,j,i)=0;
 phydro->Rey_stress(5,k,j,i)=0;

phydro->fluc_sum(0,k,j,i)=0;
phydro->fluc_sum(1,k,j,i)=0;
phydro->fluc_sum(2,k,j,i)=0;
phydro->fluc_sum(3,k,j,i)=0;
phydro->fluc_sum(4,k,j,i)=0;
phydro->fluc_sum(5,k,j,i)=0;
phydro->fluc_sum(6,k,j,i)=0;
phydro->fluc_sum(7,k,j,i)=0;

counter = 0;
    }
	 }
	 }
  }
}

void Mesh::UserWorkInLoop() {
    Real Tmax = 0, T;
    Real R_mix;
    Real tmax_global;
    Real meanT;
    Real meanY;
    Real max_x = 0; // Initialize to lowest possible value
    Real max_y = 0;
    Real max_z = 0;


    // Define the range and discretization
    const int numBins = 100;
    const double minY = 0.0, maxY = 1.0;
    const double minT = 0.0, maxT = 1000.0;

    // Flags to indicate whether the histograms have been initialized for each time step
bool histogram1Initialized = false;
bool histogram2Initialized = false;
bool histogram3Initialized = false;
bool histogram4Initialized = false;


// Track the current time
double currentTime = time;

// Check if the time has crossed a threshold
const double tolerance = 1e-10;

// Determine the file name based on the current time
std::string histogramFileName;
if (std::abs(currentTime - 3.2e-05) < tolerance) {
    histogramFileName = "histogram.txt";
    if (!histogram1Initialized) {
        // Initialize histogram1 only once at the beginning of the simulation
        local_histogram.assign(numBins, std::vector<int>(numBins, 0));
        histogram1Initialized = true;
    }
} else if (std::abs(currentTime - 6.4e-05) < tolerance) {
    histogramFileName = "histogram2.txt";
    if (!histogram2Initialized) {
        // Initialize histogram2 only once at the beginning of the simulation
        local_histogram.assign(numBins, std::vector<int>(numBins, 0));
        histogram2Initialized = true;
    }
} else if (std::abs(currentTime - 9.6e-05) < tolerance) {
    histogramFileName = "histogram3.txt";
    if (!histogram3Initialized) {
        // Initialize histogram3 only once at the beginning of the simulation
        local_histogram.assign(numBins, std::vector<int>(numBins, 0));
        histogram3Initialized = true;
    }
} else if (std::abs(currentTime - 0.0001278) < tolerance) {
    histogramFileName = "histogram4.txt";
    if (!histogram4Initialized) {
        // Initialize histogram4 only once at the beginning of the simulation
        local_histogram.assign(numBins, std::vector<int>(numBins, 0));
        histogram4Initialized = true;
    }
}

// Declare variables to store the coordinates of T_max
//int max_i = -1, max_j = -1, max_k = -1;

    for (int bn = 0; bn < nblocal; ++bn) {
        MeshBlock *pmb = my_blocks(bn);
        for (int k = pmb->ks; k <= pmb->ke; k++) {
            for (int j = pmb->js; j <= pmb->je; j++) {
                for (int i = pmb->is; i <= pmb->ie; i++) {
                    Coordinates *pco;
                    Real x =pmb->pcoord->x1v(i);
                    Real y =pmb->pcoord->x2v(j);
                    Real z =pmb->pcoord->x3v(k);
                    //std::cout << "x " << pco->x1v(i)<< std::endl;


                    Real Y = pmb->pscalars->r(0, k, j, i);
                    Real M_x = 1.0 / (((1.0 - Y) / Mw_air) + (Y / Mw_fuel));
                    Real R_mix = 8.314 / M_x;
                    T = pmb->phydro->w(IPR, k, j, i) / pmb->phydro->w(IDN, k, j, i) / R_mix;
                   // Tmax = std::max(Tmax, T);
                   if (T > Tmax) {
                    Tmax = T;
                    max_x = x;
                    max_y = y;
                    max_z = z;
                    }
                    

                    const double tolerance = 1e-10;
                    if ((std::abs(pmb->pmy_mesh->time - 3.2e-05) < tolerance) ||
   (std::abs(pmb->pmy_mesh->time - 6.4e-05) < tolerance) ||
   (std::abs(pmb->pmy_mesh->time - 9.6e-05) < tolerance) ||
   (std::abs(pmb->pmy_mesh->time - 0.0001278) < tolerance)) {

                        meanT = pmb->phydro->U_mean(3, k, j, i); // mean Temp
                        meanY = pmb->phydro->U_mean(4, k, j, i); // mean Mass fraction

                      //  std::cout << " meanY: " << meanY << ", meanT: " << meanT << std::endl;

                      // Discretize meanY and meanT using round function
                        int binMeanY = std::round((meanY - minY) / (maxY - minY) * (numBins));
                        binMeanY = std::max(0, std::min(binMeanY, numBins - 1));
                        int binMeanT = std::round((meanT - minT) / (maxT - minT) * (numBins));
                        binMeanT = std::max(0, std::min(binMeanT, numBins - 1));

                       // std::cout << "After rounding - binMeanY: " << binMeanY << ", binMeanT: " << binMeanT << std::endl;

                        // Check if the calculated indices are within bounds
                        if (binMeanY < 0 || binMeanY >= numBins || binMeanT < 0 || binMeanT >= numBins) {
                         //   std::cerr << "Error: Indices out of bounds - binMeanY: " << binMeanY << ", binMeanT: " << binMeanT << std::endl;
                            // Handle the error as needed
                        } else {
                            // Ensure local_histogram has sufficient size before accessing indices
                           // if (local_histogram.size() > binMeanY && local_histogram[binMeanY].size() > binMeanT) {
                                // Increment the corresponding entry in the histogram
                                local_histogram[binMeanY][binMeanT]++;
                            }
                           //
                             
                          
                        }
                    }
                }
            }
        }

/*Real max_x_local = max_x; // Store local max_x value
Real max_y_local = max_y; // Store local max_y value
Real max_z_local = max_z; // Store local max_z value

Real global_max_x, global_max_y, global_max_z; // Variables to store global max coordinates*/


#ifdef MPI_PARALLEL
    // Flatten local_histogram into a one-dimensional array
     MPI_Allreduce(&Tmax, &tmax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     // Perform MPI reduction to find global maximum coordinates
   //  MPI_Allreduce(&max_x_local, &global_max_x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   //  MPI_Allreduce(&max_y_local, &global_max_y, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   //  MPI_Allreduce(&max_z_local, &global_max_z, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    std::vector<int> local_histogram_flat;
    for (const auto &row : local_histogram) {
        local_histogram_flat.insert(local_histogram_flat.end(), row.begin(), row.end());
    }

    // Ensure that the size of local_histogram_flat matches the total number of elements in local_histogram
    if (local_histogram_flat.size() != numBins * numBins) {
     //   std::cerr << "Error: Flattened array size does not match local_histogram size." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1); // Abort MPI
    }

    // Create a vector to store the global histogram
    std::vector<int> global_histogram_flat(numBins * numBins, 0);

    // Use MPI_Allreduce with the flattened array
    MPI_Allreduce(local_histogram_flat.data(), global_histogram_flat.data(), numBins * numBins, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // Reshape the global_histogram_flat into a 2D vector
    std::vector<std::vector<int>> global_histogram(numBins, std::vector<int>(numBins, 0));
    for (int i = 0; i < numBins; ++i) {
        for (int j = 0; j < numBins; ++j) {
global_histogram[i][j] = global_histogram_flat[i * numBins + j];       
 }
    }
#endif // MPI_PARALLEL

    if (Tmax==tmax_global)
     {
        std::ofstream P_dis;
        // Open the file stream outside the loop
        P_dis.open("T_Max.txt", std::ios::app);
      // P_dis << time <<" "<< tmax_global<<std::endl;
      P_dis << time <<" "<< tmax_global << " at position: (" << max_x << ", " << max_y << ", " << max_z << ")" <<std::endl;
      P_dis.close(); // Corrected the closing statement
        
        // Open the file stream for the current histogram file
    std::ofstream outFile(histogramFileName);
    for (int i = 0; i < numBins; ++i) {
        for (int j = 0; j < numBins; ++j) {
            outFile << global_histogram[i][j] << " ";
        }
        outFile << std::endl;
    }

        // Close the file stream
        outFile.close();
        // you can add the time as well:
    }
}


//----------------------------------------------------------------------------------------
//! \fn void JetInnerX1()
//  \brief Sets boundary condition on left X boundary (iib) for jet problem

/*void JetInnerX2_2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set primitive variables in inlet ghost zones
  for (int k=kl; k<=ku; ++k) {
    for (int i=il; i<=iu; ++i) {


        Real x = pco->x1v(i);
        Real y = pco->x2v(jl);
        Real z = pco->x3v(k);
        Real rad = std::sqrt(SQR(x - x_0) + SQR(z - z_0));  

        if (rad <= r_jet) {

          prim(IDN,k,jl,i) =  d_jet;
          prim(IVX,k,jl,i) =  vx_jet; 
          prim(IVY,k,jl,i) =  vy_jet; 
          prim(IVZ,k,jl,i) =  vz_jet; 
          prim(IPR,k,jl,i) =  p_jet;
          if (NSCALARS > 0){
            pmb->pscalars->s(0, k, jl, i) = d_jet;
            pmb->pscalars->r(0, k, jl, i) = 1.0;
          }
	} else if (rad> r_jet && rad<= (r_jet+wall_thick))
          {
          prim(IDN,k,jl,i) = d_amb;
          prim(IVX,k,jl,i) = 0.0;
          prim(IVY,k,jl,i) = 0.0;
          prim(IVZ,k,jl,i) = 0.0;
          prim(IPR,k,jl,i) = p_amb;
          if (NSCALARS > 0){
            pmb->pscalars->s(0, k,jl,i) = pmb->pscalars->s(0, k,jl,i);
            pmb->pscalars->r(0, k,jl,i) = pmb->pscalars->r(0, k,jl,i);
          }
          } 
          else {
          prim(IDN,k,jl,i) = d_amb;
          prim(IVX,k,jl,i) = vx_amb;
          prim(IVY,k,jl,i) = vy_amb;
          prim(IVZ,k,jl,i) = vz_amb;
          prim(IPR,k,jl,i) = p_amb;
          if (NSCALARS > 0){
           pmb->pscalars->s(0, k,jl,i) = pmb->pscalars->s(0, k,jl,i);
           pmb->pscalars->r(0, k,jl,i) = pmb->pscalars->r(0, k,jl,i);
          }
	}
      }
    }
  }*/




//----------------------------------------------------------------------------------------
//! \fn void JetInnerX1()
//  \brief Sets boundary condition on left X boundary (iib) for jet problem

void JetInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  std::default_random_engine generator;
  Real var = 3.33; 
  std::normal_distribution<Real> distribution (0.0, var);
  // set primitive variables in inlet ghost zones
//if (time==0)
  // JetInnerX2_2(pmb,pco,prim,b,time,dt,il,iu,jl,ju,kl,ku,NGHOST);

  for (int k=kl; k<=ku; ++k) {
    for (int i=il; i<=iu; ++i) {
      for (int j=1; j<=ngh; ++j) {
 
        Real x = pco->x1v(i);
        Real y = pco->x2v(j);
        Real z = pco->x3v(k);
//        std::cout << "x " << pco->x1v(i)<< std::endl;

        Real rad = std::sqrt(SQR(x - x_0) + SQR(z - z_0));  

        if (rad <= r_jet) {

          Real vx_random = distribution (generator);
          Real vz_random = distribution (generator);
          prim(IDN,k,jl-j,i) =  d_jet;
          prim(IVX,k,jl-j,i) =  vx_jet+vx_random;
          prim(IVY,k,jl-j,i) =  vy_jet; 
          prim(IVZ,k,jl-j,i) =  vz_jet+vz_random; 
          prim(IPR,k,jl-j,i) =  p_jet;
          if (NSCALARS > 0){
            pmb->pscalars->s(0, k, j, i) = d_jet;
            pmb->pscalars->r(0, k, j, i) = 1.0;
          }
	} /*else if (rad> r_jet && rad<= (r_jet+wall_thick))
          {
          prim(IDN,k,jl-j,i) = prim(IDN,k,jl,i);
          prim(IVX,k,jl-j,i) = -prim(IVX,k,jl+j-1,i)+2*0.0;;
          prim(IVY,k,jl-j,i) = -prim(IVY,k,jl+j-1,i)+2*0.0;
          prim(IVZ,k,jl-j,i) = -prim(IVZ,k,jl+j-1,i)+2*0.0;
          prim(IPR,k,jl-j,i) = prim(IPR,k,jl,i);
          if (NSCALARS > 0){
            pmb->pscalars->s(0, k,jl-j,i) = pmb->pscalars->s(0, k,jl,i);
            pmb->pscalars->r(0, k,jl-j,i) = pmb->pscalars->r(0, k,jl,i);
          }
          }*/ 

  else {
          prim(IDN,k,jl-j,i) = prim(IDN,k,jl,i);
          prim(IVX,k,jl-j,i) = prim(IVX,k,jl,i);
          prim(IVY,k,jl-j,i) = -prim(IVY,k,jl,i);
          prim(IVZ,k,jl-j,i) = prim(IVZ,k,jl,i);
          prim(IPR,k,jl-j,i) = prim(IPR,k,jl,i);
          if (NSCALARS > 0){
            pmb->pscalars->s(0, k,jl-j,i) = pmb->pscalars->s(0, k,jl,i);
            pmb->pscalars->r(0, k,jl-j,i) = pmb->pscalars->r(0, k,jl,i);
          }
        }
      }
    }
  }
}



//----------------------------------------------------------------------------------------
//! \fn int RefinementCondition(MeshBlock *pmb)
//! \brief refinement condition: maximum density and pressure curvature

/*int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps = 0.0;
  //Real V;
  int k = pmb->ks;
  for (int j=pmb->js; j<=pmb->je; j++) {
    for (int i=pmb->is; i<=pmb->ie; i++) {
      Real epsr = (std::abs(w(IDN,k,j,i+1) - 2.0*w(IDN,k,j,i) + w(IDN,k,j,i-1))
                  + std::abs(w(IDN,k,j+1,i) - 2.0*w(IDN,k,j,i) + w(IDN,k,j-1,i)))
                  /w(IDN,k,j,i);
      Real epsp = (std::abs(w(IPR,k,j,i+1) - 2.0*w(IPR,k,j,i) + w(IPR,k,j,i-1))
                  + std::abs(w(IPR,k,j+1,i) - 2.0*w(IPR,k,j,i) + w(IPR,k,j-1,i)))
                  /w(IPR,k,j,i);
      Real eps = std::max(epsr, epsp);
      maxeps = std::max(maxeps, eps);
//      Real V=w(IVX,k,j,i)+w(IVY,k,j,i)+w(IVZ,k,j,i);
     //  V = sqrt(SQR((IVX,k,j,i))+SQR(w(IVY,k,j,i))+SQR(w(IVZ,k,j,i)));
    }
  }
  // refine : curvature > 0.01
  //if (std::abs(V)> 0.0) return 1;
  if (maxeps > 0.01) return 1;
  // derefinement: curvature < 0.005
  if (maxeps < 0.005) return -1;
  // otherwise, stay
  return 0;
}*/

int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps = 0.0;
  int k = pmb->ks;
  for (int j=pmb->js; j<=pmb->je; j++) {
    for (int i=pmb->is; i<=pmb->ie; i++) {
      Real epsr = (std::abs(w(IDN,k,j,i+1) - 2.0*w(IDN,k,j,i) + w(IDN,k,j,i-1))
                  + std::abs(w(IDN,k,j+1,i) - 2.0*w(IDN,k,j,i) + w(IDN,k,j-1,i)))
                  /w(IDN,k,j,i);
      Real epsp = (std::abs(w(IPR,k,j,i+1) - 2.0*w(IPR,k,j,i) + w(IPR,k,j,i-1))
                  + std::abs(w(IPR,k,j+1,i) - 2.0*w(IPR,k,j,i) + w(IPR,k,j-1,i)))
                  /w(IPR,k,j,i);
      Real eps = std::max(epsr, epsp);
      maxeps = std::max(maxeps, eps);
    }
  }
  // refine : curvature > 0.01
  if (maxeps > 0.01) return 1;
  // derefinement: curvature < 0.005
  if (maxeps < 0.005) return 0;
  // otherwise, stay
  return 0;
}


void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  // AthenaArray<Real> vol_;

  // vol_.NewAthenaArray((ie - is) + 2 * NGHOST + 1);

  const bool f2 = pmy_mesh->f2;
  const bool f3 = pmy_mesh->f3;

  int il = is - 1;
  int iu = ie + 1;
  int jl, ju, kl, ku;
  Real area_p1, area;
  Real u_p1, u, v_p1, v;
  Real u_cross_a_p1, u_cross_a;
  Real v_cross_a_p1, v_cross_a;

  Real d_inf = pin->GetReal("problem", "djet");
  Real p_inf = pin->GetReal("problem", "pjet");
  Real u_inf = pin->GetReal("problem", "vxjet");

  if (!f2) // 1D
    jl = js, ju = je, kl = ks, ku = ke;
  else if (!f3) // 2D
    jl = js - 1, ju = je + 1, kl = ks, ku = ke;
  else // 3D
    jl = js - 1, ju = je + 1, kl = ks - 1, ku = ke + 1;

  for (int k = kl; k <= ku; ++k)
  {
    for (int j = jl; j <= ju; ++j)
    {
#pragma omp simd
      for (int i = il; i <= iu; ++i)
      {
        user_out_var(0, k, j, i) = phydro->w(IPR, k, j, i);
        user_out_var(1, k, j, i) = phydro->w(IDN, k, j, i);
        user_out_var(2, k, j, i) = phydro->w(IVX, k, j, i);
        user_out_var(3, k, j, i) = phydro->w(IVY, k, j, i);
        user_out_var(4, k, j, i) = phydro->w(IVZ, k, j, i);
        if (NSCALARS > 0) {
          Real Y = pscalars->r(0, k, j, i); // Y is mass fraction H2/air
          Real M_x = 1.0 / (((1-Y) / Mw_air) + (Y / Mw_fuel));
          Real R_gas = 8.314/M_x;
          user_out_var(5, k, j, i) = phydro->w(IPR, k, j, i) / phydro->w(IDN, k, j, i) / R_gas;
        }
        Real u2 = phydro->w(IVX, k, j, i) * phydro->w(IVX, k, j, i) + phydro->w(IVY, k, j, i) * phydro->w(IVY, k, j, i) + phydro->w(IVZ, k, j, i) * phydro->w(IVZ, k, j, i);
        Real c2 = phydro->w(IPR, k, j, i) / phydro->w(IDN, k, j, i) * peos->GetGamma();
        Real mach = sqrt(u2 / c2);
        user_out_var(6, k, j, i) = mach;
        if (NSCALARS > 0)
          user_out_var(7, k, j, i) = pscalars->r(0, k, j, i);
        Real den_x = 0, den_y = 0, den_z = 0;
        den_x = (phydro->w(IDN, k, j, i + 1) - phydro->w(IDN, k, j, i - 1)) / pcoord->dx1v(i);
        if (f2)
          den_y = (phydro->w(IDN, k, j + 1, i) - phydro->w(IDN, k, j - 1, i)) / pcoord->dx2v(j);
        if (f3)
          den_z = (phydro->w(IDN, k + 1, j, i) - phydro->w(IDN, k - 1, j, i)) / pcoord->dx3v(k);
        user_out_var(8, k, j, i) = sqrt(SQR(den_x) + SQR(den_y) + SQR(den_z)); // Schlieren image
        user_out_var(9, k, j, i) = phydro->U_mean(0, k, j, i); // velocity_x mean
        user_out_var(10, k, j, i) = phydro->U_mean(1, k, j, i); // velocity_y mean
        user_out_var(11, k, j, i) = phydro->U_mean(2, k, j, i); // velocity-Z mean
        user_out_var(12, k, j, i) = phydro->U_mean(3, k, j, i); // Temperature mean
        user_out_var(13, k, j, i) = phydro->U_mean(4, k, j, i); // massfraction
        user_out_var(14, k, j, i) = phydro->U_mean(5, k, j, i); // Density
        user_out_var(15, k, j, i) = phydro->U_mean(6, k, j, i); // mach
        user_out_var(16, k, j, i) = phydro->U_mean(7, k, j, i); // Mean U_r
        user_out_var(17, k, j, i) = phydro->U_fluc(0, k, j, i); // fluct_vx
        user_out_var(18, k, j, i) = phydro->U_fluc(1, k, j, i); // fluct_vy
        user_out_var(19, k, j, i) = phydro->U_fluc(2, k, j, i); // fluct_vz
        user_out_var(20, k, j, i) = phydro->U_fluc(3, k, j, i); // fluct_u_r
        user_out_var (21,k, j, i) = phydro->Rey_stress(0, k ,j , i); // Re_xx
        user_out_var (22,k, j, i) = phydro->Rey_stress(1, k ,j , i); // Re_yy
        user_out_var (23,k, j, i) = phydro->Rey_stress(2, k ,j , i); // Re_zz
        user_out_var (24,k, j, i) = phydro->Rey_stress(3, k ,j , i); // Re_xy
        user_out_var (25,k, j, i) = phydro->Rey_stress(4, k ,j , i); // Re_xz
        user_out_var (26,k, j, i) = phydro->Rey_stress(5, k ,j , i); // Re_yz
        user_out_var (27,k, j, i) = phydro->Rey_stress(6, k ,j , i); // Re_Ur*Ur
        user_out_var (28,k, j, i) = phydro->Rey_stress(7, k ,j , i); // Re_Ur*Ux        
       // user_out_var(9, k, j, i) = phydro->U_mean(0,k,j,i)/ phydro->U_mean(2,k,j,i); // 
       // user_out_var(10, k, j, i) = phydro->U_mean(1,k,j,i)/ phydro->U_mean(2,k,j,i); // Schlieren image
       // user_out_var(11, k, j, i) = phydro->U_mean(2,k,j,i);

      }
    }
  }
}
