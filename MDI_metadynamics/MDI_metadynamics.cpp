#include <iostream>
#include <mpi.h>
#include <stdexcept>
#include <string.h>
#include <iomanip>
#include <math.h>
#include "Dihedral.h"
#include <fstream>
#include "mdi.h"

using namespace std;

int main(int argc, char **argv) {

  // Initialize the MPI environment
  MPI_Comm world_comm;
  MPI_Init(&argc, &argv);

  // Read through all the command line options
  int iarg = 1;
  bool initialized_mdi = false;
  while ( iarg < argc ) {

    if ( strcmp(argv[iarg],"-mdi") == 0 ) {

      // Ensure that the argument to the -mdi option was provided
      if ( argc-iarg < 2 ) {
	throw runtime_error("The -mdi argument was not provided.");
      }

      // Initialize the MDI Library
      world_comm = MPI_COMM_WORLD;
      int ret = MDI_Init(argv[iarg+1], &world_comm);
      if ( ret != 0 ) {
	throw runtime_error("The MDI library was not initialized correctly.");
      }
      initialized_mdi = true;
      iarg += 2;

    }
    else {
      throw runtime_error("Unrecognized option.");
    }

  }
  if ( not initialized_mdi ) {
    throw runtime_error("The -mdi command line option was not provided.");
  }

  // Metadynamics settings
  // TODO: change all arrays to vectors

  // Define collective variables 
  const int num_colvars= 1;

  array<CollectiveVariable *, num_colvars> colvars;
  colvars[0] = new Dihedral(4, 6, 8, 14);
 
  std::array<double, num_colvars> width;
  width[0] = 0.2; // radians. Gaussian width of first collective variable.

  std::array<double, num_colvars> height;
  height[0] = 0.1; // kcal/mol. Gaussian height of first collective variable.
  
  const int total_steps = 1000;  // Number of MD iterations. Note timestep = 2fs.

  const int tau_gaussian = 500; // Frequency of addition of Gaussians.

  const int output_freq = 5;

  int current_gaussians = 0; // Number of Gaussians added so far.

  const int total_gaussians = total_steps / tau_gaussian; 

  array<array2d, total_gaussians> s_of_t; // value of collective variable at time t'.
  

  // Create output file
  ofstream output_file;
  output_file.open("output.dat", fstream::app);
  if (!output_file)
  {
	  cout << "Output file could not be opened" << endl;
	  return 1;
  }

  for (int idx_cv = 0; idx_cv < num_colvars; idx_cv++) {
    output_file << setw(14) << "# Colvar " << idx_cv <<  endl;
  }

  // Create bias output file
  ofstream bias_file;
  bias_file.open("bias.dat", fstream::app);
  if (!bias_file)
  {
	  cout << "Bias file could not be opened" << endl;
	  return 1;
  }

  // Connect to the engines
  MDI_Comm comm = MDI_NULL_COMM;
  comm = MDI_Accept_Communicator();
 
  // Get engine name
  char* engine_name = new char[MDI_NAME_LENGTH];
  MDI_Send_Command("<NAME", comm); 
  MDI_Recv(engine_name, MDI_NAME_LENGTH, MDI_CHAR, comm);

  // Get number of atoms
  int natoms;
  MDI_Send_Command("<NATOMS", comm);
  MDI_Recv(&natoms, 1, MDI_INT, comm);
  
  // Initialize MD simulation
  MDI_Send_Command("MD_INIT", comm);

  cout << "MD simulation successfully initialized." << endl;

  for (int time_step = 0; time_step < total_steps; time_step++) {

    cout << "Iteration: " << time_step << " out of " << total_steps << endl;

    // Get simulation box size
    double cell_size[9];

    MDI_Send_Command("<CELL", comm);
    MDI_Recv(&cell_size,9, MDI_DOUBLE, comm);

    array3d lo, hi, box_len;

    lo[0] = cell_size[0];
    lo[1] = cell_size[1];
    lo[2] = cell_size[2];

    hi[0] = cell_size[3];
    hi[1] = cell_size[4];
    hi[2] = cell_size[5];

    box_len[0] = hi[0] - lo[0];
    box_len[1] = hi[1] - lo[1];
    box_len[2] = hi[2] - lo[2];

    cout << "  Read box size successfully. " << endl;

    // Get current Cartesian coordinates
   
    double coords[3*natoms];

    MDI_Send_Command("<COORDS", comm);
    MDI_Recv(&coords, 3*natoms, MDI_DOUBLE, comm);
    cout << "  Read coordinates successfully. " << endl;

    // This gets the current values of CVs and gradients
    for (int idx_cv = 0; idx_cv < num_colvars; idx_cv++) {
      colvars[idx_cv]->Evaluate(coords, natoms, box_len);
      
      if (time_step % output_freq == 0) {
      	output_file << setw(15) << colvars[idx_cv]->Get_Value();
      }
    }
      if (time_step % output_freq == 0) {
    output_file << endl; }

  

    // Update the bias function
    if (time_step % tau_gaussian == 0) {


      // Get the current value of CVs and save it in s_of_t array
      for (int idx_cv = 0; idx_cv < num_colvars; idx_cv++) {
      	s_of_t[idx_cv][current_gaussians] = colvars[idx_cv]->Get_Value();
      }

      double vg = 0.0;

      for (int idx_t = 0; idx_t <= current_gaussians; idx_t++) {

	double g = 0.0;

        for (int idx_cv = 0; idx_cv < num_colvars; idx_cv++) {

          double s_of_x = colvars[idx_cv]->Get_Value();
          double arg = s_of_x - s_of_t[idx_cv][idx_t];
          g = g + Gaussian(arg, width[idx_cv], height[idx_cv]);

        }

        vg = vg + g;

      }

      bias_file << "# Iteration: " << time_step << endl;
      bias_file << "# NGaussian: " << current_gaussians << endl;
      
      for (int idx_t = 0; idx_t < s_of_t.size(); idx_t++) {

	      bias_file << vg << endl;
		      
      }

      current_gaussians += 1;
      if (current_gaussians > total_gaussians) {
        return 1;
      }

    cout << "  Computed bias successfully. " << endl;
    }


    // Evaluate the CVs and their gradients wrt to Cartesian coordinates

    double dVg_ds = 0;

    for (int idx_t = 0; idx_t <= current_gaussians-1; idx_t++){
     
      double dg_ds = 0.0;

      for (int idx_cv = 0; idx_cv < num_colvars; idx_cv++) {

        double s_of_x = colvars[idx_cv]->Get_Value();
        double arg = s_of_x - s_of_t[idx_cv][idx_t];
        dg_ds = dg_ds + Gaussian_derv(arg, width[idx_cv], height[idx_cv]);
//	cout << s_of_x << "   " << s_of_t[idx_cv][idx_t] << "  " << Gaussian_derv(arg, width[idx_cv], height[idx_cv]) << "  " <<  dg_ds << endl;
      }

      dVg_ds = dVg_ds - dg_ds;

    }

    cout << "  Evaluated gradients successfully. " << endl;

    // Set the forces
   
    double forces[3*natoms];
    MDI_Send_Command("@FORCES", comm);
    MDI_Send_Command("<FORCES", comm);
    MDI_Recv(&forces, 3*natoms, MDI_DOUBLE, comm);

    array<array3d, 4> delta_force;

    for (int idx_cv = 0; idx_cv < num_colvars; idx_cv++) {
    
      array<array3d, 4> ds_dr = colvars[idx_cv] -> Get_Gradient();
      
      array4dint atoms_colvar = colvars[idx_cv]->Get_Atoms();

      for (int idx_atom = 0; idx_atom < 4; idx_atom++) {
    
        for (int idx_dir = 0; idx_dir < 3; idx_dir++) {

          delta_force[idx_atom][idx_dir] = dVg_ds * ds_dr[idx_atom][idx_dir];

	  forces[3 * atoms_colvar[idx_atom]+idx_dir] -= delta_force[idx_atom][idx_dir];

	}
      }
     
    }
    
    MDI_Send_Command(">FORCES", comm);
    MDI_Send(&forces, 3*natoms, MDI_DOUBLE, comm);
   
    cout << "Set biased forces successfully." <<  endl;
    MDI_Send_Command("@COORDS", comm);
    cout << "Moved to next step." <<  endl;


  } // Main MD loop


   // Send the "EXIT" command to each of the engines
   MDI_Send_Command("EXIT", comm);

   // Synchronize all MPI ranks
   MPI_Barrier(world_comm);

   output_file.close();
   bias_file.close();

   cout << "Finished" << endl;
   return 0;
 }
