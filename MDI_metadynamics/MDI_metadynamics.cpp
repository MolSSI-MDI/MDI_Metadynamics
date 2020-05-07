#include <iostream>
#include <mpi.h>
#include <stdexcept>
#include <string.h>
#include <iomanip>
#include <math.h>
#include <fstream>
#include "mdi.h"
#include "Distance.h"

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

  // Define collective variables 

  CollectiveVariable * colvar;
  colvar = new Distance(319, 320);

  double kcalmol_to_atomic;

  MDI_Conversion_Factor("kilocalorie_per_mol","atomic_unit_of_energy", &kcalmol_to_atomic);

  double angstrom_to_atomic;
  
  MDI_Conversion_Factor("angstrom","atomic_unit_of_length", &angstrom_to_atomic);

  double kcalmol_per_angstrom_to_atomic = kcalmol_to_atomic / angstrom_to_atomic;

  double width = 0.2 * angstrom_to_atomic; // Gaussian width of first collective variable.

//  double height = 0.02 * kcalmol_to_atomic; //Gaussian height of first collective variable.
  double height = 0.1 * kcalmol_to_atomic; //Gaussian height of first collective variable.

  const int total_steps = 100000000;  // Number of MD iterations. Note timestep = 2fs.

  const int tau_gaussian = 400; // Frequency of addition of Gaussians.

  double upper_restraint = 8.0 * angstrom_to_atomic;
  double k_restraint = 10 * kcalmol_per_angstrom_to_atomic; 
  bool verbose = false;

  const int total_gaussians = (total_steps >= tau_gaussian) ? total_steps / tau_gaussian : 1;
  
  array<double, total_gaussians> s_of_t = {0}; // value of collective variable at time t'.
  double s_of_x;

  double dVg_ds=0.0;

  double dg_ds=0.0;

  bool in_bounds = true;

  double arg;

  int current_gaussians = 0;

  double colvar_val = 0.0; 

  // Create output file
  ofstream output_file;
  output_file.open("output.dat", fstream::app);
  if (!output_file)
  {
	  cout << "Output file could not be opened" << endl;
	  return 1;
  }

  output_file << "#" << setw(9) << "Time_step";
  output_file << setw(20) << "Colvar";
  output_file << setw(20) << "Width";
  output_file << setw(20) << "Height";
  output_file << endl;

  // Connect to the engines
  MDI_Comm comm = MDI_COMM_NULL;
  MDI_Accept_Communicator(&comm);
  if (verbose)
  	cout << "Connected to engine." << endl;
 
  // Get engine name
  char* engine_name = new char[MDI_NAME_LENGTH];
  MDI_Send_Command("<NAME", comm); 
  MDI_Recv(engine_name, MDI_NAME_LENGTH, MDI_CHAR, comm);
  if (verbose)
  	cout << "Engine name: " << engine_name << endl;

  // Get number of atoms
  int natoms;
  MDI_Send_Command("<NATOMS", comm);
  MDI_Recv(&natoms, 1, MDI_INT, comm);

  double forces[3*natoms];

  // Initialize MD simulation
  MDI_Send_Command("@INIT_MD", comm);

  if (verbose)
	  cout << "MD simulation successfully initialized." << endl;

  for (int time_step = 0; time_step < total_steps + 1; time_step++) {

  	if (verbose)
   		cout << "Iteration: " << time_step << " out of " << total_steps << endl;

    // Get simulation box size
    double cell_size[9];

    MDI_Send_Command("<CELL", comm);
    MDI_Recv(&cell_size,9, MDI_DOUBLE, comm);

    array3d lo, hi, box_len;

    // NOTE: The following assumes that the cell vectors are orthogonal
    box_len[0] = cell_size[0];
    box_len[1] = cell_size[4];
    box_len[2] = cell_size[8];

    if (verbose) {
      cout << "  Read box size successfully. " << endl;
      cout << "     Box Size: " << box_len[0] << " " << box_len[1] << " " << box_len[2] << endl;
    }

    // Get current Cartesian coordinates
   
    double coords[3*natoms];

    MDI_Send_Command("<COORDS", comm);
    MDI_Recv(&coords, 3*natoms, MDI_DOUBLE, comm);

	if (verbose)
   		cout << "  Read coordinates successfully. " << endl;

    // This gets the current values of CVs and gradients
    colvar -> Evaluate(coords, natoms, box_len);
    colvar_val = colvar->Get_Value();

    // Update the bias function
    if (time_step % tau_gaussian == 0) {

        output_file << setw(10) << time_step;
        output_file << setw(20) << colvar_val / angstrom_to_atomic;
        output_file << setw(20) << width / angstrom_to_atomic;
        output_file << setw(20) << height / kcalmol_to_atomic;
        output_file << endl;

        s_of_t[current_gaussians] = colvar_val;
        current_gaussians++;
	}

    // Evaluate the derivative of Gaussians wrt to Cartesian Coordinates 
    dVg_ds = 0;
	for (int idx_t = 0; idx_t < current_gaussians; idx_t++) {
		s_of_x = colvar_val;
        arg = s_of_x - s_of_t[idx_t];
        dVg_ds = dVg_ds + Gaussian_derv(arg, width, height);
	}

	// Restraints
    if (colvar_val > upper_restraint)
         dVg_ds = k_restraint * (colvar_val - upper_restraint);
//     else if (colvar_val > 4.0)
//         dVg_ds += 200 * kcalmol_to_atomic * (colvar_val - 14.0 * angstrom_to_atomic);

	if (verbose)
	    cout << "  Evaluated gradients successfully. " << endl;

    // Set the forces

  MDI_Send_Command("@FORCES", comm);
  MDI_Send_Command("<FORCES", comm);
  MDI_Recv(&forces, 3*natoms, MDI_DOUBLE, comm);

  array<array3d, 2> delta_force;
  
  array<array3d, 2> ds_dr = colvar -> Get_Gradient(); // dimensionless 
  
  array2dint atoms_colvar = colvar->Get_Atoms();
  
  for (int idx_atom = 0; idx_atom < 2; idx_atom++) {
  
	cout << atoms_colvar[idx_atom] << endl;
    for (int idx_dir = 0; idx_dir < 3; idx_dir++) {
 
      delta_force[idx_atom][idx_dir] = dVg_ds * ds_dr[idx_atom][idx_dir];
  
      forces[3 * atoms_colvar[idx_atom]+idx_dir] -= delta_force[idx_atom][idx_dir];
  
   	}
  }
     
    
    MDI_Send_Command(">FORCES", comm);
   MDI_Send(&forces, 3*natoms, MDI_DOUBLE, comm);
  
	if (verbose)
	    cout << "Set biased forces successfully." <<  endl;

    MDI_Send_Command("@COORDS", comm);

	if (verbose)
	    cout << "Moved to next step." <<  endl;


  } // Main MD loop


   // Send the "EXIT" command to each of the engines
   MDI_Send_Command("EXIT", comm);

   // Synchronize all MPI ranks
   MPI_Barrier(world_comm);

   output_file.close();

   cout << "Finished" << endl;

   return 0;
 }
