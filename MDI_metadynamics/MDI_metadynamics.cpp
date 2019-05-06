#include <iostream>
#include <mpi.h>
#include <stdexcept>
#include <string.h>
#include "CollectiveVariable.h"
extern "C" {
#include "mdi.h"
}

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

  // Metadynamics
  
  int niterations = 10;  // Number of MD iterations
  double width_colvar1= 10.0;  // Gaussian width of dihedral phi 
  double width_colvar2= 20.0;  // Gaussian width of dihedral theta
  int freq = 10; // Frequency of addition of Gaussians
  
  // Define collective variables
  CollectiveVariable * colvar = new Dihedral(0, 1, 2, 3);

  // Connect to the engines
  MDI_Comm comm = MDI_NULL_COMM;
  comm = MDI_Accept_Communicator();
 
  // Get engine name
  char* engine_name = new char[MDI_NAME_LENGTH];
  MDI_Send_Command("<NAME", comm); 
  MDI_Recv(engine_name, MDI_NAME_LENGTH, MDI_CHAR, comm);
  cout << "Engine name: " << engine_name << endl;

  // Get number of atoms
  int natoms;
  MDI_Send_Command("<NATOMS", comm);
  MDI_Recv(&natoms, 1, MDI_INT, comm);
  
  cout << "Number of atoms: " << natoms << endl;
  
  // Initialize MD simulation
  MDI_Send_Command("MD_INIT", comm);


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


  // Get current Cartesian coordinates
  double coords[3*natoms];
  double forces[3*natoms];

  MDI_Send_Command("<COORDS", comm);
  MDI_Recv(&coords, 3*natoms, MDI_DOUBLE, comm);
  colvar->Evaluate(coords, natoms, box_len);
  colvar->Compute_gradient(coords, natoms, box_len);

  colvar -> Get_Value();
  std::array<array3d, 4> grad = colvar -> Get_Gradient();

  //MDI_Send_Command("@FORCES", comm);

  //MDI_Send_Command("<FORCES", comm);
  //MDI_Send(&forces, 3*natoms, MDI_DOUBLE, comm);

  //MDI_Send_Command("@COORDS", comm);



  // Perform the simulation
  // <YOUR CODE GOES HERE>

  // Send the "EXIT" command to each of the engines
  // <YOUR CODE GOES HERE>

  // Synchronize all MPI ranks
  MDI_Send_Command("EXIT", comm);
  MPI_Barrier(world_comm);

  return 0;
}
