#include <iostream>
#include <mpi.h>
#include <stdexcept>
#include <string.h>
#include "dihedral.h"
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
  int natoms;
  double coords[3*natoms];
  double forces[3*natoms];
  char* engine_name = new char[MDI_NAME_LENGTH];

  // Connect to the engines
  MDI_Comm comm = MDI_NULL_COMM;
  comm = MDI_Accept_Communicator();
 
  // Get engine name
  MDI_Send_Command("<NAME", comm); 
  MDI_Recv(engine_name, MDI_NAME_LENGTH, MDI_CHAR, comm);
  cout << "Engine name: " << engine_name << endl;

  // Get number of atoms
  MDI_Send_Command("<NATOMS", comm);
  MDI_Recv(&natoms, 1, MDI_INT, comm);
  
  cout << "Number of atoms: " << natoms << endl;
  
  



  // Perform the simulation
  // <YOUR CODE GOES HERE>

  // Send the "EXIT" command to each of the engines
  // <YOUR CODE GOES HERE>

  // Synchronize all MPI ranks
  MDI_Send_Command("EXIT", comm);
  MPI_Barrier(world_comm);

  return 0;
}
