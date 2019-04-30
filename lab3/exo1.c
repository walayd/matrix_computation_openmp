#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){

	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  	// int partner_rank = (rank + 1) % 2;


	MPI_Status status;
	int number_databuffer;
	
	void one_to_one(){
		// printf("Hello, World!, rank : %d \n", rank);
		if (rank == 0) {
			number_databuffer = 100; // data that will be send to rank 1
			MPI_Send(&number_databuffer, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			printf("process rank %d send number %d to process rank 1\n", rank, number_databuffer);
		}
		if (rank == 1) {
			MPI_Recv(&number_databuffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("process rank %d received number %d from process rank 0\n", rank, number_databuffer);
		}
	}

	void one_to_many(){
		// printf("Hello, World!, rank : %d \n", rank);
		if (rank == 0) {
			number_databuffer = 100; // data that will be send to rank 1
			for(int i = 1; i < numprocs; i++){
				MPI_Send(&number_databuffer, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			} 
			printf("process rank %d send number %d to process rank 1\n", rank, number_databuffer);
		}
		else {
			MPI_Recv(&number_databuffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("process rank %d received number %d from process rank 0\n", rank, number_databuffer);
		}
	}

	void one_to_many_MPI_ANY_SOURCE(){
		// printf("Hello, World!, rank : %d \n", rank);
		if (rank == 0) {
			number_databuffer = 100; // data that will be send to rank 1
			for(int i = 1; i < numprocs; i++){
				MPI_Send(&number_databuffer, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			} 
			printf("process rank %d send number %d to process rank 1\n", rank, number_databuffer);
		}
		else {
			MPI_Recv(&number_databuffer, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			printf("process rank %d received number %d from process rank %d\n", rank, number_databuffer, status.MPI_SOURCE);
		}
	}

	void ring(){
		// ici ça marche quand on fait tous send en même temps
		// parce que Open MPI a un buffer interne, tant qu'on dépasse pas 
		// la taille du buffer, le send marche, il dira qu'il va l'envoyer..(vu que c'est des petit int)
		// mais si on envoit de gros tableaux, style 4k, le send devient blocant tant qu'on a pas recieve
		// astuce pour éviter ça : deux voisins ne font pas de send/recieve en même temps
		
		if (rank != 0) {
			MPI_Recv(&number_databuffer, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
			printf("process rank %d received number %d from process rank %d\n", rank, number_databuffer, status.MPI_SOURCE);
		}
		else {
			number_databuffer = 100;
		}
		MPI_Send(&number_databuffer, 1, MPI_INT, (rank + 1) % numprocs, 0, MPI_COMM_WORLD);
		printf("process rank %d send number %d to process rank 1\n", rank, number_databuffer);
		if (rank == 0) {
		    MPI_Recv(&number_databuffer, 1, MPI_INT, numprocs - 1, 0, MPI_COMM_WORLD, &status);
		    printf("process rank %d received number %d from process rank %d\n", rank, number_databuffer, status.MPI_SOURCE);
		}

	}
	// one_to_one();
	// one_to_many();
	// one_to_many_MPI_ANY_SOURCE();
	ring();

	MPI_Finalize();

	return 0;
}
