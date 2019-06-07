#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <omp.h>

/**
  compile : mpicc projet.c -o projet
  run : mpirun --oversubscribe -np 4 ./projet data/my_a_8 data/my_b_8
**/

/*************************************************************************************************************************
*** Projet du cours Programmation parallèle. Il s'agit de calculer le produit matriceXmatrice de manière distribuée.
*** Utilisation de la bibliothèque OPENMP pour distribuer les processus.
*** @Auteur : Aourinmouche Soufiane.
*** Filière : SI.
*** Date rendu : N/D.
****************************************************************************************************************************/

/************************************************** La structure du tableau **************************************************/
struct tablo {
    int * tab;
    int size;
    int rows;
    int cols;
};

/**
  Cette fonction retourne le nombre de valeurs par ligne/colonne
  C'est la racine carré de la taille, puisque la matrice est carrée
**/
int getNbColumns_Rows(struct tablo * tmp) {
    return (int)(sqrt(tmp->size));
}

void printArray(struct tablo * tmp) {
    int size = tmp->size;
    int i;
    for (i = 0; i < size; ++i) {
        printf("%i ", tmp->tab[i]);
    }
    printf("\n");
}

void printArrayAsMatrix(struct tablo * tmp) {
    int size = tmp->size;
    int nb_rows = getNbColumns_Rows(tmp);

    int i;
    for (i = 0; i < size; ++i) {
        printf("%i ", tmp->tab[i]);
        if ( i%nb_rows == nb_rows-1 ) {
            printf("\n");
        }
    }
}

struct tablo * allocateTablo(int size) {
    struct tablo * tmp = malloc(sizeof(struct tablo));
    tmp->size = size;
    tmp->tab = malloc(size*sizeof(int));

#pragma omp parallel for
    for ( int i = 0 ; i < tmp->size ; i ++) {
        tmp->tab[i] = 0;
    }

    return tmp;
}

// retourne le nombre de colone dans une matrice
// tab est la matrice globale
int get_nb_rows_offset(struct tablo * tab) {
    return (int)(sqrt(tab->size));
}

int get_size_slice(struct tablo * tab, int numprocs) {
    return (tab->size/numprocs);
}

int get_slice_width(struct tablo * tab, int numprocs) {
    return (getNbColumns_Rows(tab)/numprocs);
}


int get_nb_cols(struct tablo * tab, int numprocs) {
    return ( ((int)sqrt(tab->size)) / numprocs );
}



// retourne le nombre de lignes dans une matrice (sous-matrice)
// tab est la matrice globale
int get_nb_rows(struct tablo * tab, int numprocs) {
    return (get_nb_rows_offset(tab) / numprocs);
}


/************************************************** Utiles **************************************************/

/*
  Cette fonction permet de remplir le tableau source par les entiers récupérés du fichier d'entrée
*/
struct tablo * getArray(FILE * file) {
    int k =-1;
    int i = 0;
    struct tablo * input = allocateTablo(2);
    while ( (fscanf(file , "%d", &k) != EOF )) {
        //printf("___%d",k);

        if ( i > input->size - 1 ) {
            //input->size *= 2;
            input->size ++;
            input->tab = realloc(input->tab, input->size * sizeof(int));
        }
        input->tab[i] = k;
        i ++;
    }

    //printf(".\n");
    return input;
}


/**
  Cette fonction retourne la tranche "numbre_slice"  dans les lignes de la  matrice.
  Cette tranche sera envoyé au process correspondant pour la traiter dans le Scatter.
**/
struct tablo * get_slice_row(struct tablo * tableau, int numbre_slice, int numbre_process) {

    assert(numbre_slice <= numbre_process);

    int size = tableau->size;
    int size_slice = size/numbre_process;
    struct tablo * slice = allocateTablo(size_slice);

    //printf("size = %d , element_by_row = %d ,  slice size = %d\n",tableau->size, element_by_row, size_slice);

    int i = 0;


    for ( i = 0 ; i < size && numbre_slice > 1 ; i += size_slice ) {
        numbre_slice --;
    }

    //printf(" i =  %d\n",i);

#pragma omp parallel for
    for ( int j = i ; j < i+size_slice ; j ++ ) {
        //printf("__ j'ajoute %d dans case %d \n",tableau->tab[j], j);
        slice->tab[j-i] = tableau->tab[j];
    }

    return slice;
}

/**
  Cette fonction retourne la tranche "numbre_slice"  dans les colones de la  matrice.
  Cette tranche sera envoyé au process correspondant pour la traiter dans le Scatter.
**/
struct tablo * get_slice_col(struct tablo * tableau, int numbre_slice, int numbre_process) {

    assert(numbre_slice <= numbre_process);

    int element_by_col = getNbColumns_Rows(tableau);
    int size = tableau->size;
    int size_slice = size/numbre_process;
    int slice_width = getNbColumns_Rows(tableau)/numbre_process;
    struct tablo * slice = allocateTablo(size_slice);

    int j = 0;

    for ( int k = 0 ; k < slice_width ; k ++ ) {
        for ( int i = (numbre_slice-1)*slice_width ; i < size ; i += element_by_col ) {
            slice->tab[j] = tableau->tab[i+k];
            j++;
        }
    }


    return slice;
}









void send_it(int * array_to_send, int size_to_send, int i) {

    MPI_Send(array_to_send, size_to_send, MPI_INT, i, 0, MPI_COMM_WORLD);

}


int * receive_it(int sender_rank, int size_to_receive) {

    //struct tab * A = allocateTablo(size_to_receive);
    int * tablo = malloc(sizeof(int)*size_to_receive);

    MPI_Recv(tablo, size_to_receive, MPI_INT, sender_rank, 0, MPI_COMM_WORLD, 0);

    return tablo;

}



/***********************************************************

  get_line = validée ! marche pour matrice et sous matrice
  get_column = marche ! pour sous matrice
**************************************************/

int * get_line(struct tablo * tableau, int num_line, int offset_rows) {


    int * line = malloc(sizeof(int) * offset_rows);

#pragma omp parallel for
    for ( int i = 0 ; i < offset_rows ; i ++  ) {
        line[i] = tableau->tab[(num_line*offset_rows) + i];
    }

    return line;
}


int * get_column(struct tablo * tableau, int num_col, int offset_cols) {

    int nb_values_in_col = tableau->size / offset_cols;

    int * col = malloc(sizeof(int) * nb_values_in_col);


#pragma omp parallel for
    for ( int i = 0 ; i < nb_values_in_col ; i ++ ) {
        col[i] = tableau->tab[((num_col*nb_values_in_col) +i)];
    }

    return col;
}


int get_product_two_vectors(int * v1, int * v2, int size) {
    int sum = 0;

    for ( int i = 0; i < size ; i ++) {
        sum += v1[i]*v2[i];
    }
    return sum;
}


struct tablo * get_matricial_product_3(struct tablo * A, struct tablo * B, int numprocs, int offset_lines, int offset_cols) {


    int product_content_size = A->rows * B->cols;

    struct tablo * product = allocateTablo(product_content_size);

    int * line_A;
    int * col_B;

#pragma omp parallel for
    for ( int i = 0 ; i < A->rows ; i ++ ) {

        for ( int j = 0 ; j < B->cols ; j ++ ) {

            line_A = get_line(A, i, offset_lines);
            col_B = get_column(B, j, offset_cols);


            product->tab[ (i*A->rows) + j] = get_product_two_vectors(line_A, col_B, A->cols);

            free(line_A);
            free(col_B);
        }

    }



    return product;
}


int getN(char * filename){
    FILE *fp = fopen(filename, "r"); //read-only
    char c;
    int p=1;
    if(!fp) {
        fprintf(stderr, "Error while processing the file \n");
        exit(EXIT_FAILURE);
    }
    // Get n
    while(fscanf(fp, "%c", &c) && (c != '\n') && (!feof(fp)) ){
        if(c == ' ') {
            p++;
        }
    }
    return p;
}



void print_Matrice(struct tablo * M) {
    for ( int i = 0 ; i < M->size ; i ++ ) {
        printf("%d  ",M->tab[i]);
        if ( (i != 0) && ((i+1)%M->cols) == 0 ) {
            printf("\n");
        }
    }
}


/*void set_matrice_in_indice(tablo * matrice_to_set, tablo * matrice_globale, int indice) {
  int offset = matrice->
}*/





int main(int argc, char ** argv) {

    // MPI attributs
    int rank, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    int number_amount;
    int dimension;

    if ( argc < 3 ) {
        printf("Erreur, argument manquant : ./aourinmouche.c A_file.txt B_file.txt\n");
        exit(1);
    }

    FILE * file_A = fopen(argv[1], "r");

    if ( file_A == NULL ) {
        printf("Erreur ouverture du fichier\n");
        exit(1);
    }

    FILE * file_B = fopen(argv[2], "r");

    if ( file_B == NULL ) {
        printf("Erreur ouvertude du fichier\n");
        exit(1);
    }

    dimension = getN(argv[1]);


    if ( rank == 0 ) {  // First one

        struct tablo * input_A = getArray(file_A);
        struct tablo * input_B = getArray(file_B);

        printf("************************** data **************************\n");
        printf("rows :\n");
        printArray(input_A);
        printf("cols :\n");
        printArray(input_B);
        printf("**************************      **************************\n\n\n");





        // Store own rows and cols
        struct tablo * own_cols;
        struct tablo * own_rows;


#pragma omp parallel sections
        {
#pragma omp section
            {
                own_rows = allocateTablo(input_A->size/numprocs);
                own_rows = get_slice_row(input_A, 0, numprocs);
            }

#pragma omp section
            {
                own_cols = allocateTablo(input_B->size/numprocs);
                own_cols = get_slice_col(input_B, 1, numprocs);
            }

        }


        own_rows->cols = dimension;
        own_rows->rows = own_rows->size/dimension;
        //printf("(rows %d, (%d,%d))\n", rank, own_rows->rows, own_rows->cols);
        printArray(own_rows);
        //print_Matrice(own_rows);


        own_cols->rows = dimension;
        own_cols->cols = own_cols->size/dimension;
        //printf("(cols %d, (%d,%d))\n", rank, own_cols->rows, own_cols->cols);
        printArray(own_cols);
        //print_Matrice(own_cols);




        struct tablo * product = get_matricial_product_3(own_rows, own_cols, numprocs, own_rows->cols, own_cols->cols);
        //printf("(produit %d)\n\n", rank);
        printArray(product);
        //print_Matrice(product);



        // Send rows
#pragma omp parallel for
        for ( int i = 1 ; i < numprocs ; i ++ ) {

            struct tablo * line_slices = allocateTablo(input_A->size/numprocs);
            line_slices = get_slice_row(input_A, i+1, numprocs);

            send_it(line_slices->tab, line_slices->size, i);



            free(line_slices->tab);
            free(line_slices);
        }


        // Send columns
#pragma omp parallel for
        for ( int i = 1 ; i < numprocs ; i ++ ) {

            struct tablo * col_slices = allocateTablo(input_B->size/numprocs);
            col_slices = get_slice_col(input_B, i+1, numprocs);

            send_it(col_slices->tab, col_slices->size, i);



            free(col_slices->tab);
            free(col_slices);
        }


        // Firstly, send zero-initialized array for each process
        struct tablo * initial_array = allocateTablo(dimension * (dimension/numprocs));
        initial_array->cols = own_cols->cols;
        initial_array->rows = own_rows->rows;


#pragma omp parallel for
        for ( int i = 1 ; i < numprocs ; i ++ ) {
            send_it(initial_array->tab, dimension * (dimension/numprocs), i);
        }

        //printf("(initial__________ %d, %d:(%d,%d))\n", rank, initial_array->size, initial_array->rows, initial_array->cols);
        //printf("dimension : %d * %d \n",initial_array->rows, initial_array->cols);
        //printArray(initial_array);
        //print_Matrice(initial_array);




        // Compute product
        //struct tablo * product = get_matricial_product_3(own_rows, own_cols, get, int numprocs, int offset_lines, int offset_cols)

        // Some useful space for lisibility
        //printf("\n\n\n\n\n\n\n");


        // Free dynamically allocated memory
        free(initial_array->tab);
        free(initial_array);

        free(input_A->tab);
        free(input_A);

        free(input_B->tab);
        free(input_B);

        free(own_rows->tab);
        free(own_rows);

        free(own_cols->tab);
        free(own_cols);
    }
    else {  // others

        int indice = rank;

        /*        receive rows              */
        // Predict size to receive for B
        MPI_Probe(0, 0, MPI_COMM_WORLD, &status);

        // When probe returns, the status object has the size and other
        // attributes of the incoming message. Get the message size
        MPI_Get_count(&status, MPI_INT, &number_amount);

        // Receive B
        struct tablo * my_part_lines = allocateTablo(number_amount);
        free(my_part_lines->tab);
        my_part_lines->tab = receive_it(0, my_part_lines->size);


        my_part_lines->cols = dimension;
        my_part_lines->rows = my_part_lines->size / my_part_lines->cols;
        //printf("(rows %d, (%d,%d) )\n", rank, my_part_lines->rows, my_part_lines->cols);
        printArray(my_part_lines);
        //print_Matrice(my_part_lines);




        /*        receive cols        */
        // Predict size to receive for B
        MPI_Probe(0, 0, MPI_COMM_WORLD, &status);

        // When probe returns, the status object has the size and other
        // attributes of the incoming message. Get the message size
        MPI_Get_count(&status, MPI_INT, &number_amount);

        // Receive B
        struct tablo * my_part_cols = allocateTablo(number_amount);
        free(my_part_cols->tab);
        my_part_cols->tab = receive_it(0, my_part_cols->size);


        my_part_cols->rows = dimension;
        my_part_cols->cols = my_part_cols->size / my_part_cols->rows;
        //printf("(cols %d,  (%d,%d))\n", rank, my_part_cols->rows, my_part_cols->cols);
        printArray(my_part_cols);
        //print_Matrice(my_part_cols);



        /*      receive initial array     */
        // Predict size to receive for B
        MPI_Probe(0, 0, MPI_COMM_WORLD, &status);

        // When probe returns, the status object has the size and other
        // attributes of the incoming message. Get the message size
        MPI_Get_count(&status, MPI_INT, &number_amount);

        //struct tablo * initial_array = allocateTablo(number_amount);
        struct tablo * initial_array = allocateTablo(dimension * (dimension/numprocs));
        free(initial_array->tab);
        initial_array->tab = receive_it(0, number_amount);
        initial_array->cols = my_part_cols->cols;
        initial_array->rows = my_part_lines->rows;

//    printf("(initial %d,  (%d,%d))\n", rank, initial_array->rows, initial_array->cols);
        //printf("dimension_S : %d * %d \n",initial_array->rows, initial_array->cols);
        //printArray(initial_array);

        //print_Matrice(initial_array);




        struct tablo * product = get_matricial_product_3(my_part_lines, my_part_cols, numprocs, my_part_lines->cols, my_part_cols->cols);
        //  printf("(produit %d) doit etre dans %d\n", rank, indice);
        indice ++;
        indice = indice % numprocs;

        printArray(product);

        //print_Matrice(product);


        // Some useful space for lisibility
        // printf("\n\n\n\n\n\n\n");



        // Free dynamically allocated memory


        free(initial_array->tab);
        free(initial_array);

        free(my_part_cols->tab);
        free(my_part_cols);

        free(my_part_lines->tab);
        free(my_part_lines);


        //free(product->tab);
        //free(product);
    }





    MPI_Finalize();


    return 0;
}