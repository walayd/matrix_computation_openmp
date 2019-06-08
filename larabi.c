#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <assert.h>

struct tablo {
    int *tab;
    int size;
    int nb_cols;
    int nb_rows;
};

struct tablo *allocateTablo(int size) {
    struct tablo *tmp = malloc(sizeof(struct tablo));
    tmp->size = size;
    tmp->nb_cols = (int) sqrt(size);
    tmp->nb_rows = (int) sqrt(size);
    tmp->tab = malloc(size * sizeof(int));
#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        tmp->tab[i] = 0;
    }
    return tmp;
}

struct tablo *reallocateTabloTwice(struct tablo *tmp) {
    int oldsize = tmp->size;
    tmp->size = oldsize * 2;
    tmp->tab = realloc(tmp->tab, tmp->size * sizeof(int));

    for (int i = oldsize; i < tmp->size; i++) {
        tmp->tab[i] = 0;
    }
    return tmp;
}

void printArray(struct tablo *tmp) {
    // printf("---- Array of size %i ---- \n", tmp->size);
    int size = tmp->size;
    int i;
    for (i = 0; i < size; ++i) {
        printf("%i ", tmp->tab[i]);
    }
    printf("\n");
}

int getNbColsWhenRead(char filename[]) {
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

int getNbRowsWhenRead(char filename[]) {
    printf("get nb rows when read ... \n");
    int nbRows = 0;
    FILE *fp;
    char c;

    fp = fopen(filename, "r");

    if (fp == NULL) {
        printf("Opening file error\n");
        return -1;
    }

    // find length of first line
    // lines end in "\n", but some malformed text files may not have this char at all
    // and whole file contents will be considered as the first line
    while ((c = fgetc(fp)) != EOF) {
        if (c == '\n') {
            nbRows++;
        }
    }

    printf("nb rows is: %d\n", nbRows);

    // seek to beginning of file
    fseek(fp, 0, SEEK_SET);
    fclose(fp);
    // allocate memory for size of first line (len)
    return nbRows;
}

int calculateMatrixCarreLength(char filename[]) {
    int nbCols = getNbColsWhenRead(filename);
    return nbCols * nbCols;
}

int calculateMatrixNotCarreLength(char filename[]) {
    int nbRow = getNbRowsWhenRead(filename);
    int nbCol = getNbColsWhenRead(filename);
    return nbRow * nbCol;
}

void gettingMatrixDataFromFile(char filename[], struct tablo *matrix) {
    FILE *fp = fopen(filename, "r"); //read-only
    int value = 0;
    int count = 0;
    if (!fp) {
        fprintf(stderr, "Error while processing the file \n");
        exit(EXIT_FAILURE);
    }
    //Going back to the start point
    fseek(fp, 0, SEEK_SET);
    // printf("%d lignes, %d colonnes \n", matrix->nb_rows , matrix->nb_cols);
    if (!matrix->tab) {
        fprintf(stderr, "Error while allocating for %d size \n", matrix->size);
        exit(EXIT_FAILURE);
    }
    //Store values in tab
    while (fscanf(fp, "%d", &value) && (!feof(fp))) {
        matrix->tab[count] = value;
        count++;
    }
    matrix->size = count;
    matrix->nb_rows = sqrt(matrix->size);
    matrix->nb_cols = sqrt(matrix->size);
    fclose(fp);
}

int getElemFromMatrix(int i, int j, struct tablo *matrix) {
    if (i > matrix->nb_rows) {
        printf("i is higher than the number of the matrix lines\n");
        exit(1);
    } else if (j > matrix->nb_cols) {
        printf("j is higher than the number of the matrix columns\n");
        exit(1);
    } else if (j < 1) {
        printf("j is less than 1, it's not possible to get\n");
        exit(1);
    } else if (i < 1) {
        printf("i is less than 1, it's not possible to get\n");
        exit(1);
    } else {
        if (i == 1 && j == 1) {
            return (matrix->tab[0]);
        }
        if (i == 1 && j != 1) {
            return (matrix->tab[j - 1]);
        } else {
            return (matrix->tab[(i - 1) * matrix->nb_cols + j - 1]);
        }
    }
}

void setElemToMatrix(int i, int j, struct tablo *matrix, int value) {
    if (i > matrix->nb_rows) {
        printf("i is higher than the number of the matrix lines\n");
        exit(1);
    } else if (j > matrix->nb_cols) {
        printf("j is higher than the number of the matrix columns\n");
        exit(1);
    } else {
        if (i == 1 && j == 1) {
            matrix->tab[0] = value;
        }
        if (i == 1 && j != 1) {
            matrix->tab[j - 1] = value;
        } else {
            matrix->tab[(i - 1) * matrix->nb_cols + j - 1] = value;
        }
    }
}

int getMatrixSizeToAllocate(struct tablo *matrixA, struct tablo *matrixB) {
    return matrixA->nb_rows * matrixB->nb_cols;
}

struct tablo *getCol(int j, struct tablo *matrixSrc, struct tablo *result) {
    for (int i = 0; i < matrixSrc->nb_rows; i++) {
        // i + 1 parce qu'on peut pas get un élément 0. ça commence par 1
        result->tab[i] = getElemFromMatrix(i + 1, j, matrixSrc);
    }
    return result;
}

struct tablo *getRow(int i, struct tablo *matrixSrc, struct tablo *result) {
    for (int j = 0; j < matrixSrc->nb_cols; j++) {
        // i + 1 parce qu'on peut pas get un élément 0. ça commence par 1
        result->tab[j] = getElemFromMatrix(i, j + 1, matrixSrc);
    }
    return result;
}

int sumMultipyTwoVectors(struct tablo *matrixRow, struct tablo *matrixCol) {
    int res = 0;
    for (int i = 0; i < matrixRow->size; i++) {
        res = res + matrixRow->tab[i] * matrixCol->tab[i];
    }
    return res;
}

// depreciated
void multiply1(struct tablo *matrixA, struct tablo *matrixB, struct tablo *matrixResult) {
    if (matrixA->nb_cols != matrixB->nb_rows) {
        printf("nb de colonnes de la première matrice est different de nombre de ligne de deuxieme matrice");
        exit(1);
    }

    int i, j, k;
    for (i = 0; i < matrixA->nb_rows; i++) {
        for (j = 0; j < matrixB->nb_cols; j++) {
            setElemToMatrix(i + 1, j + 1, matrixResult, 0);
            for (k = 0; k < matrixA->nb_rows; k++)
                setElemToMatrix(i + 1, j + 1, matrixResult, getElemFromMatrix(i + 1, j + 1, matrixResult) +
                                                            getElemFromMatrix(i + 1, k + 1, matrixA) *
                                                            getElemFromMatrix(k + 1, j + 1, matrixB));
        }
    }
}

// depreciated
void multiply2(struct tablo *matrixA, struct tablo *matrixB, struct tablo *matrixResult) {
    if (matrixA->nb_cols != matrixB->nb_rows) {
        printf("nb de colonnes de la première matrice est different de nombre de ligne de deuxieme matrice");
        exit(1);
    }
    for (int iA = 1; iA < matrixA->nb_rows + 1; iA++) {

        for (int jB = 1; jB < matrixB->nb_cols + 1; jB++) {

            int sum = 0;

            for (int iter = 1; iter < matrixA->nb_cols + 1; iter++) {
                // printf("%d * %d\n",getElemFromMatrix(iA, iter, matrixA), getElemFromMatrix(iter, jB, matrixB));
                sum = sum + getElemFromMatrix(iA, iter, matrixA), getElemFromMatrix(iter, jB, matrixB);
            }
            setElemToMatrix(iA, jB, matrixResult, sum);

        }
    }
}

// the one that i use
void multiply3(struct tablo *matrixA, struct tablo *matrixB, struct tablo *matrixResult) {

    struct tablo *matrixARow = allocateTablo(matrixA->nb_cols);
    struct tablo *matrixBCol = allocateTablo(matrixB->nb_rows);
    int value = 0;
    for (int j = 0; j < matrixB->nb_cols; j++) {
        for (int i = 0; i < matrixA->nb_rows; i++) {

            getRow(i + 1, matrixA, matrixARow);
            getCol(j + 1, matrixB, matrixBCol);

            // ok
            value = sumMultipyTwoVectors(matrixARow, matrixBCol);
            setElemToMatrix(i + 1, j + 1, matrixResult, value);
        }
    }

    free(matrixARow);
    free(matrixBCol);
}

// depreciated
void multiply_omp(struct tablo *matrixA, struct tablo *matrixB, struct tablo *matrixResult) {
    {
        int i, j, k;
        for (i = 0; i < matrixA->nb_rows; i++) {
            for (j = 0; j < matrixB->nb_cols; j++) {
                setElemToMatrix(i + 1, j + 1, matrixResult, 0);
                for (k = 0; k < matrixA->nb_rows; k++)
                    setElemToMatrix(i + 1, j + 1, matrixResult, getElemFromMatrix(i + 1, j + 1, matrixResult) +
                                                                getElemFromMatrix(i + 1, k + 1, matrixA) *
                                                                getElemFromMatrix(k + 1, j + 1, matrixB));
            }
        }
    };

}

void getSubMatrix(struct tablo *matrix, struct tablo *submatrix, int startingrow, int endingrow, int startingcol,
                  int endingcol) {
    // printf("getting the submatrix of startrow : %d, endrow : %d, startcol : %d, endcol : %d \n", startingrow, endingrow,
    // startingcol, endingcol);
    if (endingrow < startingrow || endingrow > matrix->nb_rows) {
        printf("GetSubMatrix : error int the ending row input");
        exit(1);
    }
    if (endingcol < startingcol || endingcol > matrix->nb_cols) {
        printf("GetSubMatrix : error in the ending col input");
        exit(1);
    }
    if (startingrow < 1 || startingrow > matrix->nb_rows || startingrow > endingrow) {
        printf("GetSubMatrix : error in the starting row input");
        exit(1);
    }
    if (startingcol < 1 || startingcol > matrix->nb_cols || startingcol > endingcol) {
        printf("GetSubMatrix : error in the starting col input");
        exit(1);
    }

    int k = 0;
    for (int i = startingrow; i <= endingrow; i++) {
        for (int j = startingcol; j <= endingcol; j++) {
            submatrix->tab[k] = getElemFromMatrix(i, j, matrix);
            k++;
        }
    }
    submatrix->nb_rows = endingrow - startingrow + 1;
    submatrix->nb_cols = endingcol - startingcol + 1;
}

int getSizeSubMatrix(struct tablo *matrix, int totalProcess) {
    int nbrows = totalProcess / matrix->nb_rows;
    if ((totalProcess % matrix->nb_rows) != 0) {
        nbrows = nbrows + 1;
    }
    return nbrows * matrix->nb_cols;
}

void getRowsSubMatrix(struct tablo *matrix, struct tablo *submatrix, int nbOfProcess, int totalProcess) {
    printf("getRowsSubMatrix:\n");
    int step = totalProcess / sqrt(submatrix->size);
    if (step < 1) {
        step = 1;
    }
    int startingRow = step * nbOfProcess + 1;
    int endingRow = startingRow + step - 1;
    int startingCol = 1;
    int endingCol = matrix->nb_cols;
    getSubMatrix(matrix, submatrix, startingRow, endingRow, startingCol, endingCol);
}

void getColsSubMatrix(struct tablo *matrix, struct tablo *submatrix, int nbOfProcess, int totalProcess) {
    printf("getColsSubMatrix:\n");
    int step = totalProcess / sqrt(submatrix->size);
    if (step < 1) {
        step = 1;
    }
    int startingRow = 1;
    int endingRow = matrix->nb_rows;
    int startingCol = step * nbOfProcess + 1;
    int endingCol = startingCol + step - 1;
    getSubMatrix(matrix, submatrix, startingRow, endingRow, startingCol, endingCol);
}

int getSizeSubmatrixDividedByRows(struct tablo *matrix, int nbTotalProcess) {
    // printf("matrix->nb_rows : %d, nbTotalProcess : %d", matrix->nb_rows, nbTotalProcess);
    int nbRowsByProcess = matrix->nb_rows / nbTotalProcess;
    // printf("nb rows by process : %d \n", nbRowsByProcess);
    return nbRowsByProcess * matrix->nb_cols;
}

int getSizeSubmatrixDividedByCols(struct tablo *matrix, int nbTotalProcess) {
    int nbColsByProcess = matrix->nb_cols / nbTotalProcess;
    //printf("nb cols by process : %d \n", nbColsByProcess);
    return nbColsByProcess * matrix->nb_rows;
}

void getSubmatrixDividedByRows(struct tablo *matrix, struct tablo *subMatrix, int numberOfProcess, int nbTotalProcess) {
    // printf("[getSubmatrixDividedByRows]---------------------------------------------------------------------\n");
    int nbRowsBySlice = matrix->nb_rows / nbTotalProcess;

    int startingrow = nbRowsBySlice * numberOfProcess + 1;
    int endingrow = nbRowsBySlice * (numberOfProcess + 1);
    int startingcol = 1;
    int endingcol = matrix->nb_cols;
    // printf("nbRowsBySlice %d, startingrow %d, endingrow %d, startingcol %d, endingcol %d \n", nbRowsBySlice, startingrow, endingrow, startingcol, endingcol);
    getSubMatrix(matrix, subMatrix, startingrow, endingrow, startingcol, endingcol);
    // prettyPrintMatrix(subMatrix);
    // printf("-----------------------------------------------------------------------------------------------\n");
}

struct tablo * getSubmatrixDividedByCols(struct tablo *matrix, struct tablo *subMatrix, int numberOfProcess, int nbTotalProcess) {
    // printf("[getSubmatrixDividedByRows]\n");
    // printf("matrix nb of rows : %d, col : %d..\n", matrix->nb_rows, subMatrix->nb_cols);
    int nbColsBySlice = matrix->nb_rows / nbTotalProcess;
    int startingrow = 1;
    int endingrow = matrix->nb_rows;
    int startingcol = nbColsBySlice * numberOfProcess + 1;
    int endingcol = nbColsBySlice * (numberOfProcess + 1);
    getSubMatrix(matrix, subMatrix, startingrow, endingrow, startingcol, endingcol);
    return subMatrix;
}

void send_slice(int * slice_to_send, int its_size, int i) {
    MPI_Send(slice_to_send, its_size, MPI_INT, i, 0, MPI_COMM_WORLD);
}

int * recieve_slice(int sender_rank, int size_to_receive) {

    int * tablo = malloc(sizeof(int)*size_to_receive);

    MPI_Recv(tablo, size_to_receive, MPI_INT, sender_rank, 0, MPI_COMM_WORLD, 0);

    return tablo;
}

void fillMatrix(struct tablo *bigMatrix, struct tablo *littleMatrix, int offset) {
    if (bigMatrix->nb_rows != littleMatrix->nb_rows) {
        printf("the two matrix do not have the same rows number \n");
        exit(1);
    } else {
#pragma omp parallel for
        for (int i = 1; i <= littleMatrix->nb_rows; i++) {
            for (int j = 1; j <= littleMatrix->nb_cols; j++) {
                int value = getElemFromMatrix(i, j, littleMatrix);
                //printf("littlematrix(%d, %d) value : %d\n",i, j, value);
                // printf("bigmatrix(%d, %d) value : %d\n", i, j + offset, value);
                setElemToMatrix(i, j + offset, bigMatrix, value);
            }
        }
    }
}

void fillFinalMatrix(struct tablo *finalMatrix, struct tablo *bigMatrix, int offset) {
    if (finalMatrix->nb_cols != bigMatrix->nb_cols) {
        printf("the two matrix do not have the same cols number \n");
        exit(1);
    } else {
#pragma omp parallel for
        for (int i = 1; i <= bigMatrix->nb_rows; i++) {
            for (int j = 1; j <= bigMatrix->nb_cols; j++) {
                int value = getElemFromMatrix(i, j, bigMatrix);
                //printf("littlematrix(%d, %d) value : %d\n",i, j, value);
                // printf("bigmatrix(%d, %d) value : %d\n", i, j + offset, value);
                setElemToMatrix(i+offset, j, finalMatrix, value);
            }
        }
    }
}

void prettyPrintMatrix(struct tablo *Matrix) {
    for (int i = 1; i <= Matrix->nb_rows; i++) {
        for (int j = 1; j <= Matrix->nb_cols; j++) {
            printf("%d", getElemFromMatrix(i, j, Matrix));
            if(j <= Matrix->nb_cols - 1){
                printf(" ");
            }
        }
        if(i != Matrix->nb_rows){
            printf("\n");
        }
    }
    printf("\n");
}

void rotateMatrixToTheRight(struct tablo *matrix, int k){
    // temporary array of size M
    int temp[matrix->nb_cols];

    // within the size of matrix
    k = k % matrix->nb_cols;

    for (int i = 0; i < matrix->nb_rows; i++) {

        // copy first M-k elements to temporary array
        for (int t = 0; t < matrix->nb_cols - k; t++)
            temp[t] = getElemFromMatrix(i+1, t+1, matrix);

        // copy the elements from k to end to starting
        for (int j = matrix->nb_cols - k; j < matrix->nb_cols; j++)
            setElemToMatrix(i+1, j - matrix->nb_cols + k +1, matrix, getElemFromMatrix(i+1,j+1,matrix));

        // copy elements from temporary array to end
        for (int j = k; j < matrix->nb_cols; j++)
            setElemToMatrix(i+1, j+1, matrix, temp[j - k]);
    }
}

int main(int argc, char *argv[]) {

    int dim = getNbColsWhenRead(argv[1]);


    if (argc < 3) {
        printf("Error, arguments missing\n");
        return EXIT_FAILURE;
    } else {

        // MPI attributs
        int rank, numprocs;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Status status;
        int predicted_size;

        int step = 0;

        // printf("Initialize little matrix for multiplication result\n");
        struct tablo *littleMatrix = allocateTablo((dim / numprocs) *  (dim / numprocs));
        littleMatrix->nb_cols = dim / numprocs;
        littleMatrix->nb_rows = dim / numprocs;

        // printf("Initialize big matrix for partial result\n");
        struct tablo *bigMatrix = allocateTablo(dim * dim / numprocs);
        bigMatrix->nb_rows = dim / numprocs;
        bigMatrix->nb_cols = dim;


        if (rank == 0) {  // P0

            // printf("[ PROC %d ] I have read matrix size of %d x %d\n", rank, dim, dim);

            // printf("[ PROC %d ] Initialize final matrix for the result\n", rank);
            struct tablo *finalMatrix = allocateTablo(dim * dim);
            finalMatrix->nb_rows = dim;
            finalMatrix->nb_cols = dim;

            // printf("[ PROC %d ] Initialize matrix A, B from file input\n", rank);
            // allocate space for matrix A, matrix B and the correct size for matrix result
            struct tablo *matrixA = allocateTablo(dim*dim);
            struct tablo *matrixB = allocateTablo(dim*dim);
            gettingMatrixDataFromFile(argv[1], matrixA);
            gettingMatrixDataFromFile(argv[2], matrixB);

            // printf("[ PROC %d ] i have detected %d processors (workers)\n", rank,numprocs);


            // printf("[ PROC %d ] get A rows for P0 locally\n", rank);
            // get A rows for P0 locally
            struct tablo * my_rows = allocateTablo(matrixA->nb_rows * matrixA->nb_rows/numprocs);
            my_rows->nb_rows = matrixA->nb_rows / numprocs;
            my_rows->nb_cols = matrixA->nb_cols;
            getSubmatrixDividedByRows(matrixA, my_rows, 0, numprocs);


            // printf("[ PROC %d ] P%d has recieved the following row slice... \n", rank,rank);
            // prettyPrintMatrix(my_rows);


            // printf("[ PROC %d ] sending rows to every processor\n", rank);
            // step 1 : Send rows
            for ( int i = 1; i < numprocs ; i ++ ) {
                struct tablo *SubMatrixAByRows = allocateTablo(getSizeSubmatrixDividedByRows(matrixA, numprocs));
                getSubmatrixDividedByRows(matrixA, SubMatrixAByRows, i, numprocs);
                send_slice(SubMatrixAByRows->tab, SubMatrixAByRows->size, i);
                // printf("[ PROC %d ] \t sending the following row slice to the processor %d ... \n", rank, i);
                // prettyPrintMatrix(SubMatrixAByRows);
                free(SubMatrixAByRows);
            }

            for(int rotation = 0; rotation < numprocs; rotation++){

                // step 2 : Send columns
                for ( int i = 1 ; i < numprocs ; i ++ ) {
                    struct tablo *SubMatrixAByCols = allocateTablo(getSizeSubmatrixDividedByCols(matrixB, numprocs));
                    getSubmatrixDividedByCols(matrixB, SubMatrixAByCols, i, numprocs);
                    send_slice(SubMatrixAByCols->tab, SubMatrixAByCols->size, i);
                    free(SubMatrixAByCols);
                }

                // get B columns for P0
                struct tablo * my_cols = allocateTablo(getSizeSubmatrixDividedByRows(matrixB, numprocs));
                getSubmatrixDividedByCols(matrixB, my_cols, 0, numprocs);
                my_cols->nb_rows = dim;
                my_cols->nb_cols = dim / numprocs;

                // Processing multiplication.
                struct tablo * my_product = allocateTablo(getMatrixSizeToAllocate(my_rows, my_cols));
                multiply3(my_rows,my_cols,my_product);
                fillMatrix(bigMatrix, my_product, ((rank+step)%numprocs * dim/numprocs) % dim);
                step++;
                rotateMatrixToTheRight(matrixB, (dim/numprocs * numprocs) - dim/numprocs);
                free(my_cols);

            }

            free(my_rows);


            fillFinalMatrix(finalMatrix, bigMatrix, rank * (dim / numprocs));

            //  recieve bigmatrix
            for ( int i = 1 ; i < numprocs ; i++ ) {
                struct tablo *recievedMatrix = allocateTablo(dim/numprocs * dim);
                recievedMatrix->tab = recieve_slice(i, recievedMatrix->size);
                recievedMatrix->nb_rows = dim / numprocs;
                recievedMatrix->nb_cols = dim;
                fillFinalMatrix(finalMatrix, recievedMatrix, i * dim/numprocs);
            }

            prettyPrintMatrix(finalMatrix);

            free(finalMatrix);
            free(matrixA);
            free(matrixB);
            free(my_rows);

        } else { // others processors
            MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_INT, &predicted_size);

            // Receive A rows
            struct tablo * my_rows = allocateTablo(predicted_size);
            my_rows->tab = recieve_slice(0, my_rows->size);
            my_rows->nb_rows = dim / numprocs;
            my_rows->nb_cols = dim;

            for(int rotation = 0; rotation < numprocs; rotation++) {
                // Receive B columns
                struct tablo *my_cols = allocateTablo(predicted_size);
                my_cols->tab = recieve_slice(0, my_cols->size);
                my_cols->nb_rows = dim;
                my_cols->nb_cols = dim / numprocs;

                // Processing multiplication.
                struct tablo *my_product = allocateTablo(getMatrixSizeToAllocate(my_rows, my_cols));
                multiply3(my_rows, my_cols, my_product);
                fillMatrix(bigMatrix, my_product, ((rank+step)%numprocs * dim/numprocs) % dim);
                step++;
                free(my_product);
                free(my_cols);
            }
            free(my_rows);
            send_slice(bigMatrix->tab, bigMatrix->size, 0);
        }
        free(littleMatrix);
        free(bigMatrix);
        MPI_Finalize();
    }
    return EXIT_SUCCESS;
}


