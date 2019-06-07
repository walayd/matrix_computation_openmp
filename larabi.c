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
    int nbCols = 0;
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
            break;
        }
        if ((c != ' ') && (c != '-')) {
            nbCols++;
        }
    }

    // seek to beginning of file
    fseek(fp, 0, SEEK_SET);
    fclose(fp);
    // allocate memory for size of first line (len)
    return nbCols;
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
}

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
    // printf("[getSubmatrixDividedByRows]\n");
    int nbRowsBySlice = matrix->nb_rows / nbTotalProcess;
    int startingrow = nbRowsBySlice * numberOfProcess + 1;
    int endingrow = nbRowsBySlice * (numberOfProcess + 1);
    int startingcol = 1;
    int endingcol = matrix->nb_cols;
    // printf("nbRowsBySlice %d, startingrow %d, endingrow %d, startingcol %d, endingcol %d \n", nbRowsBySlice, startingrow, endingrow, startingcol, endingcol);
    getSubMatrix(matrix, subMatrix, startingrow, endingrow, startingcol, endingcol);
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
        for (int i = 1; i <= littleMatrix->nb_rows; i++) {
            for (int j = 1; j <= littleMatrix->nb_cols; j++) {
                int value = getElemFromMatrix(i, j, littleMatrix);
                //printf("littlematrix(%d, %d) value : %d\n",i, j, value);
                printf("bigmatrix(%d, %d) value : %d\n", i, j + offset, value);

                setElemToMatrix(i, j + offset, bigMatrix, value);
            }
        }
    }
}

void prettyPrintMatrix(struct tablo *Matrix) {
    printf("------------------- \n");
    for (int i = 1; i <= Matrix->nb_rows; i++) {
        for (int j = 1; j <= Matrix->nb_cols; j++) {
            printf(" %d", getElemFromMatrix(i, j, Matrix));
            // printf("(%d,%d) : %d\n", i, j, getElemFromMatrix(i, j, Matrix));
        }
        printf("\n");
    }
    printArray(Matrix);
    printf("Matrix length : %d. Nb Cols : %d. Nb Rows : %d \n", Matrix->size, Matrix->nb_cols, Matrix->nb_rows);
    printf("------------------- \n");
}

int main(int argc, char *argv[]) {

    int dim = getNbColsWhenRead(argv[1]);

    // printf("dim = %d \n", dim);

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

        if (rank == 0) {  // P0
            // printf("first arg : %s\nsecond arg : %s\n", argv[1], argv[2]);

            // printf(" matrix length : %d \n" , matrix_length2);

            // allocate space for matrix A, matrix B and the correct size for matrix result
            struct tablo *matrixA = allocateTablo(dim*dim);
            struct tablo *matrixB = allocateTablo(dim*dim);

            printf("A : size %d, rows %d, cols %d \n", matrixA->size, matrixA->nb_rows, matrixA->nb_cols);
            printf("B : size %d, rows %d, cols %d \n", matrixB->size, matrixB->nb_rows, matrixB->nb_cols);
            // getting matrix data (of A and B) from file
            gettingMatrixDataFromFile(argv[1], matrixA);
            gettingMatrixDataFromFile(argv[2], matrixB);

            // prettyPrintMatrix(matrixA);
            // prettyPrintMatrix(matrixB);

            /*

            printf("TESTING MULTIPLY MATRICES CARRES \n");

            printf("size to allocate for result matrix  : %d\n", getMatrixSizeToAllocate(matrixA, matrixB));
            struct tablo *matrixResult = allocateTablo(getMatrixSizeToAllocate(matrixA, matrixB));
            printf("Matrix result size: %d\n", matrixResult->size);
            multiply3(matrixA, matrixB, matrixResult);
            prettyPrintMatrix(matrixResult);

            int matrix_length3 = 6;
            int matrix_length4 = 6;

            struct tablo *matrix3 = allocateTablo(matrix_length3);
            struct tablo *matrix4 = allocateTablo(matrix_length4);
            */


            /*
            printf("TESTING MULTIPLY MATRICES NON CARRES : ok \n");

            matrix3->tab = (int[6]){1,-2,3,4,5,6};
            matrix3->nb_rows = 2;
            matrix3->nb_cols = 3;

            matrix4->tab = (int[6]){1,-2,4,5,7,8};
            matrix4->nb_rows = 3;
            matrix4->nb_cols = 2;

            prettyPrintMatrix(matrix3);
            prettyPrintMatrix(matrix4);

            printf("size to allocate for result matrix  : %d\n", getMatrixSizeToAllocate(matrix3, matrix4));
            struct tablo *matrixResult2 = allocateTablo(getMatrixSizeToAllocate(matrix3, matrix4));
            printf("Matrix result size: %d\n", matrixResult2->size);
            multiply3(matrix3, matrix4, matrixResult2);
            prettyPrintMatrix(matrixResult2);
            */

            /*
            printf("TESTING GETTING SUBMATRIX DEVIDED BY ROW/COL : ok \n");

            printf("size submatrix divided by rows : %d \n", getSizeSubmatrixDividedByRows(matrixA, 3));
            printf("size submatrix divided by cols : %d \n", getSizeSubmatrixDividedByCols(matrixB, 3));

            struct tablo *SubMatrixAByRows = allocateTablo(getSizeSubmatrixDividedByRows(matrixA, 3));
            struct tablo *SubMatrixAByCols = allocateTablo(getSizeSubmatrixDividedByCols(matrixA, 3));

            getSubmatrixDividedByRows(matrixA, SubMatrixAByRows, 0, 3);
            getSubmatrixDividedByCols(matrixA, SubMatrixAByCols, 0, 3);

            prettyPrintMatrix(SubMatrixAByRows);
            prettyPrintMatrix(SubMatrixAByCols);

            */

            /*
            printf("TESTING FILLING SLICE OF BIG MATRIX BY SMALL MATRIX : ok \n");
            struct tablo * bigMatrix = allocateTablo(18);
            bigMatrix->nb_cols = 6;
            bigMatrix->nb_rows = 3;
            printf("for bigMatrix, it has %d rows and %d cols \n", bigMatrix->nb_rows, bigMatrix->nb_cols);
            printf("for matrix A, it has %d rows and %d cols \n", matrixA->nb_rows, matrixA->nb_cols);
            fillMatrix(bigMatrix, matrixA, 3);
            printf("---------\n");
            prettyPrintMatrix(bigMatrix);
            printf("for bigMatrix, it has %d rows and %d cols \n", bigMatrix->nb_rows, bigMatrix->nb_cols);
            */



            // step 1 : Send rows
            for ( int i = 0 ; i < numprocs ; i ++ ) {
                // printf("nb of process : %d \n", numprocs);
                // printf("size submatrix divided by rows : %d \n", getSizeSubmatrixDividedByRows(matrixA, numprocs));
                struct tablo *SubMatrixAByRows = allocateTablo(getSizeSubmatrixDividedByRows(matrixA, numprocs));
                getSubmatrixDividedByRows(matrixA, SubMatrixAByRows, i, numprocs);
                // printf("for proc : %d, here's the row slice: \n", i);
                // prettyPrintMatrix(SubMatrixAByRows);
                send_slice(SubMatrixAByRows->tab, SubMatrixAByRows->size, i);
            }

            // Receive A rows
            struct tablo * my_rows = allocateTablo(getSizeSubmatrixDividedByRows(matrixA, numprocs));
            getSubmatrixDividedByRows(matrixA, my_rows, 0, numprocs);
            my_rows->nb_rows = dim / numprocs;
            my_rows->nb_cols = dim;
            printf("i'm the proc : %d, and here's my row slice: \n", rank);
            prettyPrintMatrix(my_rows);

            // step 2 : Send columns
            for ( int i = 0 ; i < numprocs ; i ++ ) {
                // printf("nb of process : %d \n", numprocs);
                // printf("size submatrix divided by rows : %d \n", getSizeSubmatrixDividedByRows(matrixA, numprocs));
                struct tablo *SubMatrixAByCols = allocateTablo(getSizeSubmatrixDividedByCols(matrixB, numprocs));
                getSubmatrixDividedByCols(matrixB, SubMatrixAByCols, i, numprocs);
                // printf("for proc : %d, here's the column slice: \n",i);
                // prettyPrintMatrix(SubMatrixAByCols);
                send_slice(SubMatrixAByCols->tab, SubMatrixAByCols->size, i);
            }

            // printf("///////////////////////////////////////////////\n");
            // printf("proc rank : %d\n", rank);
            // prettyPrintMatrix(my_rows);
            // prettyPrintMatrix(my_cols)
            // printf("///////////////////////////////////////////////\n");

        } else {  // others
            MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_INT, &predicted_size);

            // Receive A rows
            struct tablo * my_rows = allocateTablo(predicted_size);
            my_rows->tab = recieve_slice(0, my_rows->size);
            my_rows->nb_rows = dim / numprocs;
            my_rows->nb_cols = dim;
            printf("i'm the proc : %d, and here's my row slice: \n", rank);
            prettyPrintMatrix(my_rows);

            // Receive B columns
            // struct tablo * my_cols = allocateTablo(predicted_size);
            // my_cols->tab = recieve_slice(0, my_cols->size);
            // my_cols->nb_rows = dim;
            // my_cols->nb_cols = dim / numprocs;

            // printf("///////////////////////////////////////////////\n");
            // printf("proc rank : %d\n", rank);
            // printf("my rows : \n");
            // prettyPrintMatrix(my_rows);
            // printf("my columns : \n");
            // prettyPrintMatrix(my_cols);
            // printf("///////////////////////////////////////////////\n");

        }

        MPI_Finalize();

    }

    return EXIT_SUCCESS;

}










