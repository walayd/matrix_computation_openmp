#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

struct tablo {
  int * tab;
  int size;
  int nb_cols;
  int nb_rows;
};

struct tablo * allocateTablo(int size) {
  struct tablo * tmp = malloc(sizeof(struct tablo));
  tmp->size = size;
  tmp->nb_cols = (int)sqrt(size);
  tmp->nb_rows = (int)sqrt(size);
  tmp->tab = malloc(size*sizeof(int));
  return tmp;
}

struct tablo * reallocateTabloTwice(struct tablo * tmp) {
  int oldsize = tmp->size;
  tmp->size = oldsize * 2; 
  tmp->tab = realloc(tmp->tab, tmp->size * sizeof(int));

  for(int i = oldsize; i<tmp->size; i++){
  	tmp->tab[i] = 0;
  }
  return tmp;
}

void printArray(struct tablo * tmp) {
  // printf("---- Array of size %i ---- \n", tmp->size);
  int size = tmp->size;
  int i;
  for (i = 0; i < size; ++i) {
    printf("%i ", tmp->tab[i]);
  }
  printf("\n");
}

int calculateMatrixLength(char filename[]){
	printf("Starting allocating matrix \n");
	int len = 0;
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
	while((c = fgetc(fp)) != EOF) {
		if (c == '\n') {
			break;
		}
		if (c != ' ' && c != '-'){
			len++;
		}
	}
 
	printf("Length of first line is: %d\n", len);

	// seek to beginning of file
	fseek(fp, 0, SEEK_SET);
	fclose(fp);
	// allocate memory for size of first line (len)
	return len*len;
}

void gettingMatrixDataFromFile(char filename[], struct tablo * matrix){
	FILE *fp;
	fp = fopen(filename, "r");
	int i=0;
	char c;
	int isNeg = 0; // false if 0
	while((c = fgetc(fp)) != EOF) {
	    if (c == '-'){
	        isNeg = 1;
	    }
		if (c != '\n' && c != ' ' && c != '-') {
			matrix->tab[i] = (int)(c - '0');
			// printf("i : %d, c : %d\n", i, (int)(c - '0'));
			if(isNeg == 1){
                matrix->tab[i] = matrix->tab[i] * -1;
                isNeg = 0;
			}
			i++;
		}
	}
}

int getElemFromMatrix(int i, int j, struct tablo * matrix){
    if(i > matrix->nb_rows){
        printf("i is higher than the number of the matrix lines\n");
        exit(1);
    }
    else if(j > matrix->nb_cols){
        printf("j is higher than the number of the matrix columns\n");
        exit(1);
    }
    else if(j < 1){
        printf("j is less than 1, it's not possible to get\n");
        exit(1);
    }
    else if(i < 1){
        printf("i is less than 1, it's not possible to get\n");
        exit(1);
    }
    else{
        return (matrix->tab[matrix->nb_rows * (i-1) + (j-1)]);
    }
}

int setElemToMatrix(int i, int j, struct tablo * matrix, int value){
    if(i > matrix->nb_rows){
        printf("i is higher than the number of the matrix lines\n");
        exit(1);
    }
    else if(j > matrix->nb_cols){
        printf("j is higher than the number of the matrix columns\n");
        exit(1);
    }
    else{
        return (matrix->tab[matrix->nb_rows * (i-1) + (j-1)] = value);
    }
}

int getMatrixSizeToAllocate(struct tablo * matrixA, struct tablo * matrixB){
    return matrixA->nb_rows * matrixB->nb_cols;
}

struct tablo * getCol(int j, struct tablo * matrixSrc, struct tablo * result){
    for(int i = 0; i < matrixSrc->nb_rows; i++) {
        // i + 1 parce qu'on peut pas get un élément 0. ça commence par 1
        result->tab[i] = getElemFromMatrix(i+1, j, matrixSrc);
    }
}

struct tablo * getRow(int i, struct tablo * matrixSrc, struct tablo * result){
    for(int j = 0; j < matrixSrc->nb_cols; j++) {
        // i + 1 parce qu'on peut pas get un élément 0. ça commence par 1
        result->tab[j] = getElemFromMatrix(i, j+1, matrixSrc);
    }
}

int sumMultipyTwoVectors(struct tablo * matrixRow, struct tablo * matrixCol){
    int res = 0;
    for(int i = 0; i < matrixRow->size; i++){
        res = res + matrixRow->tab[i] * matrixCol->tab[i];
    }
    return res;
}

// ligne fois column
// 1ere ligne, 1ere column pour (1,1) (row, col)
// 2eme ligne, 1ere column pour (2,1) (row, col)


void multiply1(struct tablo * matrixA, struct tablo * matrixB, struct tablo * matrixResult)
{
    if (matrixA->nb_cols != matrixB->nb_rows){
        printf("nb de colonnes de la première matrice est different de nombre de ligne de deuxieme matrice");
        exit(1);
    }

    int i, j, k;
    for (i = 0; i < matrixA->nb_rows; i++)
    {
        for (j = 0; j < matrixB->nb_cols; j++)
        {
            setElemToMatrix(i+1, j+1, matrixResult, 0);
            for (k = 0; k < matrixA->nb_rows; k++)
                setElemToMatrix(i+1, j+1, matrixResult, getElemFromMatrix(i+1, j+1, matrixResult) + getElemFromMatrix(i+1, k+1, matrixA) * getElemFromMatrix(k+1,j+1, matrixB));
        }
    }
}

void multiply2(struct tablo * matrixA, struct tablo * matrixB, struct tablo * matrixResult)
{
    if (matrixA->nb_cols != matrixB->nb_rows){
        printf("nb de colonnes de la première matrice est different de nombre de ligne de deuxieme matrice");
        exit(1);
    }
    for(int iA = 1; iA < matrixA->nb_rows+1; iA++){
        int sum = 0;
        for(int jB = 1; jB < matrixB->nb_cols+1; jB++){
            for(int k = 1; k < matrixB->nb_rows+1; k++){
                sum = sum + getElemFromMatrix(iA, k, matrixA) * getElemFromMatrix(k, jB, matrixB);
            }
            setElemToMatrix(iA, jB, matrixResult, sum);
            sum=0;
        }
    }
}

void multiply3(struct tablo * matrixA, struct tablo * matrixB, struct tablo * matrixResult){

    struct tablo * matrixARow = allocateTablo(matrixA->nb_cols);
    struct tablo * matrixBCol = allocateTablo(matrixB->nb_rows);
    int value = 0;
    for(int j = 0; j < matrixB->nb_cols; j++){
        for(int i = 0; i < matrixA->nb_rows; i++){

            getRow(i+1, matrixA, matrixARow);
            getCol(j+1, matrixB, matrixBCol);

            // ok
            value = sumMultipyTwoVectors(matrixARow, matrixBCol);
            setElemToMatrix(i+1,j+1, matrixResult, value);
        }
    }
}


void multiply_omp(struct tablo * matrixA, struct tablo * matrixB, struct tablo * matrixResult)
{
    #pragma omp parallel
    {
        int i, j, k;
        #pragma omp for
        for (i = 0; i < matrixA->nb_rows; i++)
        {
            for (j = 0; j < matrixB->nb_cols; j++)
            {
                setElemToMatrix(i+1, j+1, matrixResult, 0);
                for (k = 0; k < matrixA->nb_rows; k++)
                    setElemToMatrix(i+1, j+1, matrixResult, getElemFromMatrix(i+1, j+1, matrixResult) + getElemFromMatrix(i+1, k+1, matrixA) * getElemFromMatrix(k+1,j+1, matrixB));
            }
        }
    };

}

void getSubMatrix(struct tablo * matrix, struct tablo * submatrix, int startingrow, int endingrow, int startingcol, int endingcol){
    printf("getting the submatrix of startrow : %d, endrow : %d, startcol : %d, endcol : %d \n", startingrow, endingrow, startingcol, endingcol);
    if(endingrow < startingrow || endingrow > matrix->nb_rows){
        printf("GetSubMatrix : error int the ending row input");
        exit(1);
    }
    if(endingcol < startingcol || endingcol > matrix->nb_cols){
        printf("GetSubMatrix : error in the ending col input");
        exit(1);
    }
    if(startingrow < 1 || startingrow > matrix->nb_rows || startingrow > endingrow){
        printf("GetSubMatrix : error in the starting row input");
        exit(1);
    }
    if(startingcol < 1 || startingcol > matrix->nb_cols || startingcol > endingcol){
        printf("GetSubMatrix : error in the starting col input");
        exit(1);
    }

    int k = 0;
    for(int i = startingrow; i <= endingrow; i++){
        for(int j = startingcol; j<= endingcol; j++){
            submatrix->tab[k] = getElemFromMatrix(i, j, matrix);
            k++;
        }
    }
    submatrix->nb_rows = endingrow - startingrow + 1;
    submatrix->nb_cols = endingcol - startingcol + 1;
}

int getSizeSubMatrix(struct tablo * matrix, int totalProcess){
    int nbrows = totalProcess / matrix->nb_rows;
    if((totalProcess % matrix->nb_rows) != 0){
        nbrows = nbrows + 1;
    }
    return nbrows * matrix->nb_cols;
}

void getRowsSubMatrix(struct tablo * matrix, struct tablo * submatrix, int nbOfProcess, int totalProcess){
    printf("getRowsSubMatrix:\n");
    int step = totalProcess / sqrt(submatrix->size);
    if (step < 1){
        step = 1;
    }
    int startingRow = step * nbOfProcess + 1;
    int endingRow = startingRow + step - 1;
    int startingCol = 1;
    int endingCol = matrix->nb_cols;
    getSubMatrix(matrix, submatrix, startingRow, endingRow, startingCol, endingCol);
}

void getColsSubMatrix(struct tablo * matrix, struct tablo * submatrix, int nbOfProcess, int totalProcess){
    printf("getColsSubMatrix:\n");
    int step = totalProcess / sqrt(submatrix->size);
    if (step < 1){
        step = 1;
    }
    int startingRow = 1;
    int endingRow = matrix->nb_rows;
    int startingCol = step * nbOfProcess + 1;
    int endingCol = startingCol + step - 1;
    getSubMatrix(matrix, submatrix, startingRow, endingRow, startingCol, endingCol);
}


int main(int argc, char* argv[]) {

    double dtime;
	if (argc < 3) {
		printf("Error, arguments missing\n");	
		return EXIT_FAILURE;
	} else {
		printf("first arg : %s\nsecond arg : %s\n", argv[1], argv[2]);

		// calculate matrix length to know how much size to allocate
		int matrix_length = calculateMatrixLength(argv[1]);

		// allocate space for matrix A, matrix B and the correct size for matrix result
		struct tablo * matrixA = allocateTablo(matrix_length);
		struct tablo * matrixB = allocateTablo(matrix_length);
		printf("size to allocate for result matrix %d\n", getMatrixSizeToAllocate(matrixA, matrixB));
		struct tablo * matrixResult = allocateTablo(getMatrixSizeToAllocate(matrixA, matrixB));

		// getting matrix data (of A and B) from file
		gettingMatrixDataFromFile(argv[1], matrixA);
		gettingMatrixDataFromFile(argv[2], matrixB);

		printArray(matrixA);
		printArray(matrixB);

		/*
        printf("testing multiply 1 without openmp\n");
        dtime = omp_get_wtime();
        multiply1(matrixA, matrixB, matrixResult);
        dtime = omp_get_wtime() - dtime;
        printf("%f\n", dtime);
        printArray(matrixResult);

        // multiply 2 seems to be the best one.
        printf("testing multiply 2 without openmp\n");
        dtime = omp_get_wtime();
        multiply2(matrixA, matrixB, matrixResult);
        dtime = omp_get_wtime() - dtime;
        printf("%f\n", dtime);
        printArray(matrixResult);

        printf("testing multiply 3 without openmp\n");
        dtime = omp_get_wtime();
        multiply3(matrixA, matrixB, matrixResult);
        dtime = omp_get_wtime() - dtime;
        printf("%f\n", dtime);
        printArray(matrixResult);

        printf("testing multiply with openmp\n");
        dtime = omp_get_wtime();
        multiply_omp(matrixA, matrixB, matrixResult);
        dtime = omp_get_wtime() - dtime;
        printf("%f\n", dtime);
        printArray(matrixResult);
        */


		
        /*// testing submatrix (of 4 elements)
        struct tablo * subMatrix = allocateTablo(4);
        getSubMatrix(matrixA, subMatrix, 1, 2, 2, 3);
        printArray(subMatrix);
        printf("nb cols : %d\n", subMatrix->nb_cols);
        printf("nb rows : %d\n", subMatrix->nb_rows);


        struct tablo * matrixACol = allocateTablo(sqrt(matrix_length));
        struct tablo * matrixARow = allocateTablo(sqrt(matrix_length));

        int sizesubmatrix = getSizeSubMatrix(matrixA, 3);



        struct tablo * subMatrixA = allocateTablo(sizesubmatrix);
        struct tablo * subMatrixAA = allocateTablo(sizesubmatrix);
        struct tablo * multiplysubmatrix = allocateTablo(1);
        getRowsSubMatrix(matrixA, subMatrixA, 0, 1);
        getColsSubMatrix(matrixA, subMatrixAA, 0, 1);
        printArray(subMatrixA);
        printArray(subMatrixAA);
        multiply(subMatrixA, subMatrixAA , multiplysubmatrix);
        printArray(multiplysubmatrix);*/



	}

	// scatter pour A, scatter pour B
	// matrice vecteur donc il faut un autre for

	return EXIT_SUCCESS;
	
}










