#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
		if (c != ' '){
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
	while((c = fgetc(fp)) != EOF) {
		if (c != '\n' && c != ' ') {
			matrix->tab[i] = (int)(c - '0');
			// printf("i : %d, c : %d\n", i, (int)(c - '0'));
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

void multiplyMatrix(struct tablo * matrixA, struct tablo * matrixB, struct tablo * matrixResult){

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


int main(int argc, char* argv[]) {

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
		struct tablo * matrixResult = allocateTablo(matrix_length);

		// getting matrix data (of A and B) from file
		gettingMatrixDataFromFile(argv[1], matrixA);
		gettingMatrixDataFromFile(argv[2], matrixB);

        struct tablo * matrixACol = allocateTablo(sqrt(matrix_length));
        struct tablo * matrixARow = allocateTablo(sqrt(matrix_length));


		printArray(matrixA);
		printArray(matrixB);

        multiplyMatrix(matrixA, matrixB, matrixResult);
        printArray(matrixResult);
	}
	
	// scatter pour A, scatter pour B
	// matrice vecteur donc il faut un autre for

	return EXIT_SUCCESS;
	
}










