#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
   if (argc < 3){
     printf("Error, arguments missing\n");	
     return 1;
   }

   printf("first arg : %s\nsecond arg : %s\n", argv[1], argv[2]);
   return EXIT_SUCCESS;
}



