#include "digitalbrain.h"

int main(int argc, char **argv){

   if (argc >= 1){
      //printf("Will simulate on %s \n",argv[1]);
   }else{
      printf("Input Error! Use: ./ex1 keyword.k \n");
    }

   //read inputfile and initalize
   ReadInputFile(argv[1]);

   FreeArrays();
   return 0;
}
