#include "digitalbrain.h"

int main(int argc, char **argv)
{
   // printf() displays the string inside quotation
   //printf("Hello, World and Digital Brain!! This is example 1!\n\n");

   FILE *fp;

   if (argc >= 2)
       fp = fopen(argv[1], "r");
   else
        fprintf(stderr, "Input Error! Use: ./ex1 keyword.k \n");

   //read inputfile and initalize
   ReadInputFile();

   return 0;
}
