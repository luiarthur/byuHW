#include <stdio.h>      // standard input/output: printf

int main(int argc, char* argv[]){
  printf("%s","Hello");
  
  if (argc > 1){
    printf("%s",",");
    for (int i=1; i<argc; i++){
      printf("%s%s"," ",argv[i]);
    }
  }
  
  printf("%s\n",".");
}
