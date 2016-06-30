#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "interp.h"

using namespace std; 

int main(int argc, char* argv[]){ 
	string xyzfile1;
	if (argc < 2){
		cout << "Must include xyzfile" <<endl;    
		return -1; 
	}
	xyzfile1=argv[1]; 
	
	//ICoord molecule1;
	//ICoord molecule2; 
	//molecule1.isOpt=1;
	//molecule1.init(xyzfile1);
	//molecule2.init(xyzfile2);
	
	Interp try1; 
	try1.init(xyzfile1);
	try1.calc_interp();
	
	//molecule1.union_ic(molecule1,molecule2); 
	//printf("\n actual IC's \n"); 
	//molecule1.print_ic(); 
  
	return 0;
}


