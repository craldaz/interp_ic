#include "interp.h"

void Interp::print_xyz(){

 printf(" %i \n",natoms);
 printf("\n");
 for (int n=0;n<2;n++)
	{cout << "molecule" << n << endl; 
 	for (int i=0;i<natoms;i++) 
 		{
     cout << "  " << anames[i];
     printf(" %f %f %f \n",coords[n][3*i+0],coords[n][3*i+1],coords[n][3*i+2]);
 		}
	}
// printf("\n");

}

