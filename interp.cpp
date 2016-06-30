#include "interp.h"


int Interp::init(string xyzfile){


//	printf(" xyzfile: %s \n",xyzfile);
  printf("\n");
  cout << " xyzfile: " << xyzfile << endl;
  structure_read(xyzfile);
 
  print_xyz();

	alloc_mem();
	printf(" done allocating memory\n");

//	int done = ic_create();

 //printf("\n\n");

 return 1;
}

void Interp::structure_read(string xyzfile)
{

  printf("Reading and initializing string coordinates \n");
  printf("  -Opening structure file \n");

  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    printf("\n Error opening xyz file: %s \n",xyzfile.c_str());
    exit(-1);
  }

  printf("  -reading file... \n");

  string line;
  bool success=true;
  success=getline(infile, line);
  if (success){
    int length=StringTools::cleanstring(line);
    natoms=atoi(line.c_str());
  }
  cout <<"  -The number of atoms is: " << natoms << endl;

  success=getline(infile, line);

  anumbers = new int[1+natoms];
  amasses = new double[1+natoms];
  anames = new string[1+natoms];

  cout <<"  -Reading the atomic names...";
  for (int i=0;i<natoms;i++){
    success=getline(infile, line);
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    anames[i]=tok_line[0];
    anumbers[i]=PTable::atom_number(anames[i]);
    amasses[i]=PTable::atom_mass(anumbers[i]);
  }

  infile.close();

  coords = new double*[1+2];


  for (int i=0;i<2;i++){
    coords[i] = new double[1+natoms*3];
  }

  cout <<"  -Reading coordinates...";
  printf("Opening xyz file \n");
  infile.open(xyzfile.c_str());
  fflush(stdout);
  //cout << "xyzfile opened" << endl;


  for (int i=0;i<2;i++)
  {
    success=getline(infile, line);
    success=getline(infile, line);
    for (int j=0;j<natoms;j++)
    {
      if (infile.eof())
      {
        printf("   end of xyz file reached early, exiting \n");
        exit(1);
      }
      success=getline(infile, line);
      int length=StringTools::cleanstring(line);
      vector<string> tok_line = StringTools::tokenize(line, " \t");
//      cout << " i: " << i << " string: " << line << endl;
      int n;
      if (i==0) n = 0;
      else if (i==1) n = 2-1;
      coords[n][3*j+0]=atof(tok_line[1].c_str());
      coords[n][3*j+1]=atof(tok_line[2].c_str());
      coords[n][3*j+2]=atof(tok_line[3].c_str());
    }
  }

  //cout << " done" << endl;
  infile.close();

  double zero = 0.;
  for (int i=0;i<3*natoms;i++)
    zero += coords[0][i]*coords[0][i];
  if (zero<0.0001) 
  {
    printf("\n ERROR: initial.xyz has NULL coordinates \n");
    exit(1);
  }

  cout << "Finished reading information from structure file" << endl;
}

void Interp::calc_interp()
{
	cout << "Starting" << endl; 
  int N3 = natoms*3;
  icoords = new ICoord[2+1];
  for (int i=0;i<2;i++)
    icoords[i].alloc(natoms);
  icoords[0].reset(natoms,anames,anumbers,coords[0]);
  icoords[2-1].reset(natoms,anames,anumbers,coords[2-1]);
	cout << "I'm here about to print icoords[0]" << endl; 
	icoords[0].print_xyz();
	icoords[1].print_xyz();

  ICoord ic1,ic2,ic3; 
  ic1.alloc(natoms);
  ic2.alloc(natoms);
  ic3.alloc(natoms);
  ic1.reset(natoms,anames,anumbers,coords[0]);
  ic3.reset(natoms,anames,anumbers,coords[0]);
  ic2.reset(natoms,anames,anumbers,coords[2-1]);
  ic1.ic_create();
  ic2.ic_create();
	//cout << "I'm here testing copy_ic" << endl; 
 	//ic3.copy_ic(ic2); 
	//cout << "I'm here about to print ic3" << endl; 
	//ic3.print_ic();	
	//ic1.print_xyz();
//	ic1.print_ic();
	//ic2.print_ic();

  allcoords = new double*[2];
  for (int i=0;i<2;i++)
    allcoords[i] = icoords[i].coords;

// create union_ic
 // newic.alloc(natoms);
 // intic.alloc(natoms);
 // int2ic.alloc(natoms);
 // newic.reset(natoms,anames,anumbers,icoords[0].coords);
  //intic.reset(natoms,anames,anumbers,icoords[2-1].coords);
  //int2ic.reset(natoms,anames,anumbers,icoords[2-1].coords);
	
	//newic.union_ic(ic1,ic2);  
  //newic.bmat_alloc();
  //newic.bmat_create();


  //printf("\n actual IC's \n");
  //newic.print_ic();


  //int size_ic = newic.nbonds + newic.nangles + newic.ntor;
  //int len_d = newic.nicd0;

	ic3.union_ic(ic1,ic2);
  printf("\n actual IC's \n");
  ic3.print_ic();
	ic3.bmat_alloc();
	ic3.bmatp_create();
	ic3.bmatp_to_U();
	ic3.bmat_create();
//	ic3.print_q();
  printf("\n");
  for (int n=0;n<2;n++)
    icoords[n].copy_ic(ic3);
  for (int n=0;n<2;n++)
    icoords[n].bmat_alloc();

#if 0
    icoords[0].bmatp_create();
    icoords[2-1].bmatp_create();
    icoords[0].bmatp_to_U();
    icoords[2-1].bmatp_to_U();
    icoords[0].bmat_create();
    icoords[2-1].bmat_create();
#else
  for (int n=0;n<2;n++){
		icoords[n].bmatp_create();
    icoords[n].bmatp_to_U();
    icoords[n].bmat_create();
		icoords[n].print_q();}
#endif
		


}



