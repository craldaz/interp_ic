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


void Interp::com_rotate_move(int iR, int iP, int iN, double ff) {
 
 printf(" in com_rotate_move() \n");

  ic3.reset(natoms,anames,anumbers,icoords[iN].coords);
	int nn=4;
  double* xyz0 = new double[natoms*3];
  double* xyz1 = new double[natoms*3];
  double* xyz2 = new double[natoms*3];

// displace center of mass
  double mfrac = 0.5;
  if (2-nn+1!=1)
    mfrac = 1./(2-nn+1);
  mfrac *= ff;
  //printf(" dXYZ frac: %1.2f \n",mfrac);

  for (int i=0;i<3*natoms;i++) xyz0[i] = icoords[iR].coords[i];
  for (int i=0;i<3*natoms;i++) xyz2[i] = icoords[iP].coords[i];
  for (int i=0;i<3*natoms;i++) xyz1[i] = 0.;
  double mx0 = 0.;
  double my0 = 0.;
  double mz0 = 0.;
  double mx2 = 0.;
  double my2 = 0.;
  double mz2 = 0.;
  double mass = 0.;

  //for (int i=0;i<natoms;i++)
  //  printf(" amasses[%i]: %1.2f \n",i,amasses[i]);

    for (int i=0;i<natoms;i++)
    {
      mass += amasses[i];
      mx0 += amasses[i] * xyz0[3*i];
      my0 += amasses[i] * xyz0[3*i+1];
      mz0 += amasses[i] * xyz0[3*i+2];
    }
    mx0 = mx0 / mass; my0 = my0 / mass; mz0 = mz0 / mass;
    for (int i=0;i<natoms;i++)
    {
      mx2 += amasses[i] * xyz2[3*i];
      my2 += amasses[i] * xyz2[3*i+1];
      mz2 += amasses[i] * xyz2[3*i+2];
    }
    mx2 = mx2 / mass; my2 = my2 / mass; mz2 = mz2 / mass;

   // printf(" CoM1: %1.3f %1.3f %1.3f \n",mx0,my0,mz0);
   // printf(" CoM2: %1.3f %1.3f %1.3f \n",mx2,my2,mz2);
    double mx1 = mfrac*(mx2-mx0);
    double my1 = mfrac*(my2-my0);
    double mz1 = mfrac*(mz2-mz0);
    for (int i=0;i<natoms;i++)
    {
      ic3.coords[3*i] += mx1;
      ic3.coords[3*i+1] += my1;
      ic3.coords[3*i+2] += mz1;
    }

 // printf("\n  doing rotation \n");
  int natomsqm = natoms;
  double* amassesqm = new double[natomsqm+1];
  int c = 0;
  for (int i=0;i<natoms;i++)
  {
    xyz1[3*c]   = newic.coords[3*i];
    xyz1[3*c+1] = newic.coords[3*i+1];
    xyz1[3*c+2] = newic.coords[3*i+2];
    xyz2[3*c]   = icoords[iP].coords[3*i];
    xyz2[3*c+1] = icoords[iP].coords[3*i+1];
    xyz2[3*c+2] = icoords[iP].coords[3*i+2];
 //   amassesqm[c] = amasses[i];
    amassesqm[c] = 1.0;
    c++;
  }  
 
#if 0
  cout << " " << natoms << endl << endl;
  for (int i=0;i<natoms;i++)
    cout << anames[i] << " " << newic.coords[3*i+0] << " "  << newic.coords[3*i+1] << " " << newic.coords[3*i+2] << endl;
#endif

  Eckart::Eckart_align(xyz2,xyz1,amassesqm,natomsqm,mfrac);
//  Eckart::Eckart_align(xyz2,xyz1,amassesqm,natomsqm);
 
  c = 0;
  for (int i=0;i<natoms;i++)
  {
    ic3.coords[3*i]   = xyz1[3*c];
    ic3.coords[3*i+1] = xyz1[3*c+1];
    ic3.coords[3*i+2] = xyz1[3*c+2];
    c++;
  }  
#if 0
  cout << " " << natoms << endl << endl;
  for (int i=0;i<natoms;i++)
    cout << anames[i] << " " << newic.coords[3*i+0] << " "  << newic.coords[3*i+1] << " " << newic.coords[3*i+2] << endl;
#endif

  delete [] amassesqm;


  delete [] xyz0;
  delete [] xyz1;
  delete [] xyz2;

  icoords[iN].reset(natoms,anames,anumbers,newic.coords);

  return;
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

  ICoord ic1,ic2,ic3,ic4; 
  ic1.alloc(natoms);
  ic2.alloc(natoms);
  ic3.alloc(natoms);
  ic1.reset(natoms,anames,anumbers,coords[0]);
  ic3.reset(natoms,anames,anumbers,coords[0]);
  ic2.reset(natoms,anames,anumbers,coords[2-1]);
  ic1.ic_create();
  ic2.ic_create();

  allcoords = new double*[2];
  for (int i=0;i<2;i++)
    allcoords[i] = icoords[i].coords;

#if 1
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

  for (int n=0;n<2;n++){
		icoords[n].bmatp_create();
    icoords[n].bmatp_to_U();
    icoords[n].bmat_create();
		icoords[n].print_q();}
  
	int size_ic = ic3.nbonds + ic3.nangles + ic3.ntor;
  //double** ictan = new double*[2];
  //for (int i=0;i<2;i++)
 		//ictan[i] = new double[size_ic+100];
 	double* ictan = new double[size_ic+100]; 
  double* ictan0 = new double[size_ic];
#endif
#if 0
// create union_ic
  newic.alloc(natoms);
  intic.alloc(natoms);
  int2ic.alloc(natoms);
  newic.reset(natoms,anames,anumbers,icoords[0].coords);
  intic.reset(natoms,anames,anumbers,icoords[1].coords);
  int2ic.reset(natoms,anames,anumbers,icoords[1].coords);
  newic.union_ic(ic1,ic2);  
  intic.copy_ic(newic);
  printf("I'm here\n");
  int2ic.copy_ic(newic);

  printf("\n actual IC's \n");
  newic.print_ic();

  newic.bmat_alloc();
  newic.bmatp_create();
  newic.bmatp_to_U();
  newic.bmat_create();
  intic.bmat_alloc();
  intic.bmatp_create();
  intic.bmatp_to_U();
  intic.bmat_create();
  int2ic.bmat_alloc();
  int2ic.bmatp_create();
  int2ic.bmatp_to_U();
  int2ic.bmat_create();

  int size_ic = newic.nbonds + newic.nangles + newic.ntor;
  double** ictan = new double*[2];
  for (int i=0;i<2;i++)
 		ictan[i] = new double[size_ic+100];
#endif

#if 0 
cout << "About to interpolate to halfway pt." << endl; 
	//double* qhalf = new double[len_d];  
	ic4.alloc(natoms);
	ic4.bmat_alloc();
	for (int i=0;i<len_d;i++)
		{ 
		ic3.q[i] = 0.5*icoords[0].q[i] + 0.5*icoords[1].q[i]; 
	  cout << ic3.q[i] << " " ; 
		}
	ic3.update_ic();	
	ic3.ic_to_xyz();
	ic3.print_xyz();
	cout << "\n done" << endl; 
#endif	
	

  int nbonds = ic3.nbonds;
  int nangles = ic3.nangles;
  int ntor = ic3.ntor;
  int len_d = ic3.nicd0;

 //full redundant tangent
  for (int i=0;i<nbonds;i++)
    ictan[i] = ic1.bondd[i] - ic2.bondd[i];
  for (int i=0;i<nangles;i++)
    ictan[nbonds+i] = (ic1.anglev[i] - ic2.anglev[i])*3.14159/180.;
  for (int i=0;i<ntor;i++)
  {
    ictan[nbonds+nangles+i] = (ic1.torv[i] - ic2.torv[i])*3.14159/180.;
    if (ictan[nbonds+nangles+i]>3.14159)
      ictan[nbonds+nangles+i] = -1*(2*3.14159 - ictan[nbonds+nangles+i]);
    if (ictan[nbonds+nangles+i]<-3.14159)
      ictan[nbonds+nangles+i] = 2*3.14159 + ictan[nbonds+nangles+i];
  }

#if 1
    printf(" printing ictan \n");
    for (int i=0;i<nbonds;i++)
      printf(" %1.2f",ictan[i]);
    printf("\n");
    for (int i=0;i<nangles;i++)
      printf(" %1.2f",ictan[nbonds+i]);
    printf("\n");
    if (ntor>0)
    for (int i=0;i<ntor;i++)
      printf(" %1.2f",ictan[nbonds+nangles+i]);
    printf("\n");
#endif

    double dqmag = 0.;
    for (int i=0;i<size_ic;i++) ictan0[i] = ictan[i];

    ic3.opt_constraint(ictan);

    for (int j=0;j<size_ic-ntor;j++)
      dqmag += ictan0[j]*ic3.Ut[ic3.nicd*size_ic+j];
    for (int j=nbonds+nangles;j<size_ic;j++)
      dqmag += ictan0[j]*ic3.Ut[ic3.nicd*size_ic+j]; 

    printf(" dqmag: %1.2f",dqmag);

    ic3.dq0[ic3.nicd0-1] = -dqmag/2;
    printf(" dq0[constraint]: %1.2f \n",ic3.dq0[ic3.nicd0-1]);
		
		com_rotate_move(1,0,1,1.0);
    int success = ic3.ic_to_xyz();
  return;
}

