//-----------------------------------------------------------------------------
// Copyright (c) 2014, Bridget L. Falck & Mark C. Neyrinck
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file COPYING.txt, distributed with this software.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define DL for (d=0;d<3;d++)
#define BF 1e30
//#define NG 256
#define max(A,B) (((A)>(B)) ? (A):(B))
#define goodmod(A,B) (((A) >= (B)) ? (A-B):(((A) < 0) ? (A+B):(A)))
#define CATCH 209923

int isneg(int h) {
  return (int)(h < 0);
}

int posread(char *posfile, float ***p, float fact, int nps, int npd, int ndiv, int divid, int nbuf);

int par(int i, int j, int k, int ng) {
  return i + (j + k*ng)*ng;
}

int main(int argc, char *argv[]) {
  int exitcode;
  int np;
  float **r;

  time_t start_t,end_t;
  double diff_t;

  FILE *pos, *tag;
  char *posfile, tagoutfile[80], *label, *outdir;
  float predict, xmin,xmax,ymin,ymax,zmin,zmax;
  
  int isitinbuf;
  char isitinmain, d;
  int numdiv, nvp, nvpall, nvpbuf;
  float width, width2, totwidth, totwidth2, bf, s, g;
  float border, boxsize, negb2,b2;
  float c[3];
  int b[3];
  double totalvol;
  int ng, ng2,ng4, np2,np3, h, i,i2, j, x,y,z,nhalo,nhalo0,nhalo1,nhalo2,nhaloany,nps,npd,ndiv,divid,nbuf;
  unsigned char *m,*m0,*m1,*m2, mn,m0n,m1n,m2n, *M; /*Morphology tag */
  
  float dx,d1,d2;


  //if (argc != 6) {
  //  printf("Wrong number of arguments.\n");
  //  printf("arg1: position file\n");
  //  printf("arg2: label describing this run\n");
  //  printf("arg3: output directory\n");
  //  printf("arg4: box size\n");
  //  printf("arg5: 1D grid number, ngrid (nparticles = ngrid^3)\n\n");
  //  exit(0);
  //}
  posfile = argv[1];
  label = argv[2];
  outdir = argv[3];
  if (sscanf(argv[4],"%f",&boxsize) != 1) {
    printf("That's no boxsize; try again.\n");
    exit(0);
  }
  b2 = boxsize/2.;
  negb2 = -boxsize/2.;
  if (sscanf(argv[5],"%d",&ng) != 1) {
    printf("That's no ngrid; try again.\n");
    exit(0);
  }
  if (sscanf(argv[6],"%d",&nps) != 1) {
    printf("That's no ngrid; try again.\n");
    exit(0);
  }
  if (sscanf(argv[7],"%d",&npd) != 1) {
    printf("That's no ngrid; try again.\n");
    exit(0);
  }
  if (sscanf(argv[8],"%d",&ndiv) != 1) {
    printf("That's no ngrid; try again.\n");
    exit(0);
  }
  if (sscanf(argv[9],"%d",&divid) != 1) {
    printf("That's no ngrid; try again.\n");
    exit(0);
  }
  if (sscanf(argv[10],"%d",&nbuf) != 1) {
    printf("That's no ngrid; try again.\n");
    exit(0);
  }
  printf("%d\n",nps);
  printf("%d\n",npd);
  printf("%d\n",ndiv);
  printf("%d\n",divid);
  printf("%d\n",nbuf);
  //nps = argv[6];
  //npd = argv[7];
  //ndiv = argv[8];
  //divid = argv[9];
  //nbuf = argv[10];

  np2 = npd+nbuf+nbuf;
  np3 = np2*np2*np2;

  ng2=ng/2;
  ng4=ng/4;

  /* Boxsize should be the range in r, yielding a range 0-1 */
  np = posread(posfile,&r,1.,nps,npd,ndiv,divid,nbuf);
  printf("%d particles\n",np);fflush(stdout);
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  m = (unsigned char *)malloc(np3*sizeof(unsigned char));
  m0 = (unsigned char *)malloc(np3*sizeof(unsigned char)); /* for the diagonals */
  m1 = (unsigned char *)malloc(np3*sizeof(unsigned char));
  m2 = (unsigned char *)malloc(np3*sizeof(unsigned char));
  M = (unsigned char *)malloc(npd*npd*npd*sizeof(unsigned char));
  for (i=0; i<np3;i++) {
    if (r[i][0]<xmin) xmin = r[i][0]; if (r[i][0]>xmax) xmax = r[i][0];
    if (r[i][1]<ymin) ymin = r[i][1]; if (r[i][1]>ymax) ymax = r[i][1];
    if (r[i][2]<zmin) zmin = r[i][2]; if (r[i][2]>zmax) zmax = r[i][2];

    m[i] = 1;
    m0[i] = 1;
    m1[i] = 1;
    m2[i] = 1;

    //x = (i%np2)-nbuf;
    //y = ((i/np2)%np2)-nbuf;
    //z = (i/(np2*np2))-nbuf;

    if (i<(npd*npd*npd)) {
      //if (x>=0 && x<npd && y>=0 && y<npd &&) {
      M[i] = 1;
    }
  }

  if (m==NULL) {
    printf("Morphology array cannot be allocated.\n");
    exit(0);
  }
  //  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);

  printf("Calculating ORIGAMI morphology.\n");
  time(&start_t);
  for (x=0; x<ng; x++){
    //    printf("%d\n",x);fflush(stdout);
    for (y=0; y<ng; y++) {
      for (z=0; z<ng; z++) {
	i = par(x,y,z,ng);
	/* First just along the Cartesian axes */
	/* x-direction */
	for (h=1; h<ng4; h++) {
	  i2 = par((x+h)%ng,y,z,ng);
	  dx = r[i2][0]-r[i][0];
	  if (dx < negb2) dx += boxsize;
	  if (dx > b2) dx -= boxsize;
	  if (dx < 0.) {
	    /*printf("x:%d %d %d %d %f\n",x,y,z,h,dx);*/
	    if (m[i] % 2 > 0) {
	      m[i] *= 2;
	    }
	    if (m[i2] % 2 > 0){
	      m[i2] *= 2;
	    }
	    break;
	  }
	}
	for (h=1; h<ng4; h++) {
	  i2 = par(x,(y+h)%ng,z,ng);
	  dx = r[i2][1]-r[i][1];
	  if (dx < negb2) dx += boxsize;
	  if (dx > b2) dx -= boxsize;
	  if (dx < 0.) {
	    /*printf("y:%d %d %d %d %f\n",x,y,z,h,dx);*/
	    if (m[i] % 3 > 0) {
	      m[i] *= 3;
	    }
	    if (m[i2] % 3 > 0){
	      m[i2] *= 3;
	    }
	    break;
	  }
	}
	for (h=1; h<ng4; h++) {
	  i2 = par(x,y,(z+h)%ng,ng);
	  dx = r[i2][2]-r[i][2];
	  if (dx < negb2) dx += boxsize;
	  if (dx > b2) dx -= boxsize;
	  if (dx < 0.) {
	    /*printf("z:%d %d %d %d %f\n",x,y,z,h,dx);*/
	    if (m[i] % 5 > 0) {
	      m[i] *= 5;
	    }
	    if (m[i2] % 5 > 0){
	      m[i2] *= 5;
	    }
	    break;
	  }
	}
	// Now do diagonal directions 
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(x,goodmod(y+h,ng),goodmod(z+h,ng),ng);
	  d1 = r[i2][1]-r[i][1];
	  d2 = r[i2][2]-r[i][2];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 + d2)*h < 0.) {
	    m0[i] *= 2;
	    break;
	  }
	}
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(x,goodmod(y+h,ng),goodmod(z-h,ng),ng);
	  d1 = r[i2][1]-r[i][1];
	  d2 = r[i2][2]-r[i][2];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 - d2)*h < 0.) {
	    m0[i] *= 3;
	    break;
	  }
	}
	// y
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(goodmod(x+h,ng),y,goodmod(z+h,ng),ng);
	  d1 = r[i2][0]-r[i][0];
	  d2 = r[i2][2]-r[i][2];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 + d2)*h < 0.) {
	    m1[i] *= 2;
	    break;
	  }
	}
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(goodmod(x+h,ng),y,goodmod(z-h,ng),ng);
	  d1 = r[i2][0]-r[i][0];
	  d2 = r[i2][2]-r[i][2];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 - d2)*h < 0.) {
	    m1[i] *= 3;
	    break;
	  }
	}
	// z
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(goodmod(x+h,ng),goodmod(y+h,ng),z,ng);
	  d1 = r[i2][0]-r[i][0];
	  d2 = r[i2][1]-r[i][1];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 + d2)*h < 0.) {
	    m2[i] *=2;
	    break;
	  }
	}
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(goodmod(x+h,ng),goodmod(y-h,ng),z,ng);
	  d1 = r[i2][0]-r[i][0];
	  d2 = r[i2][1]-r[i][1];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 - d2)*h < 0.) {
	    m2[i] *= 3;
	    break;
	  }
	  }
      }
    }
  }

  nhalo = 0;
  nhalo0 = 0;
  nhalo1 = 0;
  nhalo2 = 0;
  nhaloany = 0;
  j = 0;
  for (i=0;i<np3;i++){
    mn = (m[i]%2 == 0) + (m[i]%3 == 0) + (m[i]%5 == 0);
    m0n = (unsigned char)(m[i]%2 == 0) + (unsigned char)(m0[i]%2 == 0) + (unsigned char)(m0[i]%3 == 0);
    m1n = (unsigned char)(m[i]%3 == 0) + (unsigned char)(m1[i]%2 == 0) + (unsigned char)(m1[i]%3 == 0);
    m2n = (unsigned char)(m[i]%5 == 0) + (unsigned char)(m2[i]%2 == 0) + (unsigned char)(m2[i]%3 == 0);
    m[i] = max(mn,max(m0n,max(m1n,m2n)));
    if (mn == 3) nhalo++;
    if (m0n == 3) nhalo0++;
    if (m1n == 3) nhalo1++;
    if (m2n == 3) nhalo2++;
    if (m[i] == 3) {
      nhaloany++;
    }
    if ((i/(np2*np2))-nbuf >= npd) {
      continue;
    }
    if ((((i/np2)%np2))-nbuf >= npd) {
      continue;
    }
    if ((i%np2)-nbuf >= npd) {
      continue;
    }
    if ((i/(np2*np2))-nbuf < 0) {
      continue;
    }
    if ((((i/np2)%np2))-nbuf < 0) {
      continue;
    }
    if ((i%np2)-nbuf < 0) {
      continue;
    }
    M[j] = m[i];
    ++j;
  }
  printf("nhalo=%d,%d,%d,%d,%d\n",nhalo,nhalo0,nhalo1,nhalo2,nhaloany);
  time(&end_t);
  diff_t = difftime(end_t, start_t);
  printf("Execution time = %f\n",diff_t);
  /* Output m */
  sprintf(tagoutfile,"%s%stag.dat",outdir,label);
  printf("Writing morphology file to %s.\n",tagoutfile);fflush(stdout);
  tag = fopen(tagoutfile,"w");
  if (tag == NULL) {
    printf("Unable to open %s\n",tagoutfile);
    exit(0);
  }
  //fwrite(&np3,1, sizeof(int),tag);
  fwrite(M,npd*npd*npd,sizeof(unsigned char),tag);

  return(1);
}
