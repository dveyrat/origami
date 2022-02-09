//-----------------------------------------------------------------------------
// Copyright (c) 2014, Bridget L. Falck & Mark C. Neyrinck
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file COPYING.txt, distributed with this software.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#define DL for (d=0;d<3;d++) /* Dimension loop */
#define BF 1e30

/* Positions */
/* Returns number of particles read */
int posread(char *posfile, float ***p, float fact, int nps, int npd, int ndiv, int divid, int nbuf) {

  FILE *pos;
  int np,dum,d,i,j,np2,np3, a,b,c,x,y,z;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float *ptemp;

  printf("fact=%f\n",fact);
  fflush(stdout);

  pos = fopen(posfile, "r");
  if (pos == NULL) {
    printf("Unable to open position file %s\n\n",posfile);
    exit(0);
  }
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Read number of particles */
  /*fread(&dum,1,4,pos); */
  fread(&np,1, sizeof(int),pos); 
  /*fread(&dum,1,4,pos);*/
  np2 = npd+nbuf+nbuf;
  np3 = np2*np2*np2;

  /* Allocate the arrays */
  (*p) = (float **)malloc(np3*sizeof(float *));
  ptemp = (float *)malloc(np*sizeof(float));

  printf("np = %d\n",np);

  j = 0;
  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) {
    if (((i/(nps*nps))+nbuf)%nps >= np2+(npd*(divid/(ndiv*ndiv)))) {
      continue;
    }
    if ((((i/nps)%nps)+nbuf)%nps >= np2+(npd*((divid/ndiv)%ndiv))) {
      continue;
    }
    if (((i%nps)+nbuf)%nps >= np2+(npd*(divid%ndiv))) {
      continue;
    }
    if (((i/(nps*nps))+nbuf)%nps < npd*(divid/(ndiv*ndiv))) {
      continue;
    }
    if ((((i/nps)%nps)+nbuf)%nps < npd*((divid/ndiv)%ndiv)) {
      continue;
    }
    if (((i%nps)+nbuf)%nps < npd*(divid%ndiv)) {
      continue;
    }
    (*p)[j] = (float *)malloc(3*sizeof(float));
    if ((*p)[j] == NULL) {
      printf("Unable to allocate particle array in readfiles!\n");
      fflush(stdout);
      exit(0);
    }
    (*p)[j][0] = ptemp[i];
    ++j;
  }
  /*fread(&dum,1,4,pos); 
    fread(&dum,1,4,pos); */
  j = 0;
  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) {
    a = ;
    b = ;
    c = ;
    a = ;
    b = ;
    x = i%nps;
    y = (i/nps)%nps;
    z = i/(nps*nps);
    if (((i/(nps*nps))+nbuf)%nps >= np2+(npd*(divid/(ndiv*ndiv)))) {
      continue;
    }
    if ((((i/nps)%nps)+nbuf)%nps >= np2+(npd*((divid/ndiv)%ndiv))) {
      continue;
    }
    if (((i%nps)+nbuf)%nps >= np2+(npd*(divid%ndiv))) {
      continue;
    }
    if (((i/(nps*nps))+nbuf)%nps < npd*(divid/(ndiv*ndiv))) {
      continue;
    }
    if ((((i/nps)%nps)+nbuf)%nps < npd*((divid/ndiv)%ndiv)) {
      continue;
    }
    if (((i%nps)+nbuf)%nps < npd*(divid%ndiv)) {
      continue;
    }
    (*p)[j][1] = ptemp[i];
    ++j;
  }
  /*fread(&dum,1,4,pos);
    fread(&dum,1,4,pos); */
  j = 0;
  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) {
    if (((i/(nps*nps))+nbuf)%nps >= np2+(npd*(divid/(ndiv*ndiv)))) {
      continue;
    }
    if ((((i/nps)%nps)+nbuf)%nps >= np2+(npd*((divid/ndiv)%ndiv))) {
      continue;
    }
    if (((i%nps)+nbuf)%nps >= np2+(npd*(divid%ndiv))) {
      continue;
    }
    if (((i/(nps*nps))+nbuf)%nps < npd*(divid/(ndiv*ndiv))) {
      continue;
    }
    if ((((i/nps)%nps)+nbuf)%nps < npd*((divid/ndiv)%ndiv)) {
      continue;
    }
    if (((i%nps)+nbuf)%nps < npd*(divid%ndiv)) {
      continue;
    }
    (*p)[j][2] = ptemp[i];
    ++j;
  }
  /*fread(&dum,1,4,pos); */

  fclose(pos);
  free(ptemp);

  /* Get into physical units (Mpc/h) */

  for (i=0; i<np3; i++) DL (*p)[i][d] *= fact;


  /* Test range -- can comment out */
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np3;i++) {
    if ((*p)[i][0]<xmin) xmin = (*p)[i][0]; if ((*p)[i][0]>xmax) xmax = (*p)[i][0];
    if ((*p)[i][1]<ymin) ymin = (*p)[i][1]; if ((*p)[i][1]>ymax) ymax = (*p)[i][1];
    if ((*p)[i][2]<zmin) zmin = (*p)[i][2]; if ((*p)[i][2]>zmax) zmax = (*p)[i][2];
  }
  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);

  return(np);
}

/* Velocities */
/* Returns number of particles read */
int velread(char *velfile, float ***v, float fact) {

  FILE *vel;
  int np,dum,d,i;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float *vtemp;

  vel = fopen(velfile, "r");
  if (vel == NULL) {
    printf("Unable to open velocity file %s\n\n",velfile);
    exit(0);
  }
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Read number of particles */
  /*fread(&dum,1,4,vel);*/ fread(&np,1, sizeof(int),vel); /*fread(&dum,1,4,vel);*/
  printf("np=%d\n",np); fflush(stdout);

  /* Allocate the arrays */
  (*v) = (float **)malloc(np*sizeof(float *));
  vtemp = (float *)malloc(np*sizeof(float));

  printf("about to read velocity array\n"); fflush(stdout);

  /* Fill the arrays */
  /*fread(&dum,1,4,vel);*/
  fread(vtemp,np,4,vel);
  for (i=0; i<np; i++) {
    (*v)[i] = (float *)malloc(3*sizeof(float));
    if ((*v)[i] == NULL) {
      printf("Unable to allocate particle array in readfiles!\n");
      fflush(stdout);
      exit(0);
    }
    (*v)[i][0] = vtemp[i];
  }
  /*fread(&dum,1,4,vel); 
    fread(&dum,1,4,vel); */
  fread(vtemp,np,4,vel);
  for (i=0; i<np; i++) (*v)[i][1] = vtemp[i];
  /*fread(&dum,1,4,vel);
    fread(&dum,1,4,vel); */
  fread(vtemp,np,4,vel);
  for (i=0; i<np; i++) (*v)[i][2] = vtemp[i];
  /*fread(&dum,1,4,vel); */
  fclose(vel);

  /* Convert from code units into physical units (km/sec) */
  
  printf("np=%d\n",np); fflush(stdout);
  for (i=0; i<np; i++) DL (*v)[i][d] *= fact;

  /* Test range -- can comment out */
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if ((*v)[i][0] < xmin) xmin = (*v)[i][0]; if ((*v)[i][0] > xmax) xmax = (*v)[i][0];
    if ((*v)[i][1] < ymin) ymin = (*v)[i][1]; if ((*v)[i][1] > ymax) ymax = (*v)[i][1];
    if ((*v)[i][2] < zmin) zmin = (*v)[i][2]; if ((*v)[i][2] > zmax) zmax = (*v)[i][2];
  }
  printf("vx: %f,%f; vy: %f,%f; vz: %f,%f\n",xmin,xmax, ymin,ymax, zmin,zmax);fflush(stdout);
  return(np);
}
