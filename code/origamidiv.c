//-----------------------------------------------------------------------------
// Copyright (c) 2016, Bridget L. Falck, Mark C. Neyrinck, & Nuala McCullagh
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "params.h"

#define DL for (d=0;d<3;d++)
#define BF 1e30
#define max(A,B) (((A)>(B)) ? (A):(B))
#define goodmod(A,B) (((A) >= (B)) ? (A-B):(((A) < 0) ? (A+B):(A)))
int isneg(int h) {
    return (int)(h < 0);
}
int par(int i, int j, int k, int ng) {
    return i + (j + k*ng)*ng;
}

void read_parameters(char *filename);
long divread(char *filename, float ***p, float fact, float boxsize, int np1d, int npd, int ndiv, int divid, int nbuf);
int readgadget(char *filename, float ***p);


int main(int argc, char *argv[]) {
    int exitcode;
    long np;
    float **r;

    FILE *pos, *tag;
    char *paramfile,tagoutfile[80];
    float xmin,xmax,ymin,ymax,zmin,zmax;
  
    int s;
    float negb2,b2;
    int ng2,ng4, np2,np3, h, i,i2, j, x,y,z,nhalo,nhalo0,nhalo1,nhalo2,nhaloany, ndiv,divid,npd;
    unsigned char *m,*m0,*m1,*m2,*M, mn,m0n,m1n,m2n; /*Morphology tag */
    int zstart, ystart, xstart, zend, yend, xend;
    float dx,d1,d2;

    if (argc < 2) {
      printf("Please specify a parameter file\n");
      exit(0);
    } else paramfile = argv[1];
    read_parameters(paramfile);

    if (sscanf(argv[2],"%d",&ndiv) != 1) {
        printf("That's no ndiv; try again.\n");
        exit(0);
    }
    if (sscanf(argv[3],"%d",&divid) != 1) {
        printf("That's no divid; try again.\n");
        exit(0);
    }
    if (np1d%ndiv != 0) {
        printf("ndiv must be factor of np1d.\n");
	exit(0);
    }
    npd = np1d/ndiv;

    b2 = boxsize/2.;
    negb2 = -boxsize/2.;

    sprintf(tagoutfile,"%s%s_%dtag.dat",outdir,taglabel,divid);

    np2 = npd+nbuf+nbuf;
    np3 = np2*np2*np2;

    ng2=nbuf*2;
    ng4=nbuf;

    if (numfiles > 0) {
      np = readgadget(posfile,&r);
    } else np = divread(posfile,&r,1.,boxsize,np1d,npd,ndiv,divid,nbuf);
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

	if (i<(npd*npd*npd)) {
	    M[i] = 1;
	}
    }

    if (m==NULL) {
        printf("Morphology array cannot be allocated.\n");
        exit(0);
    }
    //  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);
    printf("Calculating ORIGAMI morphology.\n");
#pragma omp parallel for default(none) shared(np1d, ng4, r, boxsize, negb2, b2, m, m0, m1, m2, nsplit) private (x, y, z, h, dx, i, i2, d1, d2, xstart, xend, ystart, yend, zstart, zend)
    for (s = 0; s < nsplit*nsplit*nsplit; s++) {
        zstart = s / (nsplit*nsplit);
        zend = zstart + 1;
        ystart = (s - zstart * nsplit * nsplit)/nsplit;
        yend = ystart + 1;
        xstart = (s - zstart * nsplit * nsplit - ystart * nsplit);
        xend = xstart + 1;
        
        zstart *= (np1d/nsplit);
        ystart *= (np1d/nsplit);
        xstart *= (np1d/nsplit);
        zend *= (np1d/nsplit);
        yend *= (np1d/nsplit);
        xend *= (np1d/nsplit);
        for (x=xstart; x<xend; x++){
            //    printf("%d\n",x);fflush(stdout);
            for (y=ystart; y<yend; y++) {
                for (z=zstart; z<zend; z++) {
                    i = par(x,y,z,np1d);
                    /* First just along the Cartesian axes */
                    /* x-direction */
                    for (h=1; abs(h)<ng4; h = -h + isneg(h)) {
                        i2 = par(goodmod(x+h,np1d),y,z,np1d);
                        dx = r[i2][0]-r[i][0];
                        if (dx < negb2) dx += boxsize;
                        if (dx > b2) dx -= boxsize;
                        if (dx*h < 0.) {
                            if (m[i] % 2 > 0) {
                                m[i] *= 2;
                            }
                            break;
                        }
                    }
                    for (h=1; abs(h)<ng4; h = -h + isneg(h)) {
                        i2 = par(x,goodmod(y+h, np1d),z,np1d);
                        dx = r[i2][1]-r[i][1];
                        if (dx < negb2) dx += boxsize;
                        if (dx > b2) dx -= boxsize;
                        if (dx*h < 0.) {
                            /*printf("y:%d %d %d %d %f\n",x,y,z,h,dx);*/
                            if (m[i] % 3 > 0) {
                                m[i] *= 3;
                            }
                            break;
                        }
                    }
                    for (h=1; abs(h)<ng4; h = -h + isneg(h)) {
                        i2 = par(x,y,goodmod(z+h, np1d),np1d);
                        dx = r[i2][2]-r[i][2];
                        if (dx < negb2) dx += boxsize;
                        if (dx > b2) dx -= boxsize;
                        if (dx*h < 0.) {
                            /*printf("z:%d %d %d %d %f\n",x,y,z,h,dx);*/
                            if (m[i] % 5 > 0) {
                                m[i] *= 5;
                            }
                            break;
                        }
                    }
                    // Now do diagonal directions 
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(x,goodmod(y+h,np1d),goodmod(z+h,np1d),np1d);
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
                        i2 = par(x,goodmod(y+h,np1d),goodmod(z-h,np1d),np1d);
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
                        i2 = par(goodmod(x+h,np1d),y,goodmod(z+h,np1d),np1d);
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
                        i2 = par(goodmod(x+h,np1d),y,goodmod(z-h,np1d),np1d);
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
                        i2 = par(goodmod(x+h,np1d),goodmod(y+h,np1d),z,np1d);
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
                        i2 = par(goodmod(x+h,np1d),goodmod(y-h,np1d),z,np1d);
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
	x = (i%np2)-nbuf;
	if (x < 0 || x >= npd) {
	    continue;
	}
	y = ((i/np2)%np2)-nbuf;
	if (y < 0 || y >= npd) {
	    continue;
	}
	z = (i/(np2*np2))-nbuf;
	if (z < 0 || z >= npd) {
	    continue;
	}
	M[j] = m[i];
	++j;
    }
    printf("nhalo=%d,%d,%d,%d,%d\n",nhalo,nhalo0,nhalo1,nhalo2,nhaloany);

    /* Output m */
    printf("Writing morphology file to %s.\n",tagoutfile);fflush(stdout);
    tag = fopen(tagoutfile,"w");
    if (tag == NULL) {
        printf("Unable to open %s\n",tagoutfile);
        exit(0);
    }
    fwrite(&np,1, sizeof(int),tag);
    fwrite(M,npd*npd*npd,sizeof(unsigned char),tag);

    return(1);
}
