#include <stdio.h>

void srchndpushc_(
	int *ix,int *iy,int *iz,int *MAXNOD,int *MAXXYZ,int *MXYZ,int *ncnode,
	int *nctx,int *ncty,int *nctz,int *nlooku,int *nodmax,int *nx,int *ny,int *nz,int *nxy, 
	int *noct,float *UNEST,float *xmn,float *ymn,float *zmn,float *xsiz,float *ysiz,float *zsiz,
	int icnode[],int ixnode[],int iynode[],int iznode[],
	float cnodex[],float cnodey[],float cnodez[],float cnodev[],float sim[],
	int cnodeid[]
)
{
	int i,j,k,il,ind,ncsec;
	int basex,basey,basez; 
	int nnx,nny,nnz,nnodmax,nnlooku,nncnode,nnxy;
	float UUNEST; 
/*
c
c Consider all the nearby nodes until enough have been found:
c
*/
	nncnode = 0;
	ncsec  = 0;
	basex=*ix-*nctx-1;
	basey=*iy-*ncty-1;
	basez=*iz-*nctz-1;
	nnx=*nx;
	nny=*ny;
	nnz=*nz;
	nnxy=nnx*nny;
	nnodmax=*nodmax;
	nnlooku=*nlooku;
	UUNEST=*UNEST;

	for(il=1;il<nnlooku;il++){

		if(nncnode==nnodmax){ 
			*ncnode=nncnode;
			return;
		}

		i = ixnode[il] + basex;
		if(i<1 || i>nnx) continue;
		j = iynode[il] + basey;
		if(j<1 || j>nny) continue;
		k = iznode[il] + basez;
		if(k<1 || k>nnz) continue;

/*
c
c Check this potentially informed grid node:
c
*/
		ind = (k-1)*(nnxy) + (j-1)*(nnx) + i;
		if(sim[ind-1]>UUNEST){
			nncnode         = nncnode + 1;
			cnodeid[nncnode-1] = ind;
			icnode[nncnode-1] = il+1;
		}
	}
/*
c
c Return to calling program:
c
*/
      return;
}


