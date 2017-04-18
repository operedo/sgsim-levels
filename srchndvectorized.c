#include <stdio.h>
#include <emmintrin.h>
#include <smmintrin.h>

void srchndvectorized_(
	int *ix,int *iy,int *iz,int *MAXNOD,int *MAXXYZ,int *MAXORD,int *MXYZ,int *ncnode,int *maxsec,
	int *nctx,int *ncty,int *nctz,int *nlooku,int *nodmax,int *nx,int *ny,int *nz, 
	float *UNEST,float *xmn,float *ymn,float *zmn,float *xsiz,float *ysiz,float *zsiz,
	int icnode[],int ixnode[],int iynode[],int iznode[],
	float cnodex[],float cnodey[],float cnodez[],float cnodev[],float cnodet[],float tmp[],float sim[],
	int cnodeid[]
)
{
	int i,j,k,v,il,indexloc,ncsec;
 
/*
c
c Consider all the nearby nodes until enough have been found:
c
*/
	*ncnode = 0;
	ncsec  = 0;

	__m128i ixvec = _mm_set1_epi32(*ix-*nctx-1);
	__m128i iyvec = _mm_set1_epi32(*iy-*ncty-1);
	__m128i izvec = _mm_set1_epi32(*iz-*nctz-1);

	__m128i nxvec = _mm_set1_epi32(*nx+1);
	__m128i nyvec = _mm_set1_epi32(*ny+1);
	__m128i nzvec = _mm_set1_epi32(*nz+1);

	__m128i zerovec = _mm_set1_epi32(0);
	__m128i minusvec = _mm_set1_epi32(-1);
	__m128i avec = _mm_set1_epi32((*nx)*(*ny));
	__m128i bvec   = _mm_set1_epi32(*nx);
	__m128i cvec = _mm_set1_epi32(-(*nx)*(*ny)-(*nx));

	__m128i ixnodevec, iynodevec, iznodevec;

	__m128i ivec,jvec,kvec;
	__m128i indexlocvec;

	//int iarr[4],jarr[4],karr[4];
	//int iboolarr[4],jboolarr[4],kboolarr[4];
	int indexlocarr[4];

	//for(il=0;il<*nlooku;il++){
	for(il=0;il<*nlooku-3;il=il+4){

		if(*ncnode==*nodmax) return;

		ixnodevec = _mm_add_epi32(_mm_loadu_si128((__m128i*)(ixnode + il)),ixvec);
		iynodevec = _mm_add_epi32(_mm_loadu_si128((__m128i*)(iynode + il)),iyvec);
		iznodevec = _mm_add_epi32(_mm_loadu_si128((__m128i*)(iznode + il)),izvec);

		ivec = _mm_and_si128(_mm_cmplt_epi32(zerovec,ixnodevec),_mm_cmplt_epi32(ixnodevec,nxvec));
		jvec = _mm_and_si128(_mm_cmplt_epi32(zerovec,iynodevec),_mm_cmplt_epi32(iynodevec,nyvec));
		kvec = _mm_and_si128(_mm_cmplt_epi32(zerovec,iznodevec),_mm_cmplt_epi32(iznodevec,nzvec));

		indexlocvec = 
			_mm_mullo_epi32(
				_mm_add_epi32(
					_mm_add_epi32(
						_mm_mullo_epi32(iznodevec,avec),
						_mm_mullo_epi32(iynodevec,bvec)
					),
					_mm_add_epi32(
						ixnodevec,
						cvec)
				),
				_mm_mullo_epi32(
					_mm_max_epi16(
						_mm_max_epi16(
							ivec,
							jvec
						),
						kvec
					),
					minusvec
				)
			);

 		_mm_storeu_si128((__m128i*)indexlocarr,indexlocvec);

 		//_mm_storeu_si128((__m128i*)iarr,ixnodevec);
 		//_mm_storeu_si128((__m128i*)jarr,iynodevec);
 		//_mm_storeu_si128((__m128i*)karr,iznodevec);
 		
		//_mm_store_si128((__m128i*)iboolarr,ivec);
 		//_mm_store_si128((__m128i*)jboolarr,jvec);
 		//_mm_store_si128((__m128i*)kboolarr,kvec);


//		i = *ix + (int)(ixnode[il])-*nctx-1;
//		if(i<1 || i>*nx) continue;
//		j = *iy + (int)(iynode[il])-*ncty-1;
//		if(j<1 || j>*ny) continue;
//		k = *iz + (int)(iznode[il])-*nctz-1;
//		if(k<1 || k>*nz) continue;

/*
c
c Check this potentially informed grid node:
c
*/

//		for(v=0;v<4;v++){
//			//printf("iarr[%d]=%d jarr[%d]=%d karr[%d]=%d\n",v,iboolarr[v],v,jboolarr[v],v,kboolarr[v]);
//			//if(iboolarr[v]<0 && jboolarr[v]<0 && kboolarr[v]<0){
//			if(indexlocarr[v]>0){
//				//indexloc = (karr[v]-1)*(*nx)*(*ny) + (jarr[v]-1)*(*nx) + iarr[v];
//				//indexloc = karr[v]*(*nx)*(*ny) + jarr[v]*(*nx) + iarr[v] - (*nx)*(*ny) - (*nx) ;
//				indexloc = indexlocarr[v] ;
//				if(sim[indexloc-1]>*UNEST){
//					*ncnode         = *ncnode + 1;
//					cnodeid[*ncnode-1] = indexloc;
//					icnode[*ncnode-1] = il+1+v;
//				}
//			}
//		}

		if(indexlocarr[0]>0){
			indexloc = indexlocarr[0] ;
			if(sim[indexloc-1]>*UNEST){
				*ncnode         = *ncnode + 1;
				cnodeid[*ncnode-1] = indexloc;
				icnode[*ncnode-1] = il+1+0;
			}
		}
		if(indexlocarr[1]>0){
			indexloc = indexlocarr[1] ;
			if(sim[indexloc-1]>*UNEST){
				*ncnode         = *ncnode + 1;
				cnodeid[*ncnode-1] = indexloc;
				icnode[*ncnode-1] = il+1+1;
			}
		}
		if(indexlocarr[2]>0){
			indexloc = indexlocarr[2] ;
			if(sim[indexloc-1]>*UNEST){
				*ncnode         = *ncnode + 1;
				cnodeid[*ncnode-1] = indexloc;
				icnode[*ncnode-1] = il+1+2;
			}
		}
		if(indexlocarr[3]>0){
			indexloc = indexlocarr[3] ;
			if(sim[indexloc-1]>*UNEST){
				*ncnode         = *ncnode + 1;
				cnodeid[*ncnode-1] = indexloc;
				icnode[*ncnode-1] = il+1+3;
			}
		}
	}
/*
c
c Return to calling program:
c
*/
      return;
}

void srchndvectorized2_(
	int *ix,int *iy,int *iz,int *MAXNOD,int *MAXXYZ,int *MAXORD,int *MXYZ,int *ncnode,int *maxsec,
	int *nctx,int *ncty,int *nctz,int *nlooku,int *nodmax,int *nx,int *ny,int *nz, 
	float *UNEST,float *xmn,float *ymn,float *zmn,float *xsiz,float *ysiz,float *zsiz,
	int icnode[],int ixnode[],int iynode[],int iznode[],
	float cnodex[],float cnodey[],float cnodez[],float cnodev[],float cnodet[],float tmp[],float sim[],
	int cnodeid[]
)
{
	int i,j,k,il,indexloc,ncsec;
	int basex,basey,basez; 
	int nnx,nny,nnz,nnodmax,nnlooku,nncnode;
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
	nnodmax=*nodmax;
	nnlooku=*nlooku;
	UUNEST=*UNEST;

	for(il=0;il<nnlooku;il++){

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
		indexloc = (k-1)*(nnx)*(nny) + (j-1)*(nnx) + i;
		if(sim[indexloc-1]>UUNEST){
			nncnode         = nncnode + 1;
			cnodeid[nncnode-1] = indexloc;
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

//void srchndvectorized2_(
//	int *ix,int *iy,int *iz,int *MAXNOD,int *MAXXYZ,int *MAXORD,int *MXYZ,int *ncnode,int *maxsec,
//	int *nctx,int *ncty,int *nctz,int *nlooku,int *nodmax,int *nx,int *ny,int *nz, 
//	float *UNEST,float *xmn,float *ymn,float *zmn,float *xsiz,float *ysiz,float *zsiz,
//	int icnode[],int ixnode[],int iynode[],int iznode[],
//	float cnodex[],float cnodey[],float cnodez[],float cnodev[],float cnodet[],float tmp[],float sim[],
//	int cnodeid[]
//)
//{
//	int i,j,k,il,indexloc,ncsec;
// 
///*
//c
//c Consider all the nearby nodes until enough have been found:
//c
//*/
//	*ncnode = 0;
//	ncsec  = 0;
//
//
//
//
//	for(il=0;il<*nlooku;il++){
//
//		if(*ncnode==*nodmax) return;
//
//		i = *ix + (int)(ixnode[il])-*nctx-1;
//		if(i<1 || i>*nx) continue;
//		j = *iy + (int)(iynode[il])-*ncty-1;
//		if(j<1 || j>*ny) continue;
//		k = *iz + (int)(iznode[il])-*nctz-1;
//		if(k<1 || k>*nz) continue;
//
///*
//c
//c Check this potentially informed grid node:
//c
//*/
//		indexloc = (k-1)*(*nx)*(*ny) + (j-1)*(*nx) + i;
//		if(sim[indexloc-1]>*UNEST){
//			*ncnode         = *ncnode + 1;
//			cnodeid[*ncnode-1] = indexloc;
//			icnode[*ncnode-1] = il+1;
//		}
//	}
///*
//c
//c Return to calling program:
//c
//*/
//      return;
//}
//
