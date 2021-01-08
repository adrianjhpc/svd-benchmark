
#ifndef MYSCALAPACK_H
#define MYSCALAPACK_H
#include <mkl_scalapack.h>
#include <mkl.h>

#ifdef SUP
#define pdsymm_   pdsymm
#define pdmatgen_ pdmatgen
#define pdtrsm_   pdtrsm
#define psgemm_   psgemm
#define pdgemm_   pdgemm
#define numroc_   numroc
#define pselset_  pselset
#define pdelset_  pdelset
#define indxg2p_  indxg2p
#define indxg2l_  indxg2l
#define descinit_ descinit
#define pslawrite_ pslawrite
#define pdlawrite_ pdlawrite
#define blacs_get_      blacs_get
#define blacs_pinfo_    blacs_pinfo
#define blacs_gridinit_ blacs_gridinit
#define blacs_gridinfo_ blacs_gridinfo
#define blacs_gridexit_ blacs_gridexit
#define blacs_exit_     blacs_exit
#endif

extern "C" {
 void Cblacs_pinfo( int* mypnum, int* nprocs);
 void Cblacs_get( int context, int request, int* value);
 int  Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
 void Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
 void Cblacs_gridexit( int context);
 void Cblacs_exit( int error_code);

 void blacs_pinfo_( int *mypnum, int *nprocs);
 void blacs_get_( int *context, int *request, int* value);
 void blacs_gridinit_( int* context, char *order, int *np_row, int *np_col);
 void blacs_gridinfo_( int *context, int *np_row, int *np_col, int *my_row, int *my_col);
 void blacs_gridexit_( int *context);
 void blacs_exit_( int *error_code);

 void pdtrsm_ ( char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, double *a, int *ia,
		int *ja, int *desca, double *b, int *ib, int *jb, int *descb );

 void psgemm_( char *transa, char *transb, int *M, int *N, int *K,
                                          float     *alpha,
                                          float     *A, int *ia, int *ja, int *descA,
                                          float     *B, int *ib, int *jb, int *descB,
                                          float     *beta,
                                          float     *C, int *ic, int *jc, int *descC );
 void pdgemm_( char *transa, char *transb, int *M, int *N, int *K,
                                          double    *alpha,
                                          double    *A, int *ia, int *ja, int *descA,
                                          double    *B, int *ib, int *jb, int *descB,
                                          double    *beta,
                                          double    *C, int *ic, int *jc, int *descC );


 void pselset_( float     *A, int *ia, int *ja, int *descA, float     *alpha);
 void pdelset_( double    *A, int *ia, int *ja, int *descA, double    *alpha);

 void pslawrite_( char **filenam, int *m, int *n, float  *A, int *ia, int *ja, int *descA, int *irwrit, int *icwrit, float  *work);
 void pdlawrite_( char **filenam, int *m, int *n, double *A, int *ia, int *ja, int *descA, int *irwrit, int *icwrit, double *work);

 int indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
 int indxg2l_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
 int numroc_( MKL_INT *n, MKL_INT *nb, MKL_INT *iproc, MKL_INT *isrcproc, MKL_INT *nprocs);
 void descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc,
				int *ictxt, int *lld, MKL_INT *info);

 void   pdmatgen_( int *ictxt, char *aform, char *diag, int *m, int *n, int *mb, int *nb, double *a, int *lda, int *iarow, int *iacol, int *iseed, int *iroff, int *irnum, int *icoff, int *icnum, int *myrow, int *mycol, int *nprow, int *npcol );

 void pdsymm_(char *side, char *uplo, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, double *beta, double *c, int *ic, int *jc, int *descc);
 void chk1mat_(int *ma, int *mapos0, int *na, int *napos0, int *ia,
  int *ja, int *desca, int *descapos0, int *info);
 void pchk1mat_(int *ma, int *mapos0, int *na, int *napos0, int *ia,
  int *ja, int *desca, int *descapos0, int *nextra, int *ex, int *expos,
  int *info);
  
}
#endif
