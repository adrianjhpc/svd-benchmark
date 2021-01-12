#include "mpi.h"
#ifdef MKL
#include "mkl_blacs.h"
#include "mkl.h"
#endif 

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/SVD>

#define ROOT 0

extern "C"
{
void descinit_(int *, int *, int *, int *, int *, 
               int *, int *, int *, int *, int *);
int numroc_(int *, int *, int *, int *, int *);
int indxl2g_(int *, int *, int *, int *, int *);

void pdgesvd_(char*, char*, int*, int*,double*,int*,int*,int*,double*,
                                               double*,int*,int*,int*,
                                               double*,int*,int*,int*,
                                               double*,int*,int*);
void pzgesvd_(char*, char*, int*, int*,double*,int*,int*,int*,double*,
                                               double*,int*,int*,int*,
                                               double*,int*,int*,int*,
                                               double*,int*,double*,int*);

#ifndef MKL
void blacs_pinfo_(int *, int *);

void blacs_abort_(int *, int *);

void  blacs_get_(int *, int *, int *);

void  blacs_gridinit_(int *, char *, int *, int *);

void  blacs_gridinfo_(int *, int *, int *, int *, int *);

void  blacs_barrier_(int *, char *);

void  blacs_gridexit_(int *);

void  blacs_exit_(int *);

void  dgesd2d_(int *, int *, int *, double *, int *, int *, int *);

void  dgerv2d_(int *, int *, int *, double *, int *, int *, int *);
#endif 
}
int scatter_matrix(int rank, int M, int N, int mb, int nb, int mp, int np, int nprow, int npcol, int myrow, int mycol, Eigen::VectorXd *A, Eigen::MatrixXd *TotalA, int ictxt, const int dtype);
int gather_matrix(int rank, int M, int N, int mb, int nb, int mp, int np, int nprow, int npcol, int myrow, int mycol, Eigen::VectorXd *A, Eigen::MatrixXd *TotalA, int ictxt, const int dtype);

const int   RSRC_ = 6;
const int   CSRC_ = 7;

int main(int argc, char** argv){
#ifdef MKL
   const int len = 198;
   char buf[198];
   MKL_Get_Version_String(buf,len);
   std::cout << "MKL Version: " << buf << "\n";
#endif

   char R_letter = 'R';
   char A_letter = 'A';
   char N_letter = 'N';

   int rank(0), nprocs(1);

   double start, end;

   blacs_pinfo_( &rank, &nprocs );

   const int dtype = 1;
   int i_one  = 1;
   int i_zero = 0;
   int ictxt, ictxt_some, ictxt_all;
   int myrow(0),mycol(0);
   int M = 16;
   int N = 16;
   int procrow = 4;
   int proccol = 2;
   blacs_get_( &i_zero, &i_zero, &ictxt );


   if(argc != 5){
     if(rank == ROOT){
       std:: cout << "Need to pass in the arguments M and N which are the matrix size and the number of rows and columns for the proc grid." << std::endl;
     }
     blacs_abort_(&ictxt, 0);
   }

   M = atoi(argv[1]);
   N = atoi(argv[2]);
   procrow = atoi(argv[3]);
   proccol = atoi(argv[4]);

   int nprow(procrow), npcol(proccol);

   if(rank == ROOT){
     std::cout << "Matrix size is " << M << " x " << N << std::endl;
   }
   
   if(nprocs < nprow * npcol){
     printf("Not using enough MPI processes, using %d processes, but the blacs grid is definied as %d x %d\n", nprocs, nprow, npcol);
     blacs_abort_(&ictxt, 0);
   }else{
     if(rank == ROOT){
       std::cout<< "Process grid is " << nprow << " x " << npcol << " (" << nprocs-(nprow*npcol) << " are idle)" << std::endl;
     }
   }

   ictxt_all = ictxt;
   
   blacs_gridinit_(&ictxt_all, &R_letter, &nprocs, &i_one);
   
   ictxt_some = ictxt;
   
   blacs_gridinit_(&ictxt_some, &R_letter, &nprow, &npcol);
   
   int minMN = std::min(M,N);
   std::vector<double> s(minMN*dtype);
   Eigen::MatrixXd TotalA(M,N);
   Eigen::VectorXd TotalS(minMN*dtype);
   
   TotalA.setZero();

   if(ictxt_some >= 0){

     blacs_gridinfo_( &ictxt_some, &nprow, &npcol, &myrow, &mycol );

     

     int mb = std::min(M/nprocs,N/nprocs);
     int nb = std::min(M/nprocs,N/nprocs);

     int mp = numroc_( &M, &mb, &myrow, &i_zero, &nprow); 
     int np = numroc_( &N, &nb, &mycol, &i_zero, &npcol); 

     int info(0);
     std::vector<double> a(mp*np*dtype,0.0);

     std::vector<double> pwork(mb*dtype*5);
     Eigen::VectorXd A(mp*np*dtype);
     Eigen::VectorXd newA(mp*np*dtype);
     A.setZero();

     double u[0*dtype], vt[0*dtype];
     int desca[9], descu[9], descvt[9];

     int lld = std::max(i_one,mp);
     descinit_(desca, &M, &N, &mb, &nb, &i_zero, &i_zero, &ictxt_some, &lld, &info); 

     std::vector<double> work(1*dtype);
     std::vector<double> rwork(1);
     int lwork = -1;

     // Create the distributed data initially
     for (int j = 1; j <= np; ++j){
       int gj = indxl2g_(&j,&nb,&mycol,&desca[CSRC_],&npcol);
       for (int i = 1; i <= mp; ++i){
         int gi = indxl2g_(&i,&mb,&myrow,&desca[RSRC_],&nprow);

         if(gi != gj) continue;
         A[dtype*(lld*(j-1)+(i-1))] = gi;
       }
     }

     start = MPI_Wtime();

     gather_matrix(rank, M, N, mb, nb, mp, np, nprow, npcol, myrow, mycol, &A, &TotalA, ictxt_some, dtype);

     end = MPI_Wtime();

     if(rank ==  ROOT){
            std::cout << "Blacs gather time is: " << end - start << std::endl;
     }

     blacs_barrier_(&ictxt_some, &A_letter);

     start = MPI_Wtime();

     scatter_matrix(rank, M, N, mb, nb, mp, np, nprow, npcol, myrow, mycol, &newA, &TotalA, ictxt_some, dtype);

     end = MPI_Wtime();

     if(rank ==  ROOT){
            std::cout << "Blacs scatter time is: " << end - start << std::endl;
     }

     for(int i=0; i<np; i++){
       for(int j=0; j<mp; j++){
	 if(A[i*mp+j] != newA[i*mp+j]){
	   printf("Problem comparing A and newA %d %d %d\n",rank,i,j);
	 }
       }
     }

     if(2 == dtype){
       pzgesvd_(&N_letter,&N_letter,&M,&N,&newA[0],&i_one,&i_one,desca,&s[0],
		u,&i_one,&i_one,descu,
		vt,&i_one,&i_one,descvt,
		&work[0],&lwork,&rwork[0],&info);
     }else{
        pdgesvd_(&N_letter,&N_letter,&M,&N,&newA[0],&i_one,&i_one,desca,&s[0],
         u,&i_one,&i_one,descu,
         vt,&i_one,&i_one,descvt,
         &work[0],&lwork,&info);
     }
     lwork = static_cast<int>(work[0]);
     if(2 == dtype){
       rwork.resize(static_cast<int>(rwork[0]));
     }
     if(ROOT == rank){
       printf(" lwork : %d\n", lwork);
     }
     work.resize(lwork*dtype);


     start = MPI_Wtime();
     if(2 == dtype){
       if(rank == ROOT){
	 printf("pzgesvd\n");
       }
       pzgesvd_(&N_letter,&N_letter,&M,&N,&newA[0],&i_one,&i_one,desca,&s[0],
		u,&i_one,&i_one,descu,
		vt,&i_one,&i_one,descvt,
		&work[0],&lwork,&rwork[0],&info);
     }else{
       if(rank == ROOT){
	 printf("pdgesvd\n");
       }
       pdgesvd_(&N_letter,&N_letter,&M,&N,&newA[0],&i_one,&i_one,desca,&s[0],
		u,&i_one,&i_one,descu,
		vt,&i_one,&i_one,descvt,
		&work[0],&lwork,&info);
     }
     if (ROOT == rank){ 
       std::cout << "Finished SVD.\n";
     }
     blacs_gridexit_(&ictxt_some);
   }
   blacs_barrier_(&ictxt_all, &A_letter);
   blacs_gridexit_(&ictxt_all);
   end = MPI_Wtime();

   if(ROOT == rank){
     std::cout << "Scalapack time is: " << end - start << std::endl;
     start = MPI_Wtime();
     Eigen::BDCSVD<Eigen::MatrixXd> svd(TotalA,Eigen::ComputeFullU|Eigen::ComputeFullV);
     end = MPI_Wtime();
     std::cout << "Eigen time is: " << end - start << std::endl;

     TotalS = svd.singularValues();
     int correct = 1;
     for (int i = 0; i < minMN; ++i){
       if(s[i] != TotalS[i]){
	 std::cout << i << " : " << s[i] << " : "  << TotalS[i] << "\n";
       	 correct = 0;
       }

     }

     if(correct == 0){
       std::cout << "Error in comparison betwen Eigen and Scalapack" << std::endl;
     }

   }

   blacs_exit_(&i_zero);
   
   return 0;
}


int scatter_matrix(int rank, int M, int N, int mb, int nb, int mp, int np, int nprow, int npcol, int myrow, int mycol, Eigen::VectorXd *A, Eigen::MatrixXd *TotalA, int ictxt, const int dtype){


  int i_zero = 0;
  int send_row = 0;
  int send_col = 0;
  int recv_row = 0;
  int recv_col = 0;
  int datasize = 0;

  for(int row = 0; row < M; row += mb, send_row=(send_row+1)%nprow){
    send_col = 0;
    int nr = mb;
    if(N-row < nb){
      nr = N-row;
    }
    
    for(int col = 0; col < N; col += nb, send_col=(send_col+1)%npcol){
      int nc = nb;
      if (N-col < nb){
	nc = N-col;
	}
      
      if(ROOT == rank){
	datasize = M;
	dgesd2d_(&ictxt, &nr, &nc, &(*TotalA).data()[(M*col+row)], &datasize, &send_row, &send_col);
      }
      
      if(myrow == send_row && mycol == send_col){ 
	datasize = mp;
	dgerv2d_(&ictxt, &nr, &nc, (double *)&(*A).data()[(mp*recv_col+recv_row)], &datasize, &i_zero, &i_zero);
	recv_col = (recv_col+nc)%np;
      }
    }
    
    if(myrow == send_row){
      recv_row = (recv_row+nr)%mp;
    }
    
  }
  
  return 0;
  
}



int gather_matrix(int rank, int M, int N, int mb, int nb, int mp, int np, int nprow, int npcol, int myrow, int mycol, Eigen::VectorXd *A, Eigen::MatrixXd *TotalA, int ictxt, const int dtype){
  
  
  int i_zero = 0;
  int send_row = 0;
  int send_col = 0;
  int recv_row = 0;
  int recv_col = 0;
  int datasize = 0;
  
  for(int row = 0; row < M; row += mb, send_row=(send_row+1)%nprow){
    send_col = 0;
    int nr = mb;
    if(N-row < nb){
      nr = N-row;
    }

    for(int col = 0; col < N; col += nb, send_col=(send_col+1)%npcol){
      int nc = nb;
      if (N-col < nb){
	nc = N-col;
      }
      
      if(myrow == send_row && mycol == send_col){
	datasize = mp;
	dgesd2d_(&ictxt, &nr, &nc, &(*A).data()[(mp*recv_col+recv_row)], &datasize, &i_zero, &i_zero);
	recv_col = (recv_col+nc)%np;
      }
      
      if(ROOT == rank){
	datasize = M;
	dgerv2d_(&ictxt, &nr, &nc, (double *)&(*TotalA).data()[(M*col+row)], &datasize, &send_row, &send_col);
      }
    }
    
    if(myrow == send_row){
      recv_row = (recv_row+nr)%mp;
    }
    
  }
  
  return 0;
  
}
