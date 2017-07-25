!***********************************************************************
!
  SUBROUTINE GRAD_LSQ_QR(FI,GRADFI,ISTAGE)
!
!***********************************************************************
!
!      Purpose:
!      Calculates cell-centered gradients using Least-Squares approach.
!
!      Description:
!      Uses QR decomposition of system matrix via Householder or via
!      Gramm-Schmidt.
!      QR decomposition is precomputed and R^(-1)*Q^T is stored in 
!      Dqr array for every cell.
!
!      FI - field variable which gradient we look for.
!      GradFi - cell centered gradient - a three component gradient vector.
!      ISTAGE - integer. If ISTAGE=1 calculates and stores only geometrical
!      parameters - a system matrix for least square problem at every cell. 
!      Usually it is called with ISTAGE=1 at the beggining of simulation.
!      If 2 it doesn't calculate system matrix, just RHS and solves system.
!
!      Dqr - least-squares matrix, calculated once at the beggining.
!      NXYZA - Integer the length of FI array in the calling routine.
!      NX - Length of LI
!      NZ - Length of LK
!      NKMM,NIMM,NJMM - number of inner cells in Z,X and Y direction respectively.
!      LK - How many cells is there below K-th floor.
!      LI - How many cells is there left for I-th vertical cross-section.
!      IDEW,IDNS,IDTB - Integer that can bring you to E,W,N,S,T,B cells from current cell P.
!      XC,YC,ZC - coordinates of cell centers 
!
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use gradients
  use matrix_module

  IMPLICIT NONE

  INTEGER, PARAMETER :: m=6, n=3 ! m is the number of neighbours, for structured 3D mesh it's 6

  INTEGER,  INTENT(IN) :: ISTAGE
  REAL(dp), DIMENSION(NXYZA), INTENT(IN)   :: FI
  REAL(dp), DIMENSION(n,NXYZA), INTENT(INOUT) :: GRADFI

!
!    Locals
!
  INTEGER ::  I,J,K,INP
  INTEGER :: ine, inw, inn, ins, int, inb

  REAL(dp), DIMENSION(m,n) :: Dtmp
  REAL(dp), DIMENSION(m)   :: b

  !REAL(dp), DIMENSION(m,n) :: R
  !REAL(dp), DIMENSION(m,m) :: Q
  !REAL(dp), DIMENSION(n,n) :: R1
  !REAL(dp), DIMENSION(n,m) :: Q1t

  INTEGER :: INFO
  INTEGER, DIMENSION(n) :: WORK
  REAL(dp), DIMENSION(n) :: TAU
  REAL(dp), DIMENSION(m) :: v1,v2,v3
  REAL(dp), DIMENSION(m,m) :: H1,H2,H3,Ieye
  REAL(dp), DIMENSION(n,n) :: R
  REAL(dp), DIMENSION(m,m) :: Q
      

  do k=3,nkmm
  do i=3,nimm
  do j=3,njmm

  inp=lk(k)+li(i)+j

!     Indexes of east, south, north, bottom, and top cell
  ine = inp + nj
  inw = inp - nj
  inn = inp + 1
  ins = inp - 1 
  int = inp + nij
  inb = inp - nij

  if(istage.eq.1) then

  Dtmp = 0.0d0

! Coefficient matrix - should be calculated only once 

  Dtmp(1,1) = xc(ine)-xc(inp)
  Dtmp(1,2) = yc(ine)-yc(inp)
  Dtmp(1,3) = zc(ine)-zc(inp)


  Dtmp(2,1) = xc(inw)-xc(inp)
  Dtmp(2,2) = yc(inw)-yc(inp)
  Dtmp(2,3) = zc(inw)-zc(inp)


  Dtmp(3,1) = xc(inn)-xc(inp)
  Dtmp(3,2) = yc(inn)-yc(inp)
  Dtmp(3,3) = zc(inn)-zc(inp)


  Dtmp(4,1) = xc(ins)-xc(inp)
  Dtmp(4,2) = yc(ins)-yc(inp)
  Dtmp(4,3) = zc(ins)-zc(inp)


  Dtmp(5,1) = xc(int)-xc(inp)
  Dtmp(5,2) = yc(int)-yc(inp)
  Dtmp(5,3) = zc(int)-zc(inp)


  Dtmp(6,1) = xc(inb)-xc(inp)
  Dtmp(6,2) = yc(inb)-yc(inp)
  Dtmp(6,3) = zc(inb)-zc(inp)


!1 ...Decompose A=QR using Householder
!      call householder_qr(Dtmp, m, n, Q, R)
!2 ...Decompose A=QR using Gram-Schmidt
!      call mgs_qr(Dtmp, m, n, Q, R)

!      Q = transpose(Q)
!      Q1t = Q(1:n,1:m)     ! NOTE: A=Q1R1 is so-called 'thin QR factorization' - see Golub & Van Loan
                            ! Here Q1 is actually Q1^T a transpose of Q1(thin Q - Q with m-n column stripped off)
!      R1 = R(1:n,1:n)      ! our Q1 is thin transpose of original Q.
!      R1 = inv(R1)         ! inv is a function in matrix_module, now works only for 3x3 matrices.
!      Q1t  = matmul(R1,Q1t) ! this is actually R^(-1)*Q^T - a matrix of size n x m.
!      D(:,:,INP) = Q1t     ! Store it for later.

!3....LAPACK routine DGEQRF
  CALL DGEQRF( M, N, Dtmp, M, TAU, WORK, N, INFO )

  ! Upper triangular matrix R
  R(1:n,1:n)=Dtmp(1:n,1:n)

  ! Create reflectors
  !H(i) = I - TAU * v * v'
  Ieye=eye(6)
  !v(1:i-1) = 0. and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i)
  v1(1) = 1.; v1(2:m)=Dtmp(2:m,1)
  H1 = rank_one_update(Ieye,m,m,v1,v1,-TAU(1))
  v2(1) = 0.; v2(2) = 1.; v2(3:m)=Dtmp(3:m,2)
  H2 = rank_one_update(Ieye,m,m,v2,v2,-TAU(2))
  v3(1:2) = 0.; v3(3) = 1.; v3(4:m)=Dtmp(4:m,3)
  H3 = rank_one_update(Ieye,m,m,v3,v3,-TAU(3))
  ! The matrix Q is represented as a product of elementary reflectors H1, H2, ..., Hn
  Q=matmul(H1,H2)
  Q=matmul(Q,H3)

  ! Form R_1^(-1)*Q_1^T explicitely.
  Dqr(1,1,INP) = Q(1,1)/R(1,1)-(R(1,2)*Q(1,2))/(R(1,1)*R(2,2))+(Q(1,3)*(R(1,2)*R(2,3)-R(1,3)*R(2,2)))/(R(1,1)*R(2,2)*R(3,3))
  Dqr(1,2,INP) = Q(2,1)/R(1,1)-(R(1,2)*Q(2,2))/(R(1,1)*R(2,2))+(Q(2,3)*(R(1,2)*R(2,3)-R(1,3)*R(2,2)))/(R(1,1)*R(2,2)*R(3,3))
  Dqr(1,3,INP) = Q(3,1)/R(1,1)-(R(1,2)*Q(3,2))/(R(1,1)*R(2,2))+(Q(3,3)*(R(1,2)*R(2,3)-R(1,3)*R(2,2)))/(R(1,1)*R(2,2)*R(3,3))
  Dqr(1,4,INP) = Q(4,1)/R(1,1)-(R(1,2)*Q(4,2))/(R(1,1)*R(2,2))+(Q(4,3)*(R(1,2)*R(2,3)-R(1,3)*R(2,2)))/(R(1,1)*R(2,2)*R(3,3)) 
  Dqr(1,5,INP) = Q(5,1)/R(1,1)-(R(1,2)*Q(5,2))/(R(1,1)*R(2,2))+(Q(5,3)*(R(1,2)*R(2,3)-R(1,3)*R(2,2)))/(R(1,1)*R(2,2)*R(3,3))
  Dqr(1,6,INP) = Q(6,1)/R(1,1)-(R(1,2)*Q(6,2))/(R(1,1)*R(2,2))+(Q(6,3)*(R(1,2)*R(2,3)-R(1,3)*R(2,2)))/(R(1,1)*R(2,2)*R(3,3))

  Dqr(2,1,INP) = Q(1,2)/R(2,2)-(R(2,3)*Q(1,3))/(R(2,2)*R(3,3))
  Dqr(2,2,INP) = Q(2,2)/R(2,2)-(R(2,3)*Q(2,3))/(R(2,2)*R(3,3))
  Dqr(2,3,INP) = Q(3,2)/R(2,2)-(R(2,3)*Q(3,3))/(R(2,2)*R(3,3)) 
  Dqr(2,4,INP) = Q(4,2)/R(2,2)-(R(2,3)*Q(4,3))/(R(2,2)*R(3,3))
  Dqr(2,5,INP) = Q(5,2)/R(2,2)-(R(2,3)*Q(5,3))/(R(2,2)*R(3,3))
  Dqr(2,6,INP) = Q(6,2)/R(2,2)-(R(2,3)*Q(6,3))/(R(2,2)*R(3,3))

  Dqr(3,1,INP) = Q(1,3)/R(3,3)                                                                             
  Dqr(3,2,INP) = Q(2,3)/R(3,3)                                                                             
  Dqr(3,3,INP) = Q(3,3)/R(3,3)                                                                             
  Dqr(3,4,INP) = Q(4,3)/R(3,3)                                                                             
  Dqr(3,5,INP) = Q(5,3)/R(3,3)                                                                             
  Dqr(3,6,INP) = Q(6,3)/R(3,3)

  elseif(istage.eq.2) then

  b=0.0d0

!     RHS vector
  b(1) = FI(ine)-FI(inp)
  b(2) = FI(inw)-FI(inp)
  b(3) = FI(inn)-FI(inp)
  b(4) = FI(ins)-FI(inp)
  b(5) = FI(int)-FI(inp)
  b(6) = FI(inb)-FI(inp)

!     Solve overdetermined system in least-sqare sense
!  ...using precomputed QR factorization and storing R^(-1)*Q^T in D
  GRADFI(1,INP) = sum(Dqr(1,1:m,inp)*b(1:m))
  GRADFI(2,INP) = sum(Dqr(2,1:m,inp)*b(1:m))
  GRADFI(3,INP) = sum(Dqr(3,1:m,inp)*b(1:m))

  endif

  ENDDO
  ENDDO
  ENDDO

  RETURN
  END
