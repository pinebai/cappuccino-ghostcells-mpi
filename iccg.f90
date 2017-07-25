!***********************************************************************
!
      subroutine iccg(fi,ifi)
!
!***********************************************************************
!
!    This routine incorporates the incomplete cholesky preconditioned 
!    conjugate gradient solver for symmetric matrices in 3d problems
!    with seven-diagonal matrix structure (see sect. 5.3.6). array
!    index ijk converted from 3d indices i, j, and k according to
!    table 3.1. ns is the number of inner iterations (sweeps).
!
!    Writen by ismet demirdzic, sarajevo, 1991.
!    Modified by nikola mirkov, 28.01.2014. nmirkov@vinca.rs
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use coef
      use title_mod
      use variables

      implicit none
!
!***********************************************************************
!
      integer, intent(in) :: ifi
      real(prec), dimension(nxyza) :: fi 

!
!     local variables
!
      integer :: i, j, k, ijk, ns, l
      real(prec), dimension(nxyza) :: pk,zk,d
      real(prec) :: rsm, resmax, res0, resl
      real(prec) :: s0, sk, alf, bet, pkapk

!.....max no. of inner iters
      resmax = sor(ifi)
!
!.....initalize working arrays
!
      do ijk=1,nijk
        pk(ijk)=0.0d0
        zk(ijk)=0.0d0
        d(ijk)=0.0d0
        res(ijk)=0.0d0
      end do
!
!.....calculate initial residual vector and the norm
!
      res0=0.0d0
      do k=3,nkmm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
             res(ijk)=ae(ijk)*fi(ijk+nj)+aw(ijk)*fi(ijk-nj)+an(ijk)* &
             fi(ijk+1)+as(ijk)*fi(ijk-1)+at(ijk)*fi(ijk+nij)+ &
             ab(ijk)*fi(ijk-nij)+su(ijk)-ap(ijk)*fi(ijk)
            res0=res0+abs(res(ijk))
          end do
        end do
      end do
!
!.....if ltest=true, print the norm 
!
      if(ltest) write(66,'(a,1pe10.3)') '                    res0 = ',res0
!
!.....calculate elements of diagonal preconditioning matrix
!
      do k=3,nkmm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            d(ijk)=1./(ap(ijk)-aw(ijk)**2*d(ijk-nj)-as(ijk)**2*d(ijk-1)  &
                   -ab(ijk)**2*d(ijk-nij))
          end do
        end do
      end do
!
      s0=1.e20
!
!....start inner iterations
!
      ns=nsw(ifi)
      do l=1,ns
!
!.....solve for zk(ijk) -- forward substitution
!
      do k=3,nkmm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            zk(ijk)=(res(ijk)+aw(ijk)*zk(ijk-nj)+as(ijk)*zk(ijk-1)+  &
                    ab(ijk)*zk(ijk-nij))*d(ijk)
          end do
        end do
      end do

      do k=3,nkmm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            zk(ijk)=zk(ijk)/(d(ijk)+small)
          end do
        end do
      end do
!
!..... backward substitution; calculate scalar product sk
!
      sk=0.0d0
      do k=nkmm,3,-1
        do i=nimm,3,-1
          do j=njmm,3,-1
            ijk=lk(k)+li(i)+j
            zk(ijk)=(zk(ijk)+ae(ijk)*zk(ijk+nj)+an(ijk)*zk(ijk+1)+ &
                     at(ijk)*zk(ijk+nij))*d(ijk)
            sk=sk+res(ijk)*zk(ijk)
          end do
        end do
      end do
!
!.....calculate beta
!
      bet=sk/s0
!
!.....calculate new search vector pk
!
      do k=3,nkmm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            pk(ijk)=zk(ijk)+bet*pk(ijk)
          end do
        end do
      end do
!
!.... calculate scalar product (pk.a pk) and alpha (overwrite zk)
!
      pkapk=0.0d0
      do k=3,nkmm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            zk(ijk)=ap(ijk)*pk(ijk)-ae(ijk)*pk(ijk+nj)-                &
              aw(ijk)*pk(ijk-nj)-an(ijk)*pk(ijk+1)-as(ijk)*pk(ijk-1)-  &
              at(ijk)*pk(ijk+nij)-ab(ijk)*pk(ijk-nij)
            pkapk=pkapk+pk(ijk)*zk(ijk)
          end do
        end do
      end do

      alf=sk/pkapk
!
!.....calculate variable correction, new residual vector, and norm
!
      resl=0.0d0
      do k=3,nkmm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            fi(ijk)=fi(ijk)+alf*pk(ijk)
            res(ijk)=res(ijk)-alf*zk(ijk)
            resl=resl+abs(res(ijk))
          end do
        end do
      end do

      s0=sk
!
!.....check convergence
!
      if(l.eq.1) resor(ifi)=res0
      rsm=resl/(resor(ifi)+small)
      if(ltest) write(66,'(19x,3a,i4,a,1pe10.3,a,1pe10.3)') ' fi=',chvar(ifi),' sweep = ',l,' resl = ',resl,' rsm = ',rsm
      if(rsm.lt.resmax) exit
!
!.....end of iteration loop
!
      end do

!.....Write linear solver report:
      write(66,'(3a,1PE10.3,a,1PE10.3,a,I0)') &
      'PCG(IC0):  Solving for ',trim(chvarSolver(IFI)), &
      ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L
!
      return
      end
