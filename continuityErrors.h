! Global
!     continuityErrs
!
! Description
!     Calculates and prints the continuity errors.
!---------------------------------------------------------------------------


!//     volScalarField contErr(fvc::div(phi)); ili sada radim volScalarField contErr(fvc::surfaceSum(phi))
!      // idew=nj
!      // idns=1
!      // idtb=nij

!      // do k=2,nkm
!      // do i=2,nim
!      // do j=2,njm

!      // inp=lk(k)+li(i)+j

!      // inw = inp-idew
!      // ins = inp-idns
!      // inb = inp-idtb

!      //  !// sv(inp) = -f1(inp)*ar1x(inp)+f1(inw)*ar1x(inw) &
!      //  !//           -f2(inp)*ar2y(inp)+f2(ins)*ar2y(ins) &
!      //  !//           -f3(inp)*ar3z(inp)+f3(inb)*ar3z(inb)

!      // sv(inp) = -f1(inp)+f1(inp-idew) &
!      //           -f2(inp)+f2(inp-idns) &
!      //           -f3(inp)+f3(inp-idtb)

!      // enddo
!      // enddo
!      // enddo

!//     scalar sumLocalContErr = runTime.deltaTValue()*
!//         mag(contErr)().weightedAverage(mesh.V()).value();
      !//sumLocalContErr = volumeWeightedAverage(abs(sv))!timestep*


!//     scalar globalContErr = runTime.deltaTValue()*
!//         contErr.weightedAverage(mesh.V()).value();
       !//globalContErr = volumeWeightedAverage(sv)!timestep*

!//     cumulativeContErr += globalContErr;
       !// cumulativeContErr = cumulativeContErr + globalContErr


     !//  write(66,'(3(a,es10.3))') "time step continuity errors : sum local = ", sumLocalContErr, &
     !// &                          ", global = ", globalContErr, &
     !// &                          ", cumulative = ", cumulativeContErr

!// ************************************************************************* //

      idew=nj
      idns=1
      idtb=nij

      LocalContErr = 0.0d0
      globalContErr = 0.0d0
      sumLocalContErr = 0.0d0

      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      LocalContErr = -f1(inp)+f1(inp-idew) &
                     -f2(inp)+f2(inp-idns) &
                     -f3(inp)+f3(inp-idtb)

      globalContErr = globalContErr + LocalContErr

      sumLocalContErr = sumLocalContErr + abs(LocalContErr)

      enddo
      enddo
      enddo

      cumulativeContErr = cumulativeContErr + globalContErr

      write(66,'(3(a,es10.3))') "time step continuity errors : sum local = ", sumLocalContErr, &
     &                          ", global = ", globalContErr, &
     &                          ", cumulative = ", cumulativeContErr

      !// write(66,'(20x,a,1pe10.3,1x,a,1pe10.3)') ' sum  =',globalContErr,'|sum| =',sumLocalContErr
