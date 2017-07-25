! Global
!     continuityErrs
!
! Description
!     Calculates and prints the continuity errors.
!---------------------------------------------------------------------------


!//     volScalarField contErr(fvc::div(phi));
      idew=nj
      idns=1
      idtb=nij
      do k=2,nkm; do i=2,nim; do j=2,njm
      inp=lk(k)+li(i)+j
      inw = inp-idew
      ins = inp-idns
      inb = inp-idtb
      su(inp) = -f1(inp)*ar1x(inp)+f1(inw)*ar1x(inw) &
                -f2(inp)*ar2y(inp)+f2(ins)*ar2y(ins) &
                -f3(inp)*ar3z(inp)+f3(inb)*ar3z(inb)
      enddo; enddo; enddo

!//     scalar sumLocalContErr = runTime.deltaTValue()*
!//         mag(contErr)().weightedAverage(mesh.V()).value();
      sumLocalContErr = timestep*volumeWeightedAverage(abs(su))


!//     scalar globalContErr = runTime.deltaTValue()*
!//         contErr.weightedAverage(mesh.V()).value();
       globalContErr = timestep*volumeWeightedAverage(su)

!//     cumulativeContErr += globalContErr;
       cumulativeContErr = cumulativeContErr + globalContErr


      write(66,'(3(a,es10.3))') "time step continuity errors : sum local = ", sumLocalContErr, &
     &                          ", global = ", globalContErr, &
     &                          ", cumulative = ", cumulativeContErr

!// ************************************************************************* //
