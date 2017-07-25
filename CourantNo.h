!/*
!Global
!    CourantNo
!
!Description
!    Calculates and outputs the mean and maximum Courant Numbers.
!
!// ************************************************************************* //*/

 if (ltransient) then
 CoNum = 0.0
 meanCoNum = 0.0

!/*
!if (mesh.nInternalFaces())
!{
!    scalarField sumPhi
!    (
!        fvc::surfaceSum(mag(phi))().internalField()
!    );
!
!    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
!
!    meanCoNum =
!        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
!}
!*/

 su = 0.
 suma = 0.
 do k=3,nkmm; do i=3,nimm; do j=3,njmm
      inp=lk(k)+li(i)+j

      su(inp) = abs( f1(inp) ) + abs( f1(inp-nj) ) + &
                abs( f2(inp) ) + abs( f2(inp-1) ) + &
                abs( f3(inp) ) + abs( f3(inp-nij) )

      CoNum = max( CoNum , su(inp)/Vol(inp) )

      suma = suma + Vol(inp)

      meanCoNum = meanCoNum + su(inp)

 enddo; enddo; enddo

    CoNum = 0.5*CoNum*timestep
    meanCoNum = 0.5*meanCoNum/suma*timestep

  !// If we keep the value of Courant Number fixed
  if(CoNumFix) then
      dt = timestep
      timestep = CoNumFixValue/CoNum * timestep

      CoNum = CoNum * timestep/dt
      meanCoNum  = meanCoNum * timestep/dt
  endif

  time = time + timestep

 write(66,*)
 write(66,'(a,i0,a,es10.3,a,f12.6)') " Time step no. : ",ITIME," dt : ",timestep," Time = ",TIME
 write(66,*)

 write(66,'(2(a,es10.3))') "Courant Number mean: ", meanCoNum, &
&                          " max: ", CoNum



endif
!// ************************************************************************* //
