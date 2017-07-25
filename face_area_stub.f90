

!+++++east cell faces+++++
      idew=nj
      idns=1
      idtb=nij
      do i=1,nim  
      do k=2,nkm
      do j=2,njm
      inp=lk(k)+li(i)+j
      ine=inp+idew
      ins=inp-idns
      inb=inp-idtb
      inbs=inb-idns
!.....second
      dxet=0.5*(x(inp)-x(ins)+x(inb)-x(inbs))
      dyet=0.5*(y(inp)-y(ins)+y(inb)-y(inbs))
      dzet=0.5*(z(inp)-z(ins)+z(inb)-z(inbs))
!.....third
      dxzd=0.5*(x(inp)-x(inb)+x(ins)-x(inbs))
      dyzd=0.5*(y(inp)-y(inb)+y(ins)-y(inbs))
      dzzd=0.5*(z(inp)-z(inb)+z(ins)-z(inbs))
!.....second x third define  always cf area
      ar1x(inp)=dyet*dzzd-dyzd*dzet
      ar1y(inp)=dxzd*dzet-dxet*dzzd
      ar1z(inp)=dxet*dyzd-dyet*dxzd
      enddo
      enddo
      enddo

!+++++north cell faces+++++
      idew=1
      idns=nij
      idtb=nj
      do j=1,njm 
      do k=2,nkm
      do i=2,nim
      inp=lk(k)+li(i)+j
      ine=inp+idew
      ins=inp-idns
      inb=inp-idtb
      inbs=inb-idns
!.....second
      dxet=0.5*(x(inp)-x(ins)+x(inb)-x(inbs))
      dyet=0.5*(y(inp)-y(ins)+y(inb)-y(inbs))
      dzet=0.5*(z(inp)-z(ins)+z(inb)-z(inbs))
!.....third
      dxzd=0.5*(x(inp)-x(inb)+x(ins)-x(inbs))
      dyzd=0.5*(y(inp)-y(inb)+y(ins)-y(inbs))
      dzzd=0.5*(z(inp)-z(inb)+z(ins)-z(inbs))
!.....second x third define  always cf area
      ar2x(inp)=dyet*dzzd-dyzd*dzet
      ar2y(inp)=dxzd*dzet-dxet*dzzd
      ar2z(inp)=dxet*dyzd-dyet*dxzd
      enddo
      enddo
      enddo

!+++++top cell faces+++++
      idew=nij
      idns=nj
      idtb=1
      do k=1,nkm
      do i=2,nim
      do j=2,njm
      inp=lk(k)+li(i)+j
      ine=inp+idew
      ins=inp-idns
      inb=inp-idtb
      inbs=inb-idns
!.....second
      dxet=0.5*(x(inp)-x(ins)+x(inb)-x(inbs))
      dyet=0.5*(y(inp)-y(ins)+y(inb)-y(inbs))
      dzet=0.5*(z(inp)-z(ins)+z(inb)-z(inbs))
!.....third
      dxzd=0.5*(x(inp)-x(inb)+x(ins)-x(inbs))
      dyzd=0.5*(y(inp)-y(inb)+y(ins)-y(inbs))
      dzzd=0.5*(z(inp)-z(inb)+z(ins)-z(inbs))
!.....second x third define  always cf area
      ar3x(inp)=dyet*dzzd-dyzd*dzet
      ar3y(inp)=dxzd*dzet-dxet*dzzd
      ar3z(inp)=dxet*dyzd-dyet*dxzd
      enddo
      enddo
      enddo
