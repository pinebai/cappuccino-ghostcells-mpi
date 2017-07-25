      ! Initialize ths vectors
      su = 0.0_dp
      sv = 0.0_dp 
      sw = 0.0_dp


      ! Cell-face loop:
    
      ! East Cell-Face
      call fluxuvw_explicit(nj,1,nij, &
                            ar1x,ar1y,ar1z, &
                            fx,f1)
      ! North Cell-Face
      call fluxuvw_explicit(1,nij,nj, &
                            ar2x,ar2y,ar2z, &
                            fy,f2)
      ! Top   Cell-Face
      call fluxuvw_explicit(nij,nj,1, &
                            ar3x,ar3y,ar3z, &
                            fz,f3)
