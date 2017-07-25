  !***********************************************************************
  !
  subroutine output
  !
  !***********************************************************************
  !
  use types
  use parameters
  use indexes
  use variables
  use printing

  implicit none
  !
  !***********************************************************************
  !

  ! Printing of results
  if(lpri) then
  !  if(lcal(iu)) call printi(vis,'  viscos  ')
  if(lcal(iu)) call printi(u ,'u-velocity')
  if(lcal(iv)) call printi(v ,'v-velocity')
  if(lcal(iw)) call printi(w ,'w-velocity')
  if(lcal(ip)) call printi(p ,' pressure ')
  !  if(lcal(iu)) call printi(f1,' flux-1   ')
  !  if(lcal(iu)) call printi(f2,' flux-2   ')
  !  if(lcal(iu)) call printi(f3,' flux-3   ')
  if(lcal(ien)) call printi(t,'temperatur')
  endif

  if(lprj) then
  if(lcal(iu))   call printj(u  ,'u-velocity')
  if(lcal(iv))   call printj(v  ,'v-velocity')
  if(lcal(iw))   call printj(w  ,'w-velocity')
  if(lcal(ip))   call printj(p  ,' pressure ')
  if(lcal(iu))   call printj(f1 ,' flux-1   ')
  if(lcal(iu))   call printj(f2 ,' flux-2   ')
  if(lcal(iu))   call printj(f3 ,' flux-3   ')
  if(lcal(ite))  call printj(te ,'  te      ')
  if(lcal(ied))  call printj(ed ,'  ed      ')
  !  if(lcal(ite))  call printj(uu ,'  uu      ')
  !  if(lcal(ite))  call printj(vv ,'  vv      ')
  !  if(lcal(ite))  call printj(ww ,'  ww      ')
  !  if(lcal(ite))  call printj(uv ,'  uv      ')
  !  if(lcal(ite))  call printj(uw ,'  uw      ')
  !  if(lcal(ite))  call printj(vw ,'  vw      ')
  if(lcal(ivis)) call printj(vis,'  vis     ')
  if(lturb)      call printj(gen,'  gen     ')
  if(lturb)      call printj(ret,'  ret     ')
  if(lcal(ien))  call printj(t,'temperatur')
  !  if(lcal(ien))  call printj(ret,'ret number')
  !  if(lcal(ien))  call printj(dudx,'   dudx   ')
  !  if(lcal(ien))  call printj(dudy,'   dudy   ')
  !  if(lcal(ien))  call printj(dudz,'   dudz   ')
  !  if(lcal(ien))  call printj(utt,'   utt    ')
  !  if(lcal(ien))  call printj(vtt,'   vtt    ')
  !  if(lcal(ien))  call printj(wtt,'   wtt    ')
  !  if(lcal(ien))  call printj(vart,'  vart    ')
  !  if(lcal(icon))  call printj(con,'  con     ')
  endif

  if(lprk) then
  if(lcal(iu))   call printk(u  ,'u-velocity')
  if(lcal(iv))   call printk(v  ,'v-velocity')
  if(lcal(iw))   call printk(w  ,'w-velocity')
  if(lcal(ip))   call printk(p  ,' pressure ')
  !  if(lcal(iu))   call printk(f1 ,' flux-1   ')
  !  if(lcal(iu))   call printk(f2 ,' flux-2   ')
  !     if(lcal(iu))   call printk(f3 ,' flux-3   ')
  if(lcal(ite))  call printk(te ,'turb.kin.e')
  if(lcal(ied))  call printk(ed ,'  dissip. ')
  if(lcal(ivis)) call printk(vis,'viscosity ')
  if(lturb)      call printk(gen,'   g e n  ')
  if(lcal(ien))  call printk(t,'temperatur')
  !  if(lcal(ien))  call printk(dudx,'   dudx   ')
  !  if(lcal(ien))  call printk(dudy,'   dudy   ')
  !  if(lcal(ien))  call printk(dudy,'   dudz   ')
  endif
  
  return
  end
