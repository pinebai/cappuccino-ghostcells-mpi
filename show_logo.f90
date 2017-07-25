subroutine show_logo

write(66,'(a)')'................. ..............................................................'
write(66,'(a)')'..............MM........ MMM....................................................'
write(66,'(a)')'............NMN........... MM ..................................................'
write(66,'(a)')'           MMM..   OMM.    MMM.                                                 '
write(66,'(a)')'          .MMO..   MMM      MMM                         ...                     '
write(66,'(a)')'          MMM: .   MMM..    MMM                       .. .?DMD7....             '
write(66,'(a)')'          MMM~.    MMM  ....MMM                      ..MMMMMMMMMMM~.            '
write(66,'(a)')'          MMM:.. . .M~......MMM.                    .NMMMMMMMMMMMMMM.           '
write(66,'(a)')'          .MMM .....M..    .MMM.    ... .?MMMMM .. .MMMMMMMMMMMMMMMMM,..        '
write(66,'(a)')'           MMM... .,.M ..  MMM. ....MMMMMMMMMMMM$..MMMMMM?     .7MMMMM.         '
write(66,'(a)')'              MM.  M  M   MM  ..MMMMMMM~.. .MMMMM7MMMMMM.        7MMMMN         '
write(66,'(a)')'                  .  ...    =MMMMMM7..  ..MMM.MMMMMMMMM+         .MMMMM.        '
write(66,'(a)')'                         DMMMMMD..  . . MMMM. .MMMMMMMM.        ..MMMMM.        '
write(66,'(a)')'                    ..=MMMMMZ.     ...MMMM      MMMMMMMI.        .MMMM$.        '
write(66,'(a)')'                    MMMMM8.   ..  .NMMMM..      :MMMMMMM         :MMMM.         '
write(66,'(a)')'.................:MMMMM.......  ~MMMMN.......... MMMMMMMMO...... MMMM+..........'
write(66,'(a)')'             . ?MMMM:. ....  ,MMMMM...           .MMMMMMMMM.    :MMMM.          '
write(66,'(a)')'           ..IMMMM  .  .  =MMMMMI ..             ..MMMMMMMM.    MMMM=           '
write(66,'(a)')'............MMMM. .... DMMMMM? ....................8MMMMM=.....NMMMM............'
write(66,'(a)')'          +MMM.   .~MMMMMM~                        .MMMM.     ,MMMM.            '
write(66,'(a)')'.........~MM?.~MMMMMMM$..............................MMMD.....MMMMM ............'
write(66,'(a)')'        .DMMMMMMM$  ...                              MMMM.  .MMMMM..            '
write(66,'(a)')'         .MMMM.                                      .MMMN..MMMMM,              '
write(66,'(a)')'         ..MMMM.                                     .MMMM.MMMMMM..             '
write(66,'(a)')'           =MMMZ..   .=. ..         ..      =. ,:     .MMMMMMMMM.               '
write(66,'(a)')'           .MMMM~..  7  .MM.M   M   MM.    M.  ..  M  .MMMMMMMM+.               '
write(66,'(a)')'.............MMMM....Z. .MM.M...M...MM..O:.  M . ~.M...ZMMMMMMM.................'
write(66,'(a)')'              MMMM,. ., :  ,:   :  :  :    .,  ,,  ,  . MMMMMMM                 '
write(66,'(a)')'...............MMMM7. ..................................MMMMMM..................'
write(66,'(a)')'................MMMM8 ..................................MMMMMM..................'
write(66,'(a)')'.................NMMMM  ................................MMMMD. .................'
write(66,'(a)')'..................7MMMM ............................... MMMM  ..................'
write(66,'(a)')'....................MMMMO............................. :MMM.....................'
write(66,'(a)')'                    .MMMMM                          . . MMM                     '
write(66,'(a)')'.......................MMMMM..........................MMMMM.....................'
write(66,'(a)')'........................ZMMMMD.......................MMMMM......................'
write(66,'(a)')'..........................MMMMMM.................. MMMMM,.......................'
write(66,'(a)')'........................ ...MMMMMM. ............OMMMMMZ,........................'
write(66,'(a)')'........................... ..NMMMMM8.......=MMMMMMMI...........................'
write(66,'(a)')'.............................. .=MMMMMMMMMMMMMMMMM..............................'
write(66,'(a)')'.................................. NMMMMMMMMMM,.................................'
write(66,'(a)')'                                   .. . ....                                    '
write(66,'(a)')'................................................................................'
write(66,'(a)') ' '

return 
end

subroutine timestamp

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write(66,'(a)') ' '

  write ( 66, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  write(66,'(a)') ' '
  write(66,'(a)') ' '

  return
end
