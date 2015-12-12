!-------------------------------------------------
!     --- SANDER_BOMB ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Print an error message and quit
subroutine sander_bomb(routine,string1,string2)
   implicit none
   character(len=*) routine,string1,string2

   write(6, '(1x,2a)') &
         'SANDER BOMB in subroutine ', routine
   write(6, '(1x,a)') string1
   write(6, '(1x,a)') string2
   call mexit(6,1)
end subroutine sander_bomb