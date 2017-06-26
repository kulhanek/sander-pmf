! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Run single point energy calculation
!-----------------------------------------------------------------------
!     --- RUNEXT ---
!-----------------------------------------------------------------------

subroutine runext(xx,ix,ih,ipairs,x,fg,w,ib,jb,conp, &
      winv,igrp,skips,ene,carrms, qsetup)

   use fastwt
   use constants, only : zero, one, TEN_TO_MINUS5, TEN_TO_MINUS6
   use qmmm_module, only : qmmm_nml, qmmm_struct, qmmm_mpi, qm2_struct
   use poisson_boltzmann, only: outwat, oution
   use bintraj, only: end_binary_frame
   use file_io_dat

#if defined( MPI )
   use evb_data, only: evb_frc
#endif /* MPI */

#ifdef RISMSANDER
   use sander_rism_interface, only: rismprm, rism_force
#endif

#ifdef PMFLIB
   use pmf_sander
   use nblist, only: a,b,c,alpha,beta,gamma
#endif
   
   use state
   use sebomd_module, only : sebomd_obj
   implicit none

#ifdef MPI
#  include "parallel.h"
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
   include 'mpif.h'
   integer ierr
#ifdef CRAY_PVP
#  define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
   integer ist(MPI_STATUS_SIZE), partner
#endif
#include "../include/md.h"
#include "box.h"
#include "../include/memory.h"
#include "nmr.h"
#include "extra.h"
#include "ew_cntrl.h"
#include "../pbsa/pb_md.h"
#include "def_time.h"

   ! ------ passed in variables --------------------
   _REAL_   xx(*)
   integer  ix(*), ipairs(*)
   character(len=4) ih(*)
   _REAL_   x(*),fg(*),w(*)
   integer  ib(*),jb(*)
   _REAL_   conp(*),winv(*)
   integer  igrp(*)
   logical  skips(*)
   type(state_rec) ::  ene
   _REAL_   carrms
   logical :: qsetup

   ! ------ External Functions -----------------
   _REAL_ ddot

   ! ------ local variables --------------------
   _REAL_ f
   _REAL_ rms,sum
   integer i,j,ier
   integer n_force_calls

#ifdef RISMSANDER
   _REAL_ erism
#endif /*RISMSANDER*/

#if defined PMFLIB
   _REAL_           :: pmfene
#endif

   ! kulhanek ? do_list_update - is list initialized?
   logical :: do_list_update=.false.

   !     ----- EVALUATE SOME CONSTANTS -----

   if (imin /= 5 .and. master) call amopen(7,mdinfo,'U','F','W')
   n_force_calls = 1

   ! kulhanek - SANITY CHECKS
   if ( icfe .ne. 0 )then
        call sander_bomb('runext','icfe .ne. 0','not implemented')
   end if

#ifdef MPI
   call sander_bomb('runext','not tested','it is not safe to call parallel implementation')
   if ( ifsc == 1 )then
      call sander_bomb('runext','ifsc == 1','not implemented')
   end if
#endif

   if(master)then
      if ( igb == 10 .or. ipb /= 0 ) then
         pbgrid = .true.
         pbprint = .true.
         ntnbr = 1
         ntnba = 1
         npbstep = n_force_calls
      end if
   endif

   if (sebomd_obj%do_sebomd) then
     ! write down atomic charges and density matrix if needed
     sebomd_obj%iflagch = 1
   endif
   
   ! calculate force
   iprint=1
   irespa = 0
   ! Zero the state type as done in runmd()
   ene = null_state_rec

   ! nrp - number of atoms, adjusted for LES copies
   do j = 1,nrp*3
      fg(j) = 0.0
   end do

   select case(ntmin)
        case(1)
            call force(xx,ix,ih,ipairs,x,fg,ene,ene%vir, &
                    xx(l96), xx(l97), xx(l98),xx(l99),qsetup, do_list_update,n_force_calls)
#ifdef PMFLIB
            if( ntb .ne. 0 ) then
                call pmf_sander_update_box(a,b,c,alpha,beta,gamma)
            end if
            call pmf_sander_rstforce(natom,x,fg,ene%pot%tot,pmfene)
            ene%pot%constraint = ene%pot%constraint + pmfene
            ene%pot%tot = ene%pot%tot + pmfene
#endif
        case(2)
            ene%pot%rism = 0.0
#ifdef RISMSANDER           
                if(rismprm%rism == 1) then
                    call timer_start(TIME_RISM)
                    call rism_force(x,fg,erism,irespa,imin)
                    ene%pot%rism = erism
                    call timer_stop(TIME_RISM)
                endif
                ene%pot%tot = ene%pot%rism
#endif
        case default
            call sander_bomb('runext','ntmin .ne. 1 .or. ntmin .ne. 2','not implemented')
   end select

   sum = ddot(nrp*3,fg,1,fg,1)
   rms = sqrt(sum/(nrp*3))

   !     ----- WRITE THE FINAL ENERGY AND GRADIENT -----
   ! Tomas Bouchal & Petr Kulhanek
   ! nrp - number of atoms, adjusted for LES copies

   if (master) then
        open(unit=1612,file='ene_forces_byTB.dat')
        write(1612,997) nrp

        write(1612,996) 'ENERGIES'
        write(1612,998) ene%pot%tot,'# tot'  
        write(1612,998) ene%pot%bond,'# bond'
        write(1612,998) ene%pot%angle,'# angle'
        write(1612,998) ene%pot%dihedral,'# dihedral'
        write(1612,998) ene%pot%elec,'# elec'
        write(1612,998) ene%pot%elec_14,'# elec_14'
        write(1612,998) ene%pot%vdw,'# vdw'
        write(1612,998) ene%pot%vdw_14,'# vdw_14'
        write(1612,998) ene%pot%hbond,'# hbond'
        write(1612,998) ene%pot%constraint,'# constraint'
        write(1612,998) ene%pot%scf,'# scf'
        write(1612,998) ene%pot%gb,'# gb'
        write(1612,998) ene%pot%pb,'# pb'
        write(1612,998) ene%pot%rism,'# rism'
        write(1612,998) ene%pot%surf,'# surf'
        write(1612,998) ene%pot%disp,'# disp'
        write(1612,998) ene%pot%ct,'# ct'

        write(1612,996) 'GRADIENTS'
        do j = 1,nrp*3,3
            write(1612,999) -fg(j),-fg(j+1),-fg(j+2)
        end do

        write(1612,996) 'COORDINATES'
        do j = 1,nrp*3,3
            write(1612,999) x(j),x(j+1),x(j+2)
        end do
        close(1612)
    end if
    996 format(A)
    997 format(I8)
    998 format(E23.16,1X,A)
    999 format(E23.16,1X,E23.16,1X,E23.16)

   !     ----- WRITE THE FINAL RESULTS -----
   call report_min_results( n_force_calls, rms, x, &
         fg, ene, ih(m04), xx, ix, ih )  ! ih(m04) = atom names
   carrms = rms

#ifdef MPI
   if( ievb /= 0 ) then
      call evb_dealloc 
   endif
#endif

   return

end subroutine runext

