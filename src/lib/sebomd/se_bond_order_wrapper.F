      subroutine se_bond_order_wrapper
! F77 wrapper for se_bond_order_write
      implicit none
#include "sebomd.dim"
#include "se_fulsys.h"  /* ! natoms, iatnum */
#include "se_elemnt.h"  /* ! natorb         */
#include "se_bpairs.h"  /* ! ip1            */
#include "se_glbmat1.h" /* ! pdiat          */
#include "se_glbmat2.h" /* ! ijmat          */
      double precision bocut
      call se_bond_order_write(natoms,
     &                         iatnum, 
     &                         natorb(1), 
     &                         ip1, 
     &                         ijmat, 
     &                         pdiat)
      return
      end
