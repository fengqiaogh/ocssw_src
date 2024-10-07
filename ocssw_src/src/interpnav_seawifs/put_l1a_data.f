      subroutine put_l1a_data( prod_ID, nlin, smat, att, coef, iret)

c  put_l1a_data( prod_ID, nlin, smat, att, coef, iret)
c
c  Purpose:  writes the attitude angles, sensor orientation matrix and scan
c     ellipse coefficients to a SeaWiFS L1A file 
c
c  Calling Arguments:
c
c     Name         Type    I/O     Description
c     --------     ----    ---     -----------
c     prod_ID      I*4      I      HDF file ID for open file 
c     nlin         I*4      I      number of scan lines in the L1A file
c     smat(3,3,*)  R*4      I      size 3 x 3 x nlin array of sensor matrices
c     att(3,*)     R*4      I      size 3 x nlin array of attitude angles
c     coef(6,*)    R*4      I      size 3 x nlin array of ellipse coefficients
c     iret         I*4      O      return code
c                                    =0, success
c                                    =-1, failure    
c
c  By: F. S. Patt, SAIC GSC, 25 Sep 98
c
c  Notes:  
c
c  Modification History:
c
c  Modified to reset first nflag value to 0 in each scan line.
c  F. S. Patt, SAIC GSC, 2 Feb 1999.
c
c  Modified to reset both navigation fail and warning flags
c  F. S. Patt, SAIC, March 19, 2003

#include "nav_cnst.fin"

      real*4 smat(3,3,*), att(3,*), coef(6,*)
      integer*4 prod_ID, iret
      integer*4 se_id, at_id, co_id, nf_id, ind, i
      integer*4 istart(3), istr(3), idims(3)
      integer*4 nflag(8,maxlin)
      integer sfn2index, sfselect, sfwdata, sfendacc
      data istart/3*0/, istr/3*1/

      iret = 0

      ind = sfn2index(prod_ID, 'att_ang')
      if (ind.eq.-1) then
         iret = -1
         write(*,*) 'Error getting index for att_ang'
         return
      end if
      at_id = sfselect(prod_ID, ind)
      if (at_id.eq.-1) then
         iret = -1
         write(*,*) 'Error selecting att_ang'
         return
      end if
      
      ind = sfn2index(prod_ID, 'sen_mat')
      if (ind.eq.-1) then
         iret = -1
         write(*,*) 'Error getting index for sen_mat'
         return
      end if
      se_id = sfselect(prod_ID, ind)
      if (se_id.eq.-1) then
         iret = -1
         write(*,*) 'Error selecting sen_mat'
         return
      end if
      
      ind = sfn2index(prod_ID, 'scan_ell')
      if (ind.eq.-1) then
         iret = -1
         write(*,*) 'Error getting index for scan_ell'
         return
      end if
      co_id = sfselect(prod_ID, ind)
      if (co_id.eq.-1) then
         iret = -1
         write(*,*) 'Error selecting scan_ell'
         return
      end if
        
      ind = sfn2index(prod_ID, 'nflag')
      if (ind.eq.-1) then
         iret = -1
         write(*,*) 'Error getting index for nflag'
         return
      end if
      nf_id = sfselect(prod_ID, ind)
      if (nf_id.eq.-1) then
         iret = -1
         write(*,*) 'Error selecting nflag'
         return
      end if
        
      idims(1) = 3
      idims(2) = 3
      idims(3) = nlin
      iret = sfwdata(se_id, istart, istr, idims, smat)
      if (iret.eq.-1) then
        write(*,*) 'Error writing sen_mat'
        return
      end if
     
      idims(2) = nlin
      iret = sfwdata(at_id, istart, istr, idims, att)
      if (iret.eq.-1) then
        write(*,*) 'Error writing att_ang'
        return
      end if

      idims(1) = 6
      iret = sfwdata(co_id, istart, istr, idims, coef)
      if (iret.eq.-1) then
        write(*,*) 'Error writing scan_ell'
        return
      end if

      idims(1) = 8
      iret = sfrdata(nf_id, istart, istr, idims, nflag)
      if (iret.eq.-1) then
        write(*,*) 'Error reading nflag'
        return
      end if

c    Reset navfail and navwarn flags 
      do i=1,nlin
         nflag(1,i) = 0
         if (nflag(7,i).eq.0) nflag(8,i) = 0
      end do

      iret = sfwdata(nf_id, istart, istr, idims, nflag)
      if (iret.eq.-1) then
        write(*,*) 'Error writing nflag'
        return
      end if

      iret = sfendacc(se_id) 
      iret = sfendacc(at_id)
      iret = sfendacc(co_id)
      iret = sfendacc(nf_id)

      return
      end

