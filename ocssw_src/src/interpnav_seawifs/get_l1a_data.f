      subroutine get_l1a_data( prod_ID, nlin, msec, pos, 
     1     smat, att, tilt, iret)

c  get_l1a_data( prod_ID, nlin, msec, orb_vec, sen_mat, att, tilt, iret )
c
c  Purpose: reads the scan line times and several navigation arrays from
c     an open SeaWiFS L1A file.
c
c  Calling Arguments:
c
c     Name         Type    I/O     Description
c     --------     ----    ---     -----------
c     prod_ID      I*4      I      HDF file ID
c     nlin         I*4      I      number of scan lines in the L1A file
c     msec(*)      I*4      O      size nlin array of scan line times
c     pos(3,*)     R*4      O      size 3 x nlin array of position vectors
c     smat(3,3,*)  R*4      O      size 3 x 3 x nlin array of sensor matrices
c     att(3,*)     R*4      O      size 3 x nlin array of attitude angles
c     tilt(*)      R*4      O      size nlin array of tilt angles
c     iret         I*4      O      return code
c                                    =0, success
c                                    =-1, failure    
c
c  By: F. S. Patt, SAIC GSC, 24 Sep 98
c
c  Notes:  
c
c  Modification History:
c

      real*4 pos(3,*), smat(3,3,*), att(3,*), tilt(*)
      integer*4 prod_ID, iret
      integer*4 msec(*)
      integer*4 or_id, se_id, at_id, ti_id, ms_id, ind
      integer*4 istart(3), istr(3), idims(3), nflags(8)
      integer sfn2index, sfselect, sfrdata, sfendacc
      logical first/.true./
      data istart/3*0/, istr/3*1/, idims/3*1/
      data msecday/86400000/

      iret = 0

c  Get sdid for scan time times
      ind = sfn2index(prod_ID, 'msec')
      if (ind.eq.-1) then
         iret = -1
         write(*,*) 'Error getting index for msec'
         return
      end if
      ms_id = sfselect(prod_ID, ind)
      if (ms_id.eq.-1) then
         iret = -1
         write(*,*) 'Error selecting msec'
         return
      end if

c  Get sdid for orbit vectors
      ind = sfn2index(prod_ID, 'orb_vec')
      if (ind.eq.-1) then
         iret = -1
         write(*,*) 'Error getting index for orb_vec'
         return
      end if
      or_id = sfselect(prod_ID, ind)
      if (or_id.eq.-1) then
         iret = -1
         write(*,*) 'Error selecting orb_vec'
         return
      end if

c  Get sdid for attitude angles
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

c  Get sdid for sensor matrices
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
      
c  Get sdid for tilt angles
      ind = sfn2index(prod_ID, 'tilt')
      if (ind.eq.-1) then
         iret = -1
         write(*,*) 'Error getting index for tilt'
         return
      end if
      ti_id = sfselect(prod_ID, ind)
      if (ti_id.eq.-1) then
         iret = -1
         write(*,*) 'Error selecting tilt'
         return
      end if
        

c  Read sensor matrices
      idims(1) = 3
      idims(2) = 3
      idims(3) = nlin
      iret = sfrdata(se_id, istart, istr, idims, smat)
      if (iret.eq.-1) then
        write(*,*) 'Error reading sen_mat'
        return
      end if
     
c  Read orbit vectors
      idims(2) = nlin
      iret = sfrdata(or_id, istart, istr, idims, pos)
      if (iret.eq.-1) then
        write(*,*) 'Error reading orb_vec'
        return
      end if

c  Read attitude angles
      iret = sfrdata(at_id, istart, istr, idims, att)
      if (iret.eq.-1) then
        write(*,*) 'Error reading att_ang'
        return
      end if

c  Read scan line times
      idims(1) = nlin
      iret = sfrdata(ms_id, istart, istr, idims, msec)
      if (iret.eq.-1) then
        write(*,*) 'Error reading msec'
        return
      end if

c  Read tilt angles
      idims(1) = nlin
      iret = sfrdata(ti_id, istart, istr, idims, tilt)
      if (iret.eq.-1) then
        write(*,*) 'Error reading tilt'
        return
      end if

c  Free sdid's
      iret = sfendacc(or_id)
      iret = sfendacc(se_id) 
      iret = sfendacc(at_id)
      iret = sfendacc(ti_id)
      iret = sfendacc(ms_id)

c  Check for day rollover during scene
      if ((msec(nlin)-msec(1)).lt.(-msecday/2)) then

         ind = 1
         dowhile ((msec(ind+1)-msec(ind)).gt.(-msecday/2))
            ind = ind + 1
         end do
         do i = ind+1,nlin
            msec(i) = msec(i) + msecday
         end do
      end if

      return
      end

