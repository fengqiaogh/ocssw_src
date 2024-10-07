      subroutine get_l1a_record( prod_ID, npix, irec, dat, nav_orb_vec, 
     1     nav_sun_ref, nav_sen_mat, nav_scan_ell, iret)

      real*4 nav_orb_vec(3), nav_sun_ref(3), nav_sen_mat(3,3), 
     1     nav_scan_ell(6)
      integer*4 prod_ID
      integer*2 dat(8,1285)
      integer*4 l1_id, or_id, su_id, se_id, sc_id, nf_id, ind
      integer*4 istart(3), istr(3), idims(3), nflags(8)
      integer sfn2index, sfselect, sfrdata
      logical first/.true./
      data istart/3*0/, istr/3*1/, idims/3*1/

      common /idcom/l1_id, or_id, su_id, se_id, sc_id, nf_id


      if (irec.eq.0) then

        ind = sfn2index(prod_ID, 'l1a_data')
        if (ind.eq.-1) then
          iret = -1
          write(*,*) 'Error getting index for l1a_data'
          return
        end if
        l1_id = sfselect(prod_ID, ind)
        if (l1_id.eq.-1) then
          iret = -1
          write(*,*) 'Error selecting l1a_data'
          return
        end if

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

        ind = sfn2index(prod_ID, 'sun_ref')
        if (ind.eq.-1) then
          iret = -1
          write(*,*) 'Error getting index for sun_ref'
          return
        end if
        su_id = sfselect(prod_ID, ind)
        if (su_id.eq.-1) then
          iret = -1
          write(*,*) 'Error selecting sun_ref'
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
        sc_id = sfselect(prod_ID, ind)
        if (sc_id.eq.-1) then
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
        if (sc_id.eq.-1) then
          iret = -1
          write(*,*) 'Error selecting nflag'
          return
        end if

      end if

      istart(2) = irec
      idims(1) = 8
      iret = sfrdata(nf_id, istart, istr, idims, nflags)
      if (iret.eq.-1) then
        write(*,*) 'End of l1a_data'
        return
      end if

      dowhile (nflags(1).ne.0)
        irec = irec + 1
        istart(2) = irec
        iret = sfrdata(nf_id, istart, istr, idims, nflags)
        if (iret.eq.-1) then
          write(*,*) 'End of l1a_data'
          return
        end if
      end do
 
      istart(2) = 0
      istart(3) = irec
      idims(1) = 8
      idims(2) = npix           
      iret = sfrdata(l1_id, istart, istr, idims, dat)
      if (iret.eq.-1) then
        write(*,*) 'End of l1a_data'
        return
      end if
 
      idims(1) = 3
      idims(2) = 3
      iret = sfrdata(se_id, istart, istr, idims, nav_sen_mat)
      if (iret.eq.-1) then
        write(*,*) 'Error reading sen_mat'
        return
      end if
     
      istart(2) = irec
      idims(2) = 1
      iret = sfrdata(or_id, istart, istr, idims, nav_orb_vec)
      if (iret.eq.-1) then
        write(*,*) 'Error reading orb_vec'
        return
      end if

      iret = sfrdata(su_id, istart, istr, idims, nav_sun_ref)
      if (iret.eq.-1) then
        write(*,*) 'Error reading sun_ref'
        return
      end if

      idims(1) = 6
      iret = sfrdata(sc_id, istart, istr, idims, nav_scan_ell)
      if (iret.eq.-1) then
        write(*,*) 'Error reading scan_ell'
        return
      end if

      irec = irec + 1

      return
      end

