c  This program reads a 22180-byte SeaWiFS formatted record and compresses
c  the S/C ID, time tag, instrument telemetry and scan line data from 2-byte
c  words to 10 bits.  The data are written in 13860-byte records formatted
c  as SeaWiFS minor frames.  The leading and trailing frame synch words are
c  filled with zeros.  

        integer*2 ints(4),isync
        integer*4 irec
        byte inbuf(512),inrec(22180),outrec(13860),b10(5),bint(8)
        byte frsync(7),ausync(125)
        character*50 l0file
        character*50 outfile
        data outrec/13860*0/,ints/4*0/
        data isync/5/
        data frsync/Z'a1',Z'16',Z'fd',Z'71',Z'9d',Z'83',Z'c9'/
        data ausync/
     *   Z'f8',Z'bf',Z'36',Z'd6',Z'bd',Z'a1',Z'11',Z'57',Z'11',Z'd3',
     *   Z'd0',Z'4e',Z'0a',Z'db',Z'de',Z'37',Z'19',Z'9f',
     *   Z'c9',Z'93',Z'a3',Z'86',Z'90',Z'fb',Z'63',
     *   Z'12',Z'c9',Z'55',Z'02',Z'd5',Z'a7',Z'24',Z'5d',
     *   Z'88',Z'6d',Z'29',Z'ba',Z'df',Z'f8',Z'3b',
     *   Z'f7',Z'4b',Z'67',Z'34',Z'c5',Z'bb',Z'd6',Z'7b',
     *   Z'00',Z'42',Z'60',Z'ce',Z'ed',Z'4a',Z'ea',
     *   Z'76',Z'63',Z'd4',Z'68',Z'06',Z'35',Z'0a',Z'99',
     *   Z'be',Z'f9',Z'f4',Z'd5',Z'23',Z'e5',Z'c0',
     *   Z'52',Z'f8',Z'fd',Z'56',Z'18',Z'50',Z'eb',Z'fb',
     *   Z'21',Z'72',Z'07',Z'b8',Z'48',Z'3f',Z'd1',
     *   Z'47',Z'89',Z'e0',Z'6b',Z'1c',Z'b0',Z'46',Z'46',
     *   Z'c2',Z'03',Z'9e',Z'44',Z'd1',Z'05',Z'e9',
     *   Z'2e',Z'86',Z'56',Z'5a',Z'30',Z'25',Z'16',Z'6b',
     *   Z'98',Z'71',Z'db',Z'9c',Z'57',Z'd7',Z'72',
     *   Z'83',Z'79',Z'd5',Z'e5',Z'44',Z'93',Z'65',Z'27',Z'c3',Z'cc'/

        equivalence (ints,bint)

        write(*,*) 'Enter Level 0 file name'
        read (5,'(a50)') l0file
        write(*,*) 'Enter Output file name'
        read (5,'(a50)') outfile

        open(unit=11,file=l0file,access='direct',recl=22180,
     *    status='old',ACTION='READ')
        open(unit=12,file=outfile,access='direct',recl=13860)

c  Initialize output record with frame and auxillary synchs 
        do i=1,7
          outrec(i) = frsync(i)
        end do
        do i=1,125
          outrec(i+13735) = ausync(i)
        end do

c  Skip header records
        nrec = 0
        irec = 2

c  Now loop through input records
      dowhile (.true.)

c  Read a complete L0 record
        read (11,err=999,rec=irec) inrec
        irec = irec + 1
        nrec = nrec + 1

c  Reverse bytes (not needed for linux)
c        do i=1,22180,2
c            bval       = inrec(i)
c            inrec(i)   = inrec(i+1)
c            inrec(i+1) = bval
c        enddo
        
c  Remaining 10-bit data (instrument telemetry and scan line) 
        do l=1,2772
          do i=1,8
            bint(i) = inrec(i+(l-1)*8)
          end do
          call i2to10bit(ints,b10)
          do i=1,5
            outrec(i+(l-1)*5) = b10(i)
          end do
        end do

        write (12,rec=nrec) outrec

        if (mod(nrec,100) .eq. 0) then
            write(*,*) 'Writing record ',nrec
c            goto 999
        endif

      end do

  999 continue

      stop
      end

        subroutine i2to10bit(ints,bit10)

c  This subroutine compresses a 4-element I*2 array into a 5-element byte 
c  array, taking the 10 least significant bits of each integer.  

        integer*2 ints(4)
        byte bit10(5)
        integer*2 m3fc,m003,m3f0,m3c0,m300,m00f,m03f,m0ff
        data m3fc/1020/,m003/3/,m3f0/1008/,m00f/15/
        data m3c0/960/,m03f/63/,m300/768/,m0ff/255/
        bit10(1) = iand(ints(1),m3fc)/4
        bit10(2) = iand(ints(1),m003)*64 + iand(ints(2),m3f0)/16
        bit10(3) = iand(ints(2),m00f)*16 + iand(ints(3),m3c0)/64
        bit10(4) = iand(ints(3),m03f)*4  + iand(ints(4),m300)/256
        bit10(5) = iand(ints(4),m0ff)
        return
        end
