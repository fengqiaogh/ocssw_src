c ===================================================================
c Program MKTLM
c
c Generate TLM_GAC and TLM_LAC files from L0_GAC files
c
c Written By: BA Franz, GSC, 21 Feb 96
c
c ===================================================================

      program mktlm
c
      implicit none
c
      integer*2 FREC
c#ifndef __sun
c      parameter (FREC = 4)
c#else
      parameter (FREC = 1)
c#endif
      
      integer*4 MNFLEN          ! # Bytes in L0 minor frame
      parameter (MNFLEN = 21504)!
      integer*4 HDRLEN          ! # Bytes in L0 header
      parameter (HDRLEN = 512)  !
      integer*2 MNFLUN          ! Input L0 logical unit
      parameter (MNFLUN = 10)   !
      integer*2 GTELUN          ! Output TLM_GAC logical unit
      parameter (GTELUN = 11)   !
      integer*2 LTELUN          !
      parameter (LTELUN = 12)   !
      integer*2 IDBYTE          ! Byte # in mnf where ID appears
      parameter (IDBYTE = 7)    !
      integer*2 GACID           ! GAC mnf ID
      parameter (GACID  = 15)   !
      integer*2 LACID           ! LAC mnf ID
      parameter (LACID  = 0)    !
c
      byte  mnf( MNFLEN )       ! L0 minor frame record
      byte  hdr( HDRLEN )       ! L0 header
c
      character*80 mnffile      ! L0 filename
      character*80 gtefile      ! GAC TLM filename
      character*80 ltefile      ! LAC TLM filename
      character*4  type         ! GAC or LAC or BOTH
c
      integer*4 mnfcnt          ! Minor frame counter
      integer*4 gaccnt          ! GAC mnf counter
      integer*4 laccnt          ! LAC mnf counter
      integer*2 status          ! Status flag (0=OK)
c
      integer*2 getmnf          ! Function to buffer input
      integer*2 mkgactlm        ! Function to make GAC TLM records
      integer*2 mklactlm        ! Function to make LAC TLM records

      logical*1 wantgac
      logical*1 wantlac

c     !
c     ! Initialization
c     !
      data mnfcnt  / 0 /
      data gaccnt  / 0 /
      data laccnt  / 0 /
      data wantgac /.false./
      data wantlac /.false./

c     !
c     ! Open input L0 file
c     !
      write(*,'(a)') 'Enter input L0 file name: '
      read(5,'(a)') mnffile
      open(MNFLUN,file=mnffile,access='direct',
     .            recl=512/FREC,err=900)

c     !
c     ! Get Output Type
c     !
      write(*,'(a)') 'Output TLM type (GAC, LAC, BOTH): '
      read(5,'(a)') type
 
      if ( type .eq. 'GAC' .or. type .eq. 'BOTH') then
          wantgac = .true.
c         !
c         ! Open output GAC TLM file
c         !
          write(*,'(a)') 'Enter output GAC TLM file name: '
          read(5,'(a)') gtefile
          open(GTELUN,file=gtefile,access='direct',
     .               status='unknown',recl=4/FREC,err=901)
      endif


      if ( type .eq. 'LAC' .or. type .eq. 'BOTH') then
          wantlac = .true.
c         !
c         ! Open output LAC TLM file
c         !
          write(*,'(a)') 'Enter output LAC TLM file name: '
          read(5,'(a)') ltefile
          open(LTELUN,file=ltefile,access='direct',
     .                status='unknown',recl=4/FREC,err=902)
      endif

      if (.not. wantgac .and. .not. wantlac) then
          write(*,*) 'Error: no output will be generated for type=',type
          stop
      endif

c     !
c     ! Skip header
c     !
      read(MNFLUN, rec=1, err=900) hdr

c     !
c     ! Read first record
c     !
      status = getmnf( mnflun, mnfcnt, mnf, MNFLEN)
      mnfcnt = mnfcnt + 1

c     !
c     ! Read through end of file
c     !
      dowhile( status .eq. 0 )

c         !
c         ! Check mnf type
c         !
          if ( wantgac .and. mnf(IDBYTE) .eq. GACID ) then

c             !
c             ! Write GAC TLM record
c             !
              status = mkgactlm(GTELUN,gaccnt,mnf)
              gaccnt = gaccnt + 1

          elseif ( wantlac .and. mnf(IDBYTE) .ne. GACID ) then 

c             !
c             ! Write LAC TLM record
c             !
              status = mklactlm(LTELUN,laccnt,mnf)
              laccnt = laccnt + 1

          endif

c         !
c         ! Read next record
c         !
          status = getmnf( MNFLUN, mnfcnt, mnf, MNFLEN)
          mnfcnt = mnfcnt + 1

          if (mod(mnfcnt,1000) .eq. 0) then
              write(*,*) 'Processing minor frame ',mnfcnt
          endif

      enddo
c
      write(*,*) ' '
      write(*,*) 'Total minor frames: ',mnfcnt-1
      write(*,*) 'GAC frames:         ',gaccnt
      write(*,*) 'LAC frames:         ',laccnt
      write(*,*) ' '
c
      close( MNFLUN )
      close( GTELUN )
      close( LTELUN )

      stop
c
 900  write(*,*) 'Error opening input file: ',mnffile
      goto 1000
 901  write(*,*) 'Error opening output file: ',gtefile
      goto 1000
 902  write(*,*) 'Error opening output file: ',ltefile
      goto 1000
c
 1000 continue
      end


c ----------------------------------------------------------------------
c getmnf - buffers input of L0 minor frame records
c ----------------------------------------------------------------------
      integer*2 function getmnf(lun,mnfcnt,mnf,mnflen)
c
      implicit none
c
      integer*4 BUFLEN
      parameter (BUFLEN = 512)
c
      byte      mnf(*)
      integer*2 lun
      integer*4 mnfcnt
      integer*4 mnflen
c
      byte      buf( BUFLEN )
      integer*4 i, j
      integer*4 nbuf
c
c      write(*,*) "getmnf"
c
      nbuf = mnflen/BUFLEN
c
      do i=1,nbuf
          read(lun,rec=mnfcnt*nbuf+i+1,err=10) buf
          do j=1,BUFLEN
              mnf( (i-1)*BUFLEN + j ) = buf(j)
          enddo
      enddo
c
      getmnf = 0
      return
c
 10   getmnf = 1
      return
c
 20   getmnf = 2
      return
c
      end


c ----------------------------------------------------------------------
c puttlm - buffers output of GAC or LAC TLM records
c ----------------------------------------------------------------------
      integer*2 function puttlm(lun,cnt,tlm,tlmlen)
c
      implicit none
c
      integer*2 BUFLEN
      parameter (BUFLEN = 4)
c
      integer*2 lun
      integer*4 cnt
      byte      tlm(*)
      integer*4 tlmlen
c
      byte      buf(BUFLEN)
      integer*4 nbuf
      integer*4 i, j
c
c      write(*,*) "puttlm"
c
      nbuf = tlmlen/BUFLEN
c
      do i=1,nbuf
          do j=1,BUFLEN
              buf(j) = tlm( (i-1)*BUFLEN + j )
          enddo
          write(lun,rec=cnt*nbuf+i,err=10) buf
      enddo
c
      puttlm = 0
      return
c
 10   puttlm = 1
      return
c
      end


c ----------------------------------------------------------------------
c mkgactlm - writes one GAC TLM record
c
c Synopsis:
c       status = mkgactlm( lun, mnf )
c
c       integer*2 status : 0=success, 1=error, 2=end-of-file
c       integer*2 lun    : GAC TLM file unit number
c       integer*1 mnf(*) : L0 minor frame as byte array
c
c Written By:
c       BA Franz, GSC, 21 Feb 96
c
c ----------------------------------------------------------------------
      integer*2 function mkgactlm(lun,cnt,mnf)
c
      implicit none
c
      integer*4 TLMLEN                  ! # Bytes in TLM Record
      integer*2 NFLDS                   ! # Fields in TLM Record
      
      parameter (TLMLEN = 1308)
      parameter (NFLDS  = 9)
c
      integer*2 lun                     ! TLM file unit
      integer*4 cnt                     ! Number of logical records
      byte      mnf(*)                  ! L0 mnf record as byte array

      integer*4 nbytes                  ! # Bytes for this field
      integer*4 tlmoffset               ! # Byte offset of field in tlm 
      integer*4 mnfoffset               ! # Byte offset of field in mnf 
      integer*4 i, j
      byte      tlm( TLMLEN )           ! TLM record as byte array
      save      tlm
c
      logical   firstCall
      save      firstCall
c
      integer*2 puttlm                  ! # Function to buffer output

c     !
c     ! Define mapping from input L0 file to output TLM file
c     ! Format is MNF byte offset, TLM byte offset, # bytes
c     !
      integer*2 map( 3, NFLDS )
      data map /   4,    1,   4,        ! S/C ID
     .             8,    5,   8,        ! S/C Time-Tag
     .            16,   13, 775,        ! SOH TLM
     .         20951,  789, 440,        ! INST TLM
     .           791, 1229,  16,        ! Gain & TDI, Scan #1
     .          4823, 1245,  16,        ! Gain & TDI, Scan #2
     .          8855, 1261,  16,        ! Gain & TDI, Scan #3
     .         12887, 1277,  16,        ! Gain & TDI, Scan #4
     .         16919, 1293,  16 /       ! Gain & TDI, Scan #5

      data      firstCall /.true./

c
c      write(*,*) "mkgactlm"
c

c     !
c     ! If this is the first time through, initialize output record
c     !
      if ( firstCall ) then
          firstCall = .false.
          do i=1,TLMLEN
              tlm(i) = 0
          enddo
      endif



c     !
c     ! Copy MNF Fields to TLM Record
c     !
      do i=1, NFLDS
          mnfoffset = map(1,i)
          tlmoffset = map(2,i)
          nbytes    = map(3,i)
          do j=0,nbytes-1
              tlm( tlmoffset + j ) = mnf( mnfoffset + j )
          enddo 
      enddo

c     !
c     ! Write TLM record and return status
c     !
      mkgactlm = puttlm(lun,cnt,tlm,TLMLEN)

      return
c
      end


c ----------------------------------------------------------------------
c mklactlm - writes one LAC TLM record
c
c Synopsis:
c       status = mklactlm( lun, mnf )
c
c       integer*2 status : 0=success, 1=error, 2=end-of-file
c       integer*2 lun    : LAC TLM file unit number
c       integer*1 mnf(*) : L0 minor frame as byte array
c
c Written By:
c       BA Franz, GSC, 21 Feb 96
c
c ----------------------------------------------------------------------
      integer*2 function mklactlm(lun,cnt,mnf)
c
      implicit none
c
      integer*4 TLMLEN                  ! # Bytes in TLM Record
      integer*2 NFLDS                   ! # Fields in TLM Record
      
      parameter (TLMLEN = 892)
      parameter (NFLDS  = 5)
c
      integer*2 lun                     ! TLM file unit
      integer*4 cnt                     ! Number of logical records
      byte      mnf(*)                  ! L0 mnf record as byte array

      integer*4 nbytes                  ! # Bytes for this field
      integer*4 tlmoffset               ! # Byte offset of field in tlm 
      integer*4 mnfoffset               ! # Byte offset of field in mnf 
      integer*4 i, j
      byte      tlm( TLMLEN )           ! TLM record as byte array
      save      tlm
c
      logical   firstCall
      save      firstCall
c
      integer*2 puttlm                  ! # Function to buffer output

c     !
c     ! Define mapping from input L0 file to output TLM file
c     ! Format is MNF byte offset, TLM byte offset, # bytes
c     !
      integer*2 map( 3, NFLDS )
      data map /   4,    1,   4,        ! S/C ID
     .             8,    5,   8,        ! S/C Time-Tag
     .            16,   13, 775,        ! SOH TLM
     .           791,  789,  88,        ! INST TLM
     .           879,  877,  16/        ! Gain & TDI

      data      firstCall /.true./

c
c      write(*,*) "mklactlm"
c

c     !
c     ! If this is the first time through, initialize output record
c     !
      if ( firstCall ) then
          firstCall = .false.
          do i=1,TLMLEN
              tlm(i) = 0
          enddo
      endif


c     !
c     ! Copy MNF Fields to TLM Record
c     !
      do i=1, NFLDS
          mnfoffset = map(1,i)
          tlmoffset = map(2,i)
          nbytes    = map(3,i)
          do j=0,nbytes-1
              tlm( tlmoffset + j ) = mnf( mnfoffset + j )
          enddo 
      enddo

c     !
c     ! Write TLM record and return status
c     !
      mklactlm = puttlm(lun,cnt,tlm,TLMLEN)

      return
c
      end


