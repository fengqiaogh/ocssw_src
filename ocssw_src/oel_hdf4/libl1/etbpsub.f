C
C@***s* PROJECT_PFSST/l2gen_both/etbpsub.f
C
C This header contains documentation required by the NOAA Climate Data Record
C Program (CDRP), which is managed at the NOAA National Climatic Data Center (NCDC).
C Only the code that applies to AVHRR SST data is documented in this header. 
C
C The AVHRR Pathfinder Sea Surface Temperature (PFSST) processing code was originally 
C developed at the University of Miami.  In 2010, the code was integrated into
C the multi-sensor SeaWiFS Data Analysis System (SeaDAS) obtained from NASA GSFC.
C SeaDAS was used for processing the PFSST beginning with the Pathfinder Version
C 5.2 (PFV5.2) dataset, produced jointly by the University of Miami and the NOAA
C National Oceanographic Data Center (NODC).  These data are provided to the 
C public and archived at NODC, and have been transitioned along with the production 
C code and documentation to the CDRP at NCDC.
C 
C This NOAA required header is specifically written for Pathfinder SST and may
C not be relevant to other sensors or products processed by SeaDAS. Please 
C review the SEADAS software distribution policy for public domain software 
C located at http://seadas.gsfc.nasa.gov/copying.html for more information and 
C documentation
C
C NAME
C      etbsub.f
C
C LOCATION
C      $OCSSWROOT 
C
C PURPOSE
C       etbsub.f contains a collection of functions and subroutines
C         relating to emissivity and temperature conversions, and calibration
C         of both infrared and visible AVHRR channels.
C
C DESCRIPTION
C
C       ETINVERT(ICH,RADI) - Convert radiance to temperature
C       ETINTEGRATE(ICH,ETEMP) - Convert temperature to radiance
C       ETLOADRESP(LIN,CAL) - Initialize temperature conversion routines,
C            load and normalize sensor response factors
C       ETGETRSP(LIN,ICHN,NPTS,RESP,NU0,DELNU,IOS)-Get infrared channel response
C       ETGETVIS(LIN,CAL,SLOPE,INTCP,IOS)- Get visible channel calibration
C
C       see additional details in file
C
C
C
C NOAA PFSST-SEADAS BUILD VERSION
C       Pathfinder SST V5.2 code built with SEADAS version 6.3 64 bit l2gen_both for  
C       CDR processed at University of Miami/RSMAS
C       http://seadas.gsfc.nasa.gov/seadas/doc/l2gen/l2gen.html
C
C PRIMARY SEADAS CODE DOCUMENTATION NASA
C       For complete documentation of multi sensor SEADAS code see 
C       http://seadas.gsfc.nasa.gov/seadas/doc/toplevel/sds_program.html
C
C AUTHOR
C                                      
C                              
C   PFSST project embedded code 
C     Susan Walsh
C     University of Miami/RSMAS
C
C CREATION DATE 
C  2010  
C
C COPYRIGHT
C  THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC DOMAIN AND
C  THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE.  THEY ARE FURNISHED "AS IS." THE 
C  AUTHORS, THE UNITED STATES GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, 
C  AND AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE 
C  SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
C  THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT 
C  TO USERS.
C
C MODIFICATION HISTORY
C 
C  See CVS revision control history embedded in actual file.
C               
C
C LANGUAGE
C Fortran   
C
C*****
C########################################################
C
C !F90
C
C !Description:
C    etbsub routine contains a collection of functions and subroutines
C    relating to emissivity and temperature conversions, and calibration
C    of both infrared and visible channels.
C
C Subroutines and Functions:
C    function ETINVERT(ICH,RADI)
C    !Description:
C      Convert radiance to temperature
C    !Input Parameters:
C      ICHI(integer) - infrared channel number
C      RADI(single) - radiance
C    !Output Parameters:
C      None
C    !Return(single)- temperature
C
C    function ETINTEGRATE(ICH,ETEMP)
C    !Description:
C      Convert temperature to radiance
C    !Input Parameters:
C      ICH(integer) - Infrared channel number
C      ETEMP(single) - emissivity temperature
C    !Output Parameters:
C      None
C    !Return(single) - radiance
C
C    subroutine ETLOADRESP(LIN,CAL)
C    !Description:
C      Initialize temperature conversion routines,
C       load and normalize sensor response factors
C    !Input Parameters:
C      LIN(integer) - logical unit of input file
C      CAL(char,12) - name satellite specific sensor calibration file
C    !Output Parameters:
C
C    subroutine ETGETRSP(LIN,ICHN,NPTS,RESP,NU0,DELNU,IOS)
C    !Description:
C      Get infrared channel response
C    !Input Parameters:
C      LIN(integer) - logical unit for file access
C      ICHN(integer) - infrared channel
C    !Output Parameters:
C      NPTS(integer) - number of points in the response vector
C      RESP(single array size, MAX_RSP) - response factor for channel
C      NU0(single) - base wave number
C      DELNU(double) - delta wavenumber
C      IOS(integer) - status error code
C
C
C    subroutine ETGETVIS(LIN,CAL,SLOPE,INTCP,IOS)
C    !Description:
C      Get visible channel calibration
C    !Input Parameters:
C      LIN(integer) - logical unit for input file
C       CAL(char,12) - name of satellite specific sensor calibration file
C    !Output Parameters:
C      SLOPE(single array size, 2) - calibration slope
C      INTCP(single array size, 2) - calibration intercept
C      IOS(integer) - error code
C !Revision History:
C
C
C !Team-Unique Header:
C
C Copyright 1988-2000 by Rosenstiel School of Marine and Atmospheric Science,
C University of Miami, Miami, Florida.
C
C                       All Rights Reserved
C
C Permission to use, copy, modify, and distribute this software and its
C documentation for non-commercial purposes and without fee is hereby granted,
C provided that the above copyright notice appear in all copies and that both
C that copyright notice and this permission notice appear in supporting
C documentation, and that the names of University of Miami and/or RSMAS not be
C used in advertising or publicity pertaining to distribution of the software
C without specific, written prior permission.
C
C UNIVERSITY OF MIAMI DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
C INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT

C SHALL UNIVERSITY OF MIAMI BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
C DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
C WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING
C OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
C
C !References and Credits:
C    Written BY:
C    University Of Miami
C    Rosenstiel School for Marine and Atmospheric Science
C    Division of Meteorology and Physical Oceanography
C    4600 Rickenbacker CSWY
C    Miami, Florida
C    33149
C    Contact:SWalsh@rsmas.miami.edu
C
C !Design Notes:
C
C !End###################################################################

        function etinvert(ich,radi)
        implicit none

        integer*4   ich;
        real*4      etinvert;
        real*4      radi;
       
        real*4      gues;
        real*4      etintegrate;
        real*4      step;
        real*4      temp;
        real*4      tnew, radg;

        if (ich < 3 .or. ich > 5) then
            write(*,*)
     .         '-E- : ETINVERT: ICH must be [3,5].  Is ', ich
            call exit(1)
        end if

        if (radi <= 0.) then
            etinvert=(100.)            ! Bad input value
            return
        end if

        step = 10.;
        temp = 280.;
        radg = etintegrate(ich,temp);
        if (radg > radi) then
            step = -step;
        end if

        do while (abs(radi-radg)/abs(radi+radg) >= 1e-5)
            gues = etintegrate(ich,temp+step);
            if (radg >= radi .and. gues <= radi) then
                tnew = temp+step*abs(radi-radg)/abs(radg-gues)
                temp = tnew;
                radg = etintegrate(ich,temp);
                step = step/10.;

            else if (radg <= radi .and. gues >= radi) then
                tnew = temp+step*abs(radi-radg)/abs(radg-gues)
                temp = tnew;
                radg = etintegrate(ich,temp);
                step = step/10.;

            else
                temp = temp+step;
                radg = gues;
            end if

            if (temp < 100.) then
C                write(*, *)'ETINVERT: ',RADI,RADG,TEMP,STEP
                etinvert=(100.);            ! Too small, quit now
                return

            else if (temp > 375.) then
C                write(*, *)'ETINVERT: ',RADI,RADG,TEMP,STEP
                etinvert=(375.);            ! Too large, quit now
                return
            end if

            if (radg < radi) then
                step = abs(step);
            else
                step = -abs(step);
            end if

        end do

        etinvert=(temp);
        return

        end


        function etintegrate(ich,etemp)
        implicit none


C      size of sensor response table
        integer, parameter :: MAX_RSP = 250

        integer*4   ich;
        real*4      etintegrate;
        real*4      etemp;
       
        integer*4   npts(3:5)
        integer*4   jj, ierr;
        real*4      nu0(3:5);
        real*8      delnu(3:5);
        real*4      resp(max_rsp,3:5);
        common /resp/ npts, nu0, delnu, resp;

        real*4      temp;
        real*4      y(max_rsp);
        real*4      simpsn, energy;
        real*4      v;
        real*8      c1,c2,nu,lam,t,w
        real*8      tt,nnu
       
C       data      C1/3.74125d-9/
C       data      C2/14.380d-4/

        data      c1 /1.1910659d-5/
        data      c2 /1.438833d0/

C       The function
C        (units are: milliWatts per sq meter sterradian inverse-cm)

C      Large T
        w(tt,nnu) = (c1*nnu*nnu*nnu)/((exp(c2*nnu/tt)-1.d0))
C      Small T
        v(tt,nnu)=sngl((c1*nnu*nnu*nnu*exp(-c2*nnu/tt))/
     .                 (1.d0-(exp(-c2*nnu/tt))))

        if (ich < 3 .or. ich > 5) then
            write(*,*)
     .        '-E- : ETINTEGRATE: ICH must be [3,5].  Is ', ich
            call exit(1)
        end if

        if (npts(ich) > max_rsp) then
            write(*,*) '-E- : ETINTEGRATE: NPTS > MAX_RSP'
            call exit(1)

        else if (npts(ich) <= 0) then
            write(*,*) '-E- : ETINTEGEMIS: NPTS <= 0'
            call exit(1)
        end if

        temp = etemp;
        if (temp < 100.) then
C          write(*,*)'ETINTEGRATE: ',ICH,TEMP
            temp = 100.;

        else if (temp > 375.) then
C          write(*,*) 'ETINTEGRATE: ',ICH,TEMP
            temp = 375.;
        end if

        if (nu0(ich) > 100.) then
            nu = dble(nu0(ich));            ! From left to right
            if (c2*nu/dble(temp) > 20.0d0) then
                do jj=1,npts(ich)

                    y(jj) = v(dble(temp),nu) * resp(jj,ich);
                    nu = nu + delnu(ich);
                end do

            else
                do jj=1,npts(ich)

                    y(jj) = sngl(w(dble(temp),nu)) * resp(jj,ich);
                    nu = nu + delnu(ich);
                end do

            end if

        else
            lam = dble(nu0(ich));            ! From left to right
            nu = 10000.0d0/lam;
            if (c2*nu/dble(temp) > 20.0d0) then
                do jj=1,npts(ich)

                    nu = 10000.0d0/lam;
                    y(jj) = v(dble(temp),nu) * resp(jj,ich);
                    lam = lam + delnu(ich);
                end do

            else
                do jj=1,npts(ich)

                    nu = 10000.0d0/lam;
                    y(jj) = sngl(w(dble(temp),nu)) * resp(jj,ich);
                    lam = lam + delnu(ich);
                end do

            end if

        end if

C      Integrate over this interval [NU0,NUN] or [LAM0,LAMN]

        energy = simpsn(sngl(delnu(ich)),y,npts(ich),ierr);
        if (ierr .ne. 0) then
            write(*,*) '-E- : SIMPSN: Integration error..'
            call exit(1)
        end if

        etintegrate=(energy);
        return

        end

        subroutine etloadresp(lin, cal)
        implicit none

C      size of sensor response table
        integer, parameter :: MAX_RSP = 250

        integer*4   lin
        character   cal*12;

        integer*4   npts(3:5)
        real*4      nu0(3:5);
        real*8      delnu(3:5);
        real*4      resp(max_rsp,3:5);
        common /resp/ npts, nu0, delnu, resp;

        logical     loaded;
        character   filepath*255;
        integer*4   ii, jj, ios, ier;
        integer*4   lenstr;
        integer*4   trname
        integer     slen, flen
        real*4      factor;
        real*4      simpsn;
        character filedir*255

        data      loaded /.false./;

       !  Read in response function

        if (loaded) then
            return

        end if

        call getenv('OCDATAROOT', filedir)
        flen = lenstr(filedir)
        if (flen .eq. 0) then
            write(*,*)
     .         '-E- ETLOADRESP:Environment variable OCDATAROOT undefined'
            call exit(1)
        end if

        slen = lenstr(cal)
        filepath = filedir(1:flen)//'/avhrr/cal/'//cal(1:slen)//'.cal'
        write(*,*) 'Loading '//filepath(1:lenstr(filepath))
        open(unit=lin,file=filepath, status='OLD');
        do ii=3,5

            call etgetrsp(lin,ii,npts(ii),resp(1,ii),nu0(ii),
     .                       delnu(ii),ios);
            if (ios .ne. 0) then
                call exit(1);
            end if

            !  Calculate normalizing factor re-scale response function

            factor = simpsn(sngl(delnu(ii)),resp(1,ii),npts(ii),ier);
            if (ier .ne. 0) then
                write(*,*) '-E- : SIMPSN: Normalizing factor.'
                call exit(1);
            end if

            !  Re-scale response function

            do jj=1,npts(ii)

                resp(jj,ii) = resp(jj,ii) / factor;
            end do

        end do

        close(unit=lin);

        loaded = .true.;

        end

        subroutine etgetrsp(lin,ichn,npts,resp,nu0,delnu,ios)
        implicit none


C size of sensor response table
        integer, parameter :: MAX_RSP = 250

        integer*4   lin
        integer*4   ichn
        integer*4   npts
        integer*4   ios
        real*4      resp(max_rsp)
        real*4      nu0
        real*8      delnu

        character*(128) name
        integer*4   iostat;
        integer*4   ftrim;
        integer*4   ln, ic, i;

C  Find appropriate channel

        rewind lin;
        do while (.TRUE.)
            read(lin, '(a)', end=1044) name
            ln = ftrim(name, 128);
            if (ln .eq. 0) then
                cycle

            else if (name(1:3) == "CHN") then
                read(name(1:5), 44, iostat=iostat) ic
44              format(3x,i2)
                if (ic .eq. ichn) then
C#                write(0, *)'Found response data for channel ',ic
                    exit

                end if

            end if

        end do

C  Input NU range

C      write(0, 51) (NAME(I:I),I=1,LN)
C      51 format('Decoding: ',80a1)
        read(name(1:ln), 45, iostat=iostat) nu0, delnu, npts
        if (npts > max_rsp) then
            write(*,*)
     .         '-E- : ETGETRSP: Too many points required (',npts,
     .         '), maximum is ', max_rsp
            call exit(1)
        end if

C  Input response function

        read(lin, '(x)')                   ! Skip a line
        read(lin,47) (resp(i),i=1,npts)
        ios = 0;
        return

C  No such channel

1044    continue
        write(*,*) '-E- : ETGETRSP: No such channel = ',ichn
        ios = -26
        return

45      format(11x,f14.0,4x,f10.0,5x,i3)
47      format(5(e11.5,4x))

        end


        subroutine etgetvis(lin,cal,slope,intcp,ios)
        implicit none

        integer*4      lin;
        character      cal*12;
        real*4      slope(2), intcp(2);
        integer*4      ios;

        character*(128) name;
        character      filepath*255;
        character      filedir*255;
        integer*4      ic;
        integer*4      ftrim;
        integer*4      lenstr;
        integer*4      trname
        integer*4      ln;
        integer            slen, flen;

C  Find appropriate channel

        call getenv('OCDATAROOT', filedir)
        flen = lenstr(filedir)
        if (flen .eq. 0) then
            write(*,*)
     .      '-E- : ETGETVIS: Environment variable OCDATAROOT undefined'
            call exit(1)
        end if

        slen = lenstr(cal)
        filepath = filedir(1:flen)//'/avhrr/cal/'//cal(1:slen)//'.cal'
        write(*,*) 'Loading '//filepath(1:lenstr(filepath))
        open(unit=lin,file=filepath, status='OLD');
        do while (.TRUE.)
            read(lin, '(a)',end=1044) name
        ln = ftrim(name, 128);
        if (ln .eq. 0) then
            cycle;
        else if (name(1:8) == "VISIBLE") then
C#            write(0, *)'Found visible calibration data';
            exit;
        end if

        end do

C  Input space radiance data

C      write(0, 51) (NAME(I:I),I=1,LN)
C      51 format('Decoding: ',80a1)
        ios = 0;
        read(lin,42);
        read(lin,42);
        read(lin,42);
42      format()
        read(lin,43) ic,slope(1),intcp(1);
        if (ic .ne. 1) then
        write(*,*)
     .     '-E- : ETGETVIS: Invalid visible entry for channel 1'
        ios = -1;
        go to 1043;
        end if

        read(lin,43) ic,slope(2),intcp(2);
43      format(i3,2f16.0)
        if (ic .ne. 2) then
        write(*,*)
     .     '-E- : ETGETVIS: Invalid visible entry for channel 2'
        ios = -2;
        go to 1043;
        end if

1043    continue
        close(unit=lin);
        return

C  No such channel

1044    continue
        write(*,*) '"Visible calibration data not found!'
        ios = -26;
        close(unit=lin);
        return

        end
