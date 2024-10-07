      subroutine ssssss
c**********************************************************************c
c                                                                      c
c                                                                      c
c       ********************************************************       c
c       *           second simulation of satellite signal      *       c
c       *                 in the solar spectrum                *       c
c       *           ... (6s) ....... (6s) ...... (6s) ...      *       c
c       *                 .... (6s) ...... (6s)...             *       c
c       *                     ...... (6s) ......               *       c
c       *                        version  4.1                  *       c
c       *                                                      *       c
c       *  this code predicts the satellite signal from 0.25   *       c
c       *  to 4.0 microns assuming cloudless atmosphere.       *       c
c       *  the main atmospheric effects (gaseous absorption    *       c
c       *  by water vapor,carbon dioxyde,oxygen and ozone;     *       c
c       *  scattering by molecules and aerosols) are taken     *       c
c       *  into account. non-uniform surfaces may be           *       c
c       *  considered,as well as bidirectional reflectances    *       c
c       *            as boundary conditions                    *       c
c       *                                                      *       c
c       *   the following input parameters are needed          *       c
c       *         geometrical conditions                       *       c
c       *         atmospheric model for gaseous components     *       c
c       *         aerosol model (type and concentration)       *       c
c       *         spectral condition                           *       c
c       *         ground reflectance (type and spectral var.)  *       c
c       *   at each step, you can either select some proposed  *       c
c       *  standard conditions (for example,spectral bands of  *       c
c       *  satellite for spectral conditions) or define your   *       c
c       *  own conditions(in the example,you have to define    *       c
c       *  the assumed spectral response).                     *       c
c       *                                                      *       c
c       *   more details are given at each data input step     *       c
c       *                                                      *       c
c       ********************************************************       c
c                                                                      c
c                                                                      c
c**********************************************************************c
 
c**********************************************************************c
c                                                                      c
c                                                                      c
c       ********************************************************       c
c       *             authors of this code are                 *       c
c       *                                                      *       c
c       *  (1) Vermote E.; (2) Tanre D.;(2) Deuze J.L.         *       c
c       *           (2) Herman M.,(3) MOrcrette J.J..          *       c
c       *                                                      *       c
c       *                       from                           *       c
c       *                                                      *       c
c       *     (1) Affiliation: Department of Geography         *       c
c       *         University of Maryland                       *       c
c       *         address: Goddard Space Flight Center	       *       c
c       *         Code 923 		      		       *       c
c       *         greenbelt, md 20771                          *       c
c       *         USA                                          *       c
c       *                                                      *       c
c       *     (2) laboratoire d' optique atmospherique         *       c
c       *         universite des sciences et techniques        *       c
c       *         de lille                                     *       c
c       *         u.e.r. de physique fondamentale              *       c
c       *         59655 villeneuve d' ascq cedex               *       c
c       *         france                                       *       c
c       *                                                      *       c
c       *     (3) e.c.m.w.f.                                   *       c
c       *                                                      *       c
c       ********************************************************       c
c                                                                      c
c                                                                      c
c**********************************************************************c
 
c**********************************************************************c
c       ********************************************************       c
c       *                limits of validity                    *       c
c       *                                                      *       c
c       *   geometrical parameters    solar zenith angle and   *       c
c       *                             satellite zenith angle   *       c
c       *                             must be less than 60 and *       c
c       *                             50 degrees respectively. *       c
c       *                                                      *       c
c       *   atmospheric model         no limitations           *       c
c       *                                                      *       c
c       *   aerosol model             the visibility must be   *       c
c       *                             better than 5.0km        *       c
c       *                             for smaller values       *       c
c       *                             calculations might be    *       c
c       *                             no more valid.           *       c
c       *                                                      *       c
c       *   spectral conditions       the gaseous transmittance*       c
c       *                             and the scattering func  *       c
c       *                             tions are valid from 0.25*       c
c       *                             to 4.0 micron. but the   *       c
c       *                             treatment of interaction *       c
c       *                             between absorption and   *       c
c       *                             scattering is correct for*       c
c       *                             not too large absorption *       c
c       *                             if you want to compute   *       c
c       *                             signal within absorption *       c
c       *                             bands,this interaction   *       c
c       *                             ought to be reconsidered *       c
c       *                                                      *       c
c       *   ground reflectance (type) you can consider a patchy*       c
c       *                             structure:that is a circu*       c
c       *                             lar target of radius rad *       c
c       *                             and of reflectance roc,  *       c
c       *                             within an environnement  *       c
c       *                             of reflectance roe.      *       c
c       *                                                      *       c
c       *   ground reflectance (type continued): for uniform   *       c
c       *                             surface conditions only, *       c
c       *                             you may consider directio*       c
c       *                             nal reflectance as bounda*       c
c       *                             ry conditions.           *       c
c       *                             some analytical model are*       c
c       *                             proposed, the user can   *       c
c       *                             specify his own values.  *       c
c       *                             the code assumes that the*       c
c       *                             brdf is spectrally inde- *       c
c       *                             pendent                  *       c
c       *                                                      *       c
c       *   ground reflectance (spectral variation) four typi  *       c
c       *                             cal reflectances are pro *       c
c       *                             posed, defined within    *       c
c       *                             given spectral range.    *       c
c       *                             this range differs accor *       c
c       *                             ding to the selected case*       c
c       *                             the reflectance is set to*       c
c       *                             0 outside this range,due *       c
c       *                             to the deficiency of data*       c
c       *                             user must verify these   *       c
c       *                             limits. that is obviously*       c
c       *                             irrelevant for brdf      *       c
c       *                                                      *       c
c       ********************************************************       c
c**********************************************************************c
 
c****************************************************************************c
c  for considering brdf< we have to compute the downward radiance in the     c
c  whole hemisphere. to perform such computions, we selected the successive  c
c  orders of scattering method. that method requires numerical integration   c
c  over angles and optical depth. the integration method is the gauss method,c
c  mu is the number of angles nmu+1, nmu is settled to 24. the accuracy of   c
c  the computations is obviously depending on the nmu value. this value      c
c  can be easily changed as a parameter as well as the nt value which        c
c  is the number of layers for performing the vertical integration. the      c
c  downward radiance is computed for nmu values of the zenith angle and np   c
c  values of the azimuth angle. the integration of the product of the        c
c  radiance by the brdf is so performed over the nmu*np values. np is settledc
c  to 13, that value can be also changed. mu2 is equal to 2 times nmu.       c
c  xlmus is the downward radiance, xf the downward irradiance, rm and gb     c
c  the angles and the weights for the gauss integration over the zenith, rp  c
c  and bp respectively for the azimuth integration.                          c
c****************************************************************************c
C
C Added by B.-C. Gao in June 1996
C  Declarations for common variables
      DIMENSION HHH(25),TTT(25),PPP(25),VMRR(25)
      DIMENSION WLTEMP(1050),ROTEMP(1050),DTTEMP(1050),ASTEMP(1050)
      DIMENSION WAVOBS(1024),FWHM(1024)
      DIMENSION ROTOT(1050), TTOT(1050), STOT(1050)

C  Common variables used by 6S
      COMMON /GETINPUT3/ HHH,TTT,PPP,VMRR,NB,NL,MODEL,IAER,V,TAER55,
     & VRTO3,SNO2
      COMMON /GETINPUT4/ WAVOBS,FWHM
      COMMON /GETINPUT8/ IMNN,IDYY,IYRR,IHH,IMM,ISS
      COMMON /GEOMETRY1/ SOLZNI,SOLAZ,OBSZNI,OBSPHI,IDAY
      COMMON /SIXS1/ ROTOT, TTOT, STOT
C
C Parameters for plane observation and for elevated surface scenes,
C  surface and plane altitudes are in units of km 
      REAL XPSS, XPPP
      COMMON /GETINPUT14/ XPSS, XPPP
C
C  Local indecies used in 6S
      INTEGER JINDEX, NELEM
 
C--      parameter(nt_p=13,mu_p=13,mu2_p=24,np_p=25)
      parameter(nt_p=26,mu_p=25,mu2_p=48,np_p=49)
      dimension anglem(mu2_p),weightm(mu2_p),
     s   rm(-mu_p:mu_p),gb(-mu_p:mu_p),rp(np_p),gp(np_p)
      dimension  xlmus(-mu_p:mu_p,np_p),xlmuv(-mu_p:mu_p,np_p)
      dimension angmu(10),angphi(13),brdfints(-mu_p:mu_p,np_p)
     s    ,brdfdats(10,13),
     s    brdfintv(-mu_p:mu_p,np_p),brdfdatv(10,13),robar(1501),
     s    robarp(1501),robard(1501),xlm1(-mu_p:mu_p,np_p),
     s    xlm2(-mu_p:mu_p,np_p)
 
        real anglem,weightm,rm,gb,accu2,accu3
        real rp,gp,xlmus,xlmuv,angmu,angphi,brdfints,brdfdats
        real brdfintv,brdfdatv,robar,robarp,robard,xlm1,xlm2
        real c,wldisc,ani,anr,aini,ainr,rocl,roel,zpl,ppl,tpl,whpl
        real wopl,xacc,phasel,pdgs,cgaus,pha,betal,s,wlinf,wlsup,delta
        real sigma,z,p,t,wh,wo,ext,ome,gasym,phase,roatm,dtdir
        real dtdif,utdir,utdif,sphal,wldis,trayl,traypl,pi,pi2,step
        real asol,phi0,avis,phiv,tu,xlon,xlat,xlonan,hna,dsol,campm
        real phi,phirad,xmus,xmuv,xmup,xmud,adif,uw,uo3,taer55
        real taer,v,xps,uwus,uo3us,xpp,taer55p,puw,puo3,puwus
        real puo3us,wl,wlmoy,tamoy,tamoyp,pizmoy,pizmoyp,trmoy
        real trmoyp,fr,rad,spalt
        real albbrdf,par1,par2,par3,par4,robar1,xnorm1,rob,xnor,rodir
        real rdown,rdir,robar2,xnorm2,ro,roc,roe,rapp,rocave,roeave
        real seb,sbor,swl,sb,refet,refet1,refet2,refet3,alumet
        real tgasm,rog,dgasm,ugasm,sdwava,sdozon,sddica,sdoxyg
        real sdniox,sdmoca,sdmeth,suwava,suozon,sudica,suoxyg
        real suniox,sumoca,sumeth,stwava,stozon,stdica,stoxyg,stniox
        real stmoca,stmeth,sodray,sodaer,sodtot,fophsr,fophsa,sroray
        real sroaer,srotot,ssdaer,sdtotr,sdtota,sdtott,sutotr,sutota
        real sutott,sasr,sasa,sast,dtozon,dtdica,dtoxyg
        real dtniox,dtmeth,dtmoca,utozon,utdica,utoxyg,utniox
        real utmeth,utmoca,attwava,ttozon,ttdica,ttoxyg,ttniox
        real ttmeth,ttmoca,dtwava,utwava,ttwava,coef,romix,rorayl
        real roaero,phaa,phar,tsca,tray,trayp,taerp,dtott,utott
        real astot,asray,asaer,utotr,utota,dtotr,dtota,dgtot,tgtot
        real tgp1,tgp2
        real ugtot,edifr,edifa,tdird,tdiru,tdifd,tdifu,fra
        real fae,avr,romeas1,romeas2,romeas3,alumeas,sodrayp
        real ratm1,ratm2,ratm3,rsurf
        real sodaerp,sodtotp,tdir,tdif,etn,esn,es,ea0n,ea0,ee0n
        real ee0,tmdir,tmdif,xla0n,xla0,xltn,xlt,xlen,xle,pizera
        real fophst,pizerr,pizert,xrad,xa,xb,xc
        real sha,sham
        integer nt,mu,mu2,np,k,iwr,mum1,idatmp
        integer j,iread,l,igeom,month,jday,nc,nl,idatm,iaer,iaerp,n
        integer iwave,iinf,isup,ik,i,inhomo,idirec,ibrdf,igroun
        integer igrou1,igrou2,isort
c***********************************************************************
c                             return to 6s
c***********************************************************************
      dimension c(4),wldisc(10),ani(2,3),anr(2,3),aini(2,3),ainr(2,3)
      dimension rocl(1501),roel(1501)
      real rn,ri,x1,x2,x3,cij,rsunph,nrsunph,rmax,rmin
      integer icp,irsunph,i1,i2
      character etiq1(8)*60,nsat(47)*17,atmid(7)*51,reflec(8)*71
      character FILE*80,FILE2*80
 
      logical ier
      common/sixs_ier/iwr,ier
      common /mie_in/ rmax,rmin,icp,rn(10,4),ri(10,4),x1(4),x2(4),
     s x3(4),cij(4),irsunph,rsunph(50),nrsunph(50)
c***********************************************************************
c     for considering pixel and sensor  altitude
c***********************************************************************
      real pps,palt,ftray
      common /sixs_planesim/zpl(34),ppl(34),tpl(34),whpl(34),wopl(34)
      common /sixs_test/xacc
c***********************************************************************
c     for considering brdf
c***********************************************************************
      common /sixs_sos/phasel(10,83),cgaus(83),pdgs(83)
      common /sixs_trunc/pha(83),betal(0:80)
      real optics(3),struct(4)
      integer options(5)
      integer pild,pihs
      real pxLt,pc,pRl,pTl,pRs
      real pws,phi_wind,xsal,pcl,paw
      real uli,eei,thmi,sli,cabi,cwi,vaii,rnci,rsl1i
c***********************************************************************
c                             return to 6s
c***********************************************************************
      common /sixs_ffu/s(1501),wlinf,wlsup
      common /sixs_del/ delta,sigma
      common /sixs_atm/z(34),p(34),t(34),wh(34),wo(34)
      common /sixs_aer/ext(10),ome(10),gasym(10),phase(10)
      common /sixs_disc/ roatm(3,10),dtdir(3,10),dtdif(3,10),
     s utdir(3,10),utdif(3,10),sphal(3,10),wldis(10),trayl(10),
     s traypl(10)
 
 
c****************************************************************************c
c   angmu and angphi are the angles were the brdf is measured. these values  c
c   can be changed as soon as they are well distributed over the whole space c
c   before the gauss integration, these values are interpolated to the gauss c
c   angles                                                                   c
c****************************************************************************c
      data angmu /85.0,80.0,70.0,60.0,50.0,40.0,30.0,20.0,10.0,0.00/
      data angphi/0.00,30.0,60.0,90.0,120.0,150.0,180.0,
     s          210.0,240.0,270.0,300.0,330.0,360.0/
 
c***********************************************************************
c                             return to 6s
c***********************************************************************
      data wldisc /0.400,0.488,0.515,0.550,0.633,
     s             0.694,0.860,1.536,2.250,3.750/
 
      data etiq1/
     s '(1h*,22x,34h user defined conditions          ,t79,1h*)',
     s '(1h*,22x,24h meteosat observation   ,t79,1h*)          ',
     s '(1h*,22x,25h goes east observation   ,t79,1h*)         ',
     s '(1h*,22x,25h goes west observation   ,t79,1h*)         ',
     s '(1h*,22x,30h avhrr (AM noaa) observation  ,t79,1h*)    ',
     s '(1h*,22x,30h avhrr (PM noaa) observation  ,t79,1h*)    ',
     s '(1h*,22x,24h h.r.v.   observation   ,t79,1h*)          ',
     s '(1h*,22x,24h t.m.     observation   ,t79,1h*)          '/
 
       data nsat/
     s ' constant        ',' user s          ',
     s ' meteosat        ',' goes east       ',' goes west       ',
     s ' avhrr 1 (noaa6) ',' avhrr 2 (noaa6) ',
     s ' avhrr 1 (noaa7) ',' avhrr 2 (noaa7) ',
     s ' avhrr 1 (noaa8) ',' avhrr 2 (noaa8) ',
     s ' avhrr 1 (noaa9) ',' avhrr 2 (noaa9) ',
     s ' avhrr 1 (noaa10)',' avhrr 2 (noaa10)',
     s ' avhrr 1 (noaa11)',' avhrr 2 (noaa11)',
     s ' hrv1 1          ',' hrv1 2          ',' hrv1 3          ',
     s ' hrv1 pan        ',
     s ' hrv2 1          ',' hrv2 2          ',' hrv2 3          ',
     s ' hrv2 pan        ',
     s '  tm  1          ','  tm  2          ','  tm  3          ',
     s '  tm  4          ','  tm  5          ','  tm  7          ',
     s '  mss 4          ','  mss 5          ',
     s '  mss 6          ','  mss 7          ',
     s '  mas 1          ','  mas 2          ','  mas 3          ',
     s '  mas 4          ','  mas 5          ','  mas 6          ',
     s '  mas 7          ','  modis 3        ','  modis 5        ',
     s '  modis 6        ',
     s ' avhrr 1 (noaa14)',' avhrr 2 (noaa14)'/
 
      data atmid /
     s 'no absorption computed                             ',
     s 'tropical            (uh2o=4.12g/cm2,uo3=.247cm-atm)',
     s 'midlatitude summer  (uh2o=2.93g/cm2,uo3=.319cm-atm)',
     s 'midlatitude winter  (uh2o=.853g/cm2,uo3=.395cm-atm)',
     s 'subarctic  summer   (uh2o=2.10g/cm2,uo3=.480cm-atm)',
     s 'subarctic  winter   (uh2o=.419g/cm2,uo3=.480cm-atm)',
     s 'us  standard 1962   (uh2o=1.42g/cm2,uo3=.344cm-atm)'/
 
      data  reflec /
     & '(1h*,12x,39h user defined spectral reflectance     ,f6.3,t79
     & ,1h*) ',
     & '(1h*,12x,27h monochromatic reflectance ,f6.3,t79,1h*)',
     & '(1h*,12x,39h constant reflectance over the spectra ,f6.3,t79
     & ,1h*) ',
     & '(1h*,12x,39h spectral vegetation ground reflectance,f6.3,t79
     & ,1h*) ',
     & '(1h*,12x,39h spectral clear water reflectance      ,f6.3,t79
     & ,1h*) ',
     & '(1h*,12x,39h spectral dry sand ground reflectance  ,f6.3,t79
     & ,1h*) ',
     & '(1h*,12x,39h spectral lake water reflectance       ,f6.3,t79
     & ,1h*) ',
     & '(1h*,12x,39h spectral volcanic debris reflectance  ,f6.3,t79
     & ,1h*) '/

      FILE='  '
      FILE2='  '

c***********************************************************************
c   Parameters  initialization
c***********************************************************************
      nt=nt_p
      mu=mu_p
      mu2=mu2_p
      np=np_p
      iwr=6
      ier=.FALSE.
      iinf=1
      isup=1501
c***********************************************************************
c  preliminary computations for gauss integration
c***********************************************************************
      pi=acos(-1.)
      pi2=2*pi
      accu2=1.E-03
      accu3=1.E-07
      do k=1,13
       angphi(k)=angphi(k)*pi/180.
      enddo
      do k=1,10
       angmu(k)=cos(angmu(k)*pi/180.)
      enddo
      call gauss(-1.,1.,anglem,weightm,mu2)
      call gauss(0.,pi2,rp,gp,np)
      mum1=mu-1
      do 581 j=-mum1,-1
       k=mu+j
       rm(-j-mu)=anglem(k)
       gb(-j-mu)=weightm(k)
  581 continue
      do 582 j=1,mum1
       k=mum1+j
       rm(mu-j)=anglem(k)
       gb(mu-j)=weightm(k)
  582 continue
      gb(-mu)=0.
      gb(0)=0.
      gb(mu)=0.
 
c***********************************************************************
c                             return to 6s
c***********************************************************************
c constantes values
      sigma=0.056032
      delta=0.0279
      xacc=1.e-06
      iread=5
      step=0.0025
      do 1111 l=1,10
       wldis(l)=wldisc(l)
 1111 continue
 
c**********************************************************************c
c       igeom               geometrical conditions                     c
c               --------------------------------------                 c
c                                                                      c
c                                                                      c
c   you choose your own conditions; igeom=0                            c
c         0     enter solar zenith angle   (in degrees )               c
c                     solar azimuth angle        "                     c
c                     satellite zenith angle     "                     c
c                     satellite azimuth angle    "                     c
c                     month                                            c
c                     day of the month                                 c
c**********************************************************************c
C---      read(iread,*) igeom
      igeom = 0

c   igeom=0.....
C---      read(iread,*) asol,phi0,avis,phiv,month,jday
      ASOL = SOLZNI * 57.2958
      PHI0 = SOLAZ  * 57.2958
      AVIS = OBSZNI * 57.2958
      PHIV = OBSPHI * 57.2958
C---
C---      print*,'ASOL,PHI0,AVIS,PHIV=',ASOL,PHI0,AVIS,PHIV
C--- month and jday are used in 6S, but the output calculated
C     using these quantities are not used in ATREM. 
      month = IMNN  
      jday =  IDYY 

      if(ier) stop
      dsol=1.
      call varsol(jday,month,
     1                           dsol)
C
C--      print*,'dsol = ',dsol
 
c**********************************************************************c
c                                                                      c
c                                 / scattered direction                c
c                               /                                      c
c                             /                                        c
c                           / adif                                     c
c    incident   + + + + + + + + + + + + + + +                          c
c    direction                                                         c
c                                                                      c
c**********************************************************************c
      phi=abs(phiv-phi0)
      phirad=(phi0-phiv)*pi/180.
      if (phirad.lt.0.) phirad=phirad+2.*pi
      if (phirad.gt.(2.*pi)) phirad=phirad-2.*pi
      xmus=cos(asol*pi/180.)
      xmuv=cos(avis*pi/180.)
      xmup=cos(phirad)
      xmud=-xmus*xmuv-sqrt(1.-xmus*xmus)*sqrt(1.-xmuv*xmuv)*xmup
c test vermote bug
      if (xmud.gt.1.) xmud=1.
      if (xmud.lt.-1.) xmud=-1.
      adif=acos(xmud)*180./pi
 
c**********************************************************************c
c       idatm      atmospheric model                                   c
c                 --------------------                                 c
c                                                                      c
c                                                                      c
c  you select one of the following standard atmosphere: idatm=0 to 6   c
c         0    no gaseous absorption                                   c
c         1    tropical                )                               c
c         2    midlatitude summer      )                               c
c         3    midlatitude winter      )                               c
c         4    subarctic summer        )      from lowtran             c
c         5    subarctic winter        )                               c
c         6    us standard 62          )                               c
c                                                                      c
c  or you define your own atmospheric model idatm=7 or 8               c
c         7    user profile  (radiosonde data on 34 levels)            c
c              enter altitude       (  in km )                         c
c                    pressure       (  in mb )                         c
c                    temperature    (  in k  )                         c
c                    h2o density    (in  g/m3)                         c
c                    o3  density    (in  g/m3)                         c
c                                                                      c
c           for example, altitudes are  from  0 to 25km step of 1km    c
c                        from 25 to 50km step of 5km                   c
c                        and two values at 70km and 100km              c
c                        so you have 34*5 values to input.             c
c         8    enter water vapor and ozone contents                    c
c                 uw  (in  g/cm2 )                                     c
c                 uo3 (in  cm-atm)                                     c
c                 profil is taken from us62                            c
c                                                                      c
c**********************************************************************c
C---      uw=0.
C---      uo3=0.
C---      read(iread,*) idatm
C---      if(idatm.eq.0) go to 5
C---      if(idatm.eq.8) read(iread,*) uw,uo3

          IDATM = MODEL

C--uw and uo3 are not used in computing scattering quantities, simply
C--   assign values to let the program run.
          uw    = 3.0 
          uo3   = 0.35

C***Note: in the final modified code. uw = 0.0, uo3 = 0.0, because 
C           we use 6S only for scattering effect calculations, no
C           gas absorption is needed.

C--      if(idatm.ne.7) go to 6
C--      do 7 k=1,34
C--       read(iread,*) z(k),p(k),t(k),wh(k),wo(k)
C--    7 continue
C--      go to 5
      if(idatm.EQ.7) idatm = 6

    6 if(idatm.eq.1)  call tropic
      if(idatm.eq.2)  call midsum
      if(idatm.eq.3)  call midwin
      if(idatm.eq.4)  call subsum
      if(idatm.eq.5)  call subwin
      if(idatm.eq.6)  call us62
c     we have to define an atmosphere to compute rayleigh optical depth
    5 if(idatm.eq.0.or.idatm.eq.8)  call us62
 
 
c**********************************************************************c
c                                                                      c
c       iaer       aerosol model(type)                                 c
c                  --------------                                      c
c                                                                      c
c                                                                      c
c  you select one of the following standard aerosol models:            c
c         0  no aerosols                                               c
c         1  continental model  )                                      c
c         2  maritime model     )  according to sra models             c
c         3  urban model        )                                      c
c         5  shettle model for background desert aerosol               c
c         6  biomass burning                                           c
c         7  stratospheric model                                       c
c                                                                      c
c  or you define your own model using basic components: iaer=4         c
c         4 enter the volumic percentage of each component             c
c                 c(1) = volumic % of dust-like                        c
c                 c(2) = volumic % of water-soluble                    c
c                 c(3) = volumic % of oceanic                          c
c                 c(4) = volumic % of soot                             c
c                   between 0 to 1                                     c
c                                                                      c
c  or you define your own model using size distribution function:      c 
c         8  Multimodal Log Normal distribution (up to 4 modes)        c
c         9  Modified gamma  distribution                              c
c        10  Junge Power-Law distribution                              c
c                                                                      c
c  or you define a model using sun-photometer measurements:            c
c        11  Sun Photometer  distribution (50 values max)              c
c             you have to enter:  r and d V / d (logr)                 c
c                  where r is the radius (in micron) and V the volume  c
c                  and d V / d (logr) in (cm3/cm2/micron)              c
c             and then you have to enter: nr and ni for each wavelengthc
c                  where nr and ni are respectively the real and       c
c                  imaginary part of the refractive index              c
c                                                                      c
c  or you can use results computed and previously saved                c
c        12  Reading of data previously saved into FILE                c
c             you have to enter the identification name FILE in the    c
c             next line of inputs.                                     c
c                                                                      c
c                                                                      c
c  iaerp and FILE  aerosol model(type)-Printing of results             c
c                  ---------------------------------------             c
c                                                                      c
c For iaer=8,9,10,and 11:                                              c
c    results from the MIE subroutine may be saved into the file        c
c    FILE.mie (Extinction and scattering coefficients, single          c
c    scattering albedo, Asymmetry parameter, phase function at         c
c    predefined wavelengths) and then can be re-used with the          c 
c    option iaer=12 where FILE is an identification name you           c
c    have to enter.                                                    c
c                                                                      c
c    So, if you select iaer=8,9,10,or 11, next line following the      c
c    requested inputs by the options 8,9,10, or 11 you have to enter   c
c    iaerp                                                             c
c                                                                      c
c        iaerp=0    results will not be saved                          c
c        iaerp=1    results will be saved into the file FILE.mie       c
c                    next line enter FILE                              c
c                                                                      c
c                                                                      c
c   example for iaer and iaerp                                         c
c 8                      Multimodal Log-Normale distribution selected  c
c 0.0001 100.0 3         Rmin, Rmax, 3 components                      c
c 0.5000 2.99 1.66E-7    Rmean, Sigma, percentage density-1st componentc
c 1.53 1.53 1.53 1.53 1.53 1.53 1.52 1.40 1.22 1.27  nr-10 wavelengths c 
c .008 .008 .008 .008 .008 .008 .008 .008 .009 .011  ni-10 wavelengths c
c 0.0050 2.99 0.5945     Rmean, Sigma, percentage density-2nd componentc
c 1.53 1.53 1.53 1.53 1.53 1.53 1.52 1.51 1.42 1.452 nr-10 wavelengths c
c .005 .005 .005 .005 .006 .007 .012 .023 .010 .004  ni-10 wavelengths c
c 0.0118 2.00 0.4055     Rmean, Sigma, percentage density-3rd componentc
c 1.75 1.75 1.75 1.75 1.75 1.75 1.75 1.77 1.81 1.90  nr-10 wavelengths c
c .46  .45  .45  .44  .43  .43  .43  .46  .50  .57   ni-10 wavelengths c
c 1                      Results will be saved into FILE.mie           c
c URBAN-WCP112           Identification of the output file called FILE c
c                    -> results will be saved into URBAN-WCP112.mie    c
c                                                                      c
c**********************************************************************c
      rmin=0.
      rmax=0.
      icp=1
      do i=1,4
       x1(i)=0.0
       x2(i)=0.0
       x3(i)=0.0
       do l=1,10
        rn(l,i)=0.0
        ri(l,i)=0.0
       enddo
      enddo
      do i=1,50
       rsunph(i)=0.
       nrsunph(i)=0.
      enddo
      cij(1)=1.00
 
      if(iaer.eq.4) read(iread,*) (c(n),n=1,4)
      goto(49,40,41,42,49,49,49,49,43,44,45,46,47),iaer+1
 
   40 c(1)=0.70
      c(2)=0.29
      c(3)=0.00
      c(4)=0.01
      go to 49
   41 c(1)=0.00
      c(2)=0.05
      c(3)=0.95
      c(4)=0.00
      go to 49
   42 c(1)=0.17
      c(2)=0.61
      c(3)=0.00
      c(4)=0.22
      go to 49
   43 read(iread,*) rmin,rmax,icp
      do i=1,icp
       read(5,*)x1(i),x2(i),cij(i)
       read(5,*)(rn(l,i),l=1,10)
       read(5,*)(ri(l,i),l=1,10)
      enddo
      go to 49
   44 read(iread,*) rmin,rmax
      read(iread,*) x1(1),x2(1),x3(1)
      read(5,*)(rn(l,1),l=1,10)
      read(5,*)(ri(l,1),l=1,10)
      go to 49
   45 read(iread,*) rmin,rmax
      read(iread,*) x1(1)
      read(5,*)(rn(l,1),l=1,10)
      read(5,*)(ri(l,1),l=1,10)
      go to 49
   46 read(5,*)irsunph
      do i=1,irsunph
       read(5,*)rsunph(i),nrsunph(i)
       nrsunph(i)=nrsunph(i)/(rsunph(i)**4.)/alog(10.0)
      enddo
      rmin=rsunph(1)
      rmax=rsunph(irsunph)+1e-07
      read(5,*)(rn(l,1),l=1,10)
      read(5,*)(ri(l,1),l=1,10)
      go to 49
   47 read(5,'(A80)')FILE2
      i2=index(FILE2,' ')-1
      go to 49

   49 continue
      if (iaer.ge.8.and.iaer.le.11)then
       read(5,*)iaerp
       if (iaerp.eq.1)read(5,'(A80)')FILE
       i1=index(FILE,' ')-1
       FILE2=FILE(1:I1)//'.mie'
       i2=index(FILE2,' ')-1
      endif

      call aeroso(iaer,c,xmud,wldis,FILE2)
 
c**********************************************************************c
c              aerosol model (concentration)                           c
c              ----------------------------                            c
c                                                                      c
c                                                                      c
c   you have an estimate of the meteorological parameter: the visibi   c
c   lity v,                                                            c
c            enter directly the value of v in km(the aerosol optical   c
c            depth will be computed from a standard aerosol profile)   c
c                                                                      c
c   or you have an estimate of aerosol optical depth ,enter v=0 for    c
c            the visibility and enter the aerosol optical depth at 550 c
c                                                                      c
c   warning:  if iaer=0, enter v=-1                                    c
c**********************************************************************c
      taer55=0.
      taer=0.
C----Temp note:   Here another index should be used in the final code.
C      the current logic is:
C        READ(iread,*), ID_Taer, Taer55
C        IF(ID_Taer) 71, 10, 11
C in the users manual we should state:
C       IF(ID_Taer.LT.0), Taer55 must = 0.0
C       IF(ID_Taer.EQ.0), then v = exp(-log(taer55/2.7628)/0.79902)
C       IF(ID_Taer.GT.0), then assign v = taer55 and call oda550(iaer,v,taer55)
C---End of Temp Note****
C--      read(iread,*) v
C--       v = 23.
      if(v) 71,10,11
   10 read(iread,*) taer55
      v=exp(-log(taer55/2.7628)/0.79902)
      goto 71
   11 call oda550(iaer,v,taer55)
   71 continue
 
c**********************************************************************c
c xps is the parameter to express the  altitude of target              c
c                                                                      c
c                                                                      c
c                  xps >0. means you know the altitude of the target   c
c                        expressed in km and you put that value as xps c
c                                                                      c
c                                                                      c
c**********************************************************************c
            xps = XPSS
C--         print*,'XPSS =', XPSS, xps
C---       read(iread,*) xps
C---       if (xps.ge.0.) then
       if (xps.le.0.) then
        xps=0.
        uwus=1.424
        uo3us=0.344
       else
        if (idatm.ne.8) then
         call pressure(uw,uo3,xps)
        else
         call pressure(uwus,uo3us,xps)
        endif
       endif
 
c**********************************************************************c
c                                                                      c
c  xpp is the parameter to express the sensor altitude                 c
c                                                                      c
c                                                                      c
c           xpp= 100  means that the sensor is a board a satellite     c
c           xpp=    0 means that the sensor is at the ground level     c
c                                                                      c
c                                                                      c
c     for aircraft simulations                                         c
c     0< xpp <100  means you know the altitude of the sensor expressed c
c                  in kilometers units      			       c
C     this altitude is the ABSOLUTE target altitude (relative to the   c
C                  SEA level.                                          c
c                                                                      c
c     for aircraft simulations only, you have to give                  c
c	puw,po3   (water vapor content,ozone content between the       c
c                  aircraft and the surface)                           c
c	taerp     (the aerosol optical thickness at 550nm between the  c
c                  aircraft and the surface)                           c
c    if these data are not available, enter negative values for all    c
c    of them, puw,po3 will then be interpolated from the us62 standard c
C    profile according to the values at ground level. Taerp will be    c
c    computed according to a 2km exponential profile for aerosol.      c
c**********************************************************************c
c
            xpp = XPPP
C--         print*, 'XPPP = ', XPPP, xpp
C---        read(iread,*) xpp
C---        xpp=-xpp
C---Note because of changing xpp to absolute plane altitude (relative
C        to the sea level) and in order to conform internal 6S 
C        computation scheme, we have to set xpp = xpp - xps
            xpp = xpp - xps
C---Now xpp is relative to the "ELEVATED" ground target level --------
        if (xpp.lt.0.0) then
c          ground measurement option        
           palt=0.
           pps=p(1)
	   idatmp=0
	   taer55p=0.
	   puw=0.
	   puoz=0.
           else
	   if (xpp.gt.100.) then
c	       satellite case of equivalent	   
	      palt=1000.
	      pps=0.
	      taer55p=taer55
	      ftray=1.
	      idatmp=4
	      else
c	      "real" plane case	      
C--Temp Code (can remain here, because output not depend on puw & puo3):
                 puw  = -2.5
                 puo3 = -0.03
C--            read(iread,*) puw,puo3
C--End of Temp Code   
	      if (puw.lt.0.) then
                 call presplane(puw,puo3,xpp,ftray)
	         idatmp=2
	         if (idatm.eq.8) then
	            puwus=puw
	            puo3us=puo3
	            puw=puw*uw/uwus
	            puo3=puo3*uo3/uo3us
	            idatmp=8
	         endif
	         else
	         call presplane(puwus,puo3us,xpp,ftray)
	         idatmp=8
              endif
              if(ier) stop
              palt=zpl(34)
	      pps=ppl(34)
C--Temp Code (can remain here to automatically distribute aerosol amounts
C             below and above the airplane):
              taer55p = -0.2
C---              read(iread,*) taer55p
C--End of Temp Code   
	    if ((taer55p.lt.0.).or.((taer55-taer55p).lt.accu2)) then
c a scale heigh of 2km is assumed in case no value is given for taer55p
               taer55p=taer55*(1.-exp(-palt/2.))
            else
C compute effective scale heigh
               sham=exp(-palt/4.)
               sha=1.-(taer55p/taer55)
               if (sha.ge.sham) then
                  taer55p=taer55*(1.-exp(-palt/4.))
               else
                  sha=-palt/log(sha)
                  taer55p=taer55*(1.-exp(-palt/sha))
               endif
            endif
         endif
      endif
 
c**********************************************************************c
c      iwave input of the spectral conditions                          c
c            --------------------------------                          c
c                                                                      c
c  you choose to define your own spectral conditions: iwave=-1,0 or 1  c
c                   (three user s conditions )                         c
c        -2  enter wlinf, wlsup, the filter function will be equal to 1c
c            over the whole band (as iwave=0) but step by step output  c
c            will be printed                                           c
c        -1  enter wl (monochr. cond,  gaseous absorption is included) c
c                                                                      c
c         0  enter wlinf, wlsup. the filter function will be equal to 1c
c            over the whole band.                                      c
c                                                                      c
c         1  enter wlinf, wlsup and user's filter function s(lambda)   c
c                          ( by step of 0.0025 micrometer).            c
c                                                                      c
c  note: wl has to be in micrometer                                    c
c**********************************************************************c
      do 38 l=iinf,isup
       s(l)=1.
   38 continue
c
C--Temp code:
      iwave = -2
      wlinf = 0.30
      wlsup = 2.90
C--End of Temp Code.
      iinf=(wlinf-.25)/0.0025+1.5
      isup=(wlsup-.25)/0.0025+1.5
 
c**********************************************************************c
c here, we first compute an equivalent wavelenght which is the input   c
c value for monochromatic conditions or the integrated value for a     c
c filter functionr (vall equivwl) then, the atmospheric properties are c
c computed for that wavelength (call discom then call specinterp)      c
c molecular optical thickness is computed too (call odrayl). lastly    c
c the successive order of scattering code is called three times.       c
c first for a sun at thetas with the scattering properties of aerosols c
c and molecules, second with a pure molecular atmosphere, then with thec
c actual atmosphere for a sun at thetav. the iso code allows us to     c
c compute the scattering transmissions and the spherical albedo. all   c
c these computations are performed for checking the accuracy of the    c
c analytical expressions and in addition for computing the averaged    c
c directional reflectances                                             c
c**********************************************************************c
      if(iwave.ne.-1) then
        call equivwl(iinf,isup,step,
     s               wlmoy)
      else
        wlmoy=wl
      endif
c     write(6,*) "wlmoy: ",wlmoy
      call discom (idatmp,iaer,xmus,xmuv,phi
     a            ,taer55,taer55p,palt,
     a           phirad,nt,mu,np,rm,gb,rp
     s         ,ftray,xlm1,xlm2)
      if(iaer.ne.0) then
        call specinterp(wlmoy,taer55,taer55p,
     s     tamoy,tamoyp,pizmoy,pizmoyp)
      endif
      call odrayl(wlmoy,
     s                   trmoy)
      trmoyp=trmoy*ftray
      if (idatmp.eq.4) then
          trmoyp=trmoy
          tamoyp=tamoy
      endif
      if (idatmp.eq.0) then
         trmoyp=0.
         tamoyp=0.
      endif
 
c*********************************************************************c
c     inhomo        ground reflectance (type)                          c
c                   ------------------                                 c
c                                                                      c
c  you consider an homogeneous surface:                                c
c     enter - inhomo=0                                                 c
c                you may consider directional surface  effects         c
c                  idirec=0 (no directional effect)                    c
c                          you have to specify the surface reflectance:c
c                          igroun  (see note1) which is uniform and    c
c                          lambertian                                  c
c                  idirec=1 ( directional effect)                      c
c                          you have to specify the brdf of the surface c
c                           for the actual solar illumination you  are c
c                           considering as well as the brdf for a sun  c
c                           which would be at an angle thetav, in      c
c                           addition you have to give the surface      c
c                           albedo (spherical albedo). you can also    c
c                           select one of the selected model from the  c
c                           ibrdf value (see note2). 3 reflectances    c
c                           are computed, robar,robarp and robard      c
c                                                                      c
c                            ****tree****                              c
c                                                                      c
c                               inhomo                                 c
c                             /          \                             c
c                            /            \                            c
c                           /              \                           c
c                          /                \                          c
c                 ------- 0 -------       -----1 -----                 c
c                        /               /   \       \                 c
c                    idirec             /     \       \                c
c                    /  \              /       \       \               c
c                   /    \            /         \       \              c
c                  /      \       igrou1       igrou2    rad           c
c                 0        1        roc          roe     f(r)          c
c                /          \                                          c
c               /            \                                         c
c           igroun          ibrdf                                      c
c        (roc = roe)        (roc)                                      c
c                           (robar)                                    c
c                           (robarp)                                   c
c                           (robard)                                   c
c                                                                      c
c                   ground reflectance (spectral variation)            c
c                   ---------------------------------------            c
c note1: values of the reflectance selected by igroun,igrou1 or igrou2 c
c        may correspond to the following cases,                        c
c         0  constant value of ro (or roc,or roe) whatever the wavelen c
c            gth. you enter this constant value of ro (or roc or roe). c
c        -1  you have to enter the value of ro (or roc,or roe) by step c
c            of 0.0025 micron from wlinf to wlsup (if you have used thec
c            satellite bands,see implicit values for these limits).    c
c         1  mean spectral value of green vegetation                   c
c         2  mean spectral value of clear water                        c
c         3  mean spectral value of sand                               c
c         4  mean spectral value of lake water                         c
c**********************************************************************c
									
      fr=0.
      rad=0.
      do 1116 ik=iinf,isup
        rocl(ik)=0.
        roel(ik)=0.
 1116 continue
 
c**********************************************************************c
c     uniform or non-uniform surface conditions                        c
c**********************************************************************c
C---      read(iread,*) inhomo
      inhomo = 0
C---                 if(inhomo) 30,30,31
C---             30  read(iread,*) idirec
      idirec = 0
 
c**********************************************************************c
c     uniform surface with lambertian conditions                       c
c**********************************************************************c
C---  21  read(iread,*) igroun
      igroun = 0
      do 35 l=iinf,isup
        rocl(l) = 0.
   35 continue
c
      do 39 l=iinf,isup
        roel(l)=rocl(l)
   39 continue
 
c**********************************************************************c
 
 
c**********************************************************************c
c                     print of initial conditions                      c
c                                                                      c
c**********************************************************************c
 
c ---- geometrical conditions ----
      write(iwr, 98)
      write(iwr, etiq1(igeom+1))
      if(igeom.eq.0) then
	 write(iwr, 1401)
	 write(iwr, 103)month,jday
      endif
      if(igeom.ne.0) write(iwr, 101)month,jday,tu,xlat,xlon
      write(iwr, 102)asol,phi0
      write(iwr, 1110)avis,phiv,adif,phi
 
c --- atmospheric model ----
      write(iwr, 1119)
      if(idatm-7)226,227,228
  228 write(iwr, 1281)uw,uo3
      goto 219
  227 write(iwr, 1272)
      do 229 i=1,34
        write(iwr, 1271)z(i),p(i),t(i),wh(i),wo(i)
  229 continue
      goto 219
  226 write(iwr, 1261)atmid(idatm+1)
 
c --- aerosols model (type) ----
 219  if (iaer.lt.4) then
        goto(230,231,232,233),iaer+1
      else
        if (iaer.ge.5.and.iaer.le.7) goto(234,235,236),iaer-4
        if (iaer.eq.4)write(iwr,133)(c(i),i=1,4)
        if (iaer.eq.8)then
          write(iwr,134)icp
          do i=1,icp
            write(iwr,135)x1(i),x2(i),cij(i)
          enddo
        endif
        if (iaer.eq.9)write(iwr,136)x1(1),x2(1),x3(1)
        if (iaer.eq.10)write(iwr,137)x1(1)
        if (iaer.eq.11)write(iwr, 131)' Sun Photometer'
        if (iaer.eq.12)write(iwr,138)FILE2(1:i2)
	if (iaerp.eq.1)write(iwr,139)FILE2(1:i2)
        goto 249
      endif
  234 write(iwr, 131)'       Desertic'
      goto 249
  235 write(iwr, 131)'          Smoke'
      goto 249
  236 write(iwr, 131)'  Stratospheric'
      goto 249
  233 write(iwr, 131)'          Urban'
      go to 249
  232 write(iwr, 131)'       Maritime'
      goto 249
  231 write(iwr, 131)'    Continental'
      goto 249
  230 write(iwr, 1301)
  249 continue
 
c --- aerosol model (concentration) ----
      if(iaer.eq.0) write(iwr, 1401)
      if(iaer.eq.0) goto 1112
      if(abs(v).le.xacc) write(iwr, 140)taer55
      if(abs(v).gt.xacc) write(iwr, 141)v,taer55
 
c --- spectral condition ----
 1112 CONTINUE
c--- 1112 write(iwr, 148)
c---      if(iwave.eq.-2) write(iwr, 1510) nsat(1),wlinf,wlsup
c---      if(iwave.eq.-1) write(iwr, 149) wl
c---      if(iwave.ge.0) write(iwr, 1510) nsat(iwave+1), wlinf,wlsup
 
c --- ground reflectance (type and spectral variation) ----
      if(idirec.eq.0) then
        rocave=0.
        roeave=0.
        seb=0.
 
        do 264 i=iinf,isup
          sbor=s(i)
          if(i.eq.iinf.or.i.eq.isup) sbor=sbor*0.5
          wl=.25+(i-1)*step
C---
C---          call solirr(wl,
C---         1            swl)
          swl = 1.0
C---
          swl=swl*dsol
          rocave=rocave+rocl(i)*sbor*swl*step
          roeave=roeave+roel(i)*sbor*swl*step
          seb=seb+sbor*swl*step
  264   continue
        rocave=rocave/seb
        roeave=roeave/seb
        isort=0
        ro=rocave
 
        if(inhomo.eq.0) goto 260
        write(iwr, 169)rad
        igroun=igrou1
        ro=rocave
        write(iwr, 170)
        goto 261
 
  262   igroun=igrou2
        ro=roeave
        write(iwr, 171)
        goto 261
 
  260   CONTINUE
C---  260   write(iwr, 168)
  261   if (igroun.gt.0)write(iwr, reflec(igroun+3))ro
        if (igroun.gt.0)goto 158
        if(igroun.eq.-1) write(iwr, reflec(1))ro
        if(igroun.eq.-1) goto 158
        if(iwave.eq.-1)  write(iwr, reflec(2))ro
C--        if(iwave.ne.-1)  write(iwr, reflec(3))ro
 158    isort=isort+1
        if(inhomo.eq.0) goto 999
        if(isort.eq.2) goto 999
        goto 262
      else
        write(iwr, 168)
      endif

c --- pressure at ground level (174) and altitude (175) ----
  999 write(iwr, 173)
      write(iwr, 174)p(1)
      write(iwr, 175)xps
C---      if (xps.gt.0.) write(iwr, 176)uw,uo3
 
c --- plane simulation output if selected ----
      if (palt.lt.1000.) then
C---       write(iwr, 178)
       write(iwr, 179)pps
       write(iwr, 180)palt
C---       write(iwr, 181)
C---       write(iwr, 182)puo3
C---       write(iwr, 183)puw
C---       write(iwr, 184)taer55p
      endif
 
c**********************************************************************c
c                                                                      c
c                                                                      c
c                     start of computations                            c
c                                                                      c
c                                                                      c
c                                                                      c
c**********************************************************************c

c ---- spectral loop ----
C---      if (iwave.eq.-2) write(iwr,1500)
      do 51 l=iinf,isup
        sbor=s(l)
        if(l.eq.iinf.or.l.eq.isup) sbor=sbor*0.5
        if(iwave.eq.-1) sbor=1.0/step
        roc=rocl(l)
        roe=roel(l)
        wl=.25+(l-1)*step
c---
C---        call solirr(wl,
C---     s            swl)
C---
        call interp (iaer,idatmp,wl,taer55,taer55p,xmud,
     s             romix,rorayl,roaero,phaa,phar,tsca,
     s             tray,trayp,taer,taerp,dtott,utott,
     s             astot,asray,asaer,
     s             utotr,utota,dtotr,dtota)
c
        if (iwave.eq.-2) then
C--          write(iwr,1501) wl,tgtot,dtott,utott,astot,ratm2,swl,step,
C--     s            sbor,dsol,romeas2
C---        write(iwr,1501) wl,tray,dtott,utott,astot,trayp,romix,rorayl,
C---     s            roaero,taer,taerp
C
C--Modified by B.-C. Gao
      JINDEX=L-IINF+1
      WLTEMP(JINDEX) = WL
      ROTEMP(JINDEX) = romix
      DTTEMP(JINDEX) = DTOTT*UTOTT
      ASTEMP(JINDEX) = ASTOT
C---      WRITE(6,*) WLTEMP(JINDEX),ROTEMP(JINDEX),DTTEMP(JINDEX),
C---     &           ASTEMP(JINDEX),JINDEX
C---End of modifications ---
        endif

 51   continue
C
C--Modified by B.-C. Gao ---
C
      NELEM = ISUP-IINF+1
C
C--      DO L = 1, NELEM
C--         WRITE(6,*) WLTEMP(L),ROTEMP(L),DTTEMP(L),
C--     &              ASTEMP(L),L
C--      END DO
        
      CALL CUBSPLN(NELEM,WLTEMP,ROTEMP,WAVOBS,ROTOT)
      CALL CUBSPLN(NELEM,WLTEMP,DTTEMP,WAVOBS,TTOT)
      CALL CUBSPLN(NELEM,WLTEMP,ASTEMP,WAVOBS,STOT)
      WRITE(*,*) 'I,WAVOBS(I), ROTOT(I), TTOT(I), STOT(I)'
      DO 122 I=1,128
 122  WRITE(*,*) I,WAVOBS(I), ROTOT(I), TTOT(I), STOT(I)
C--      DO 123 I=1,ISUP,10
C--  123    WRITE(*,*)'I=',I,' TTOT(I)=',TTOT(I)
C--      DO 124 I=1,ISUP,10
C--  124    WRITE(*,*)'I=',I,' STOT(I)=',STOT(I)
******End of addition

      WRITE(*,9257)
 9257 format(79(1h*),/)

      return
 
c**********************************************************************c
c                                                                      c
c                   output editing formats                             c
c                                                                      c
c                                                                      c
c**********************************************************************c
   98 format(///,1h*,30(1h*),16h 6s version 4.1 ,30(1h*),t79
     s       ,1h*,/,1h*,t79,1h*,/,
     s       1h*,22x,34h geometrical conditions identity  ,t79,1h*,/,
     s       1h*,22x,34h -------------------------------  ,t79,1h*)
  101 format(1h*,15x,7h month:,i3,7h day : ,i3,
     s                 16h universal time:,f6.2,
     s                 10h (hh.dd)  ,t79,1h*,/,
     s   1h*, 15x,10hlatitude: ,f7.2,5h deg ,6x,
     s                 12h longitude: ,f7.2,5h deg ,t79,1h*)
  102 format(1h*,2x,22h solar zenith angle:  ,f6.2,5h deg ,
     s     29h solar azimuthal angle:      ,f6.2,5h deg ,t79,1h*)
  103 format(1h*,2x,7h month:,i3,7h day : ,i3,t79,1h*)
 1110 format(1h*,2x,22h view zenith angle:   ,f6.2,5h deg ,
     s       29h view azimuthal angle:       ,f6.2,5h deg ,
     s      t79,1h*,/,
     s       1h*,2x,22h scattering angle:    ,f6.2,5h deg ,
     s           29h azimuthal angle difference: ,f6.2,5h deg ,
     s      t79,1h*)
 1119 format(1h*,t79,1h*,/,
     s       1h*,22x,31h atmospheric model description ,t79,1h*,/,
     s       1h*,22x,31h ----------------------------- ,t79,1h*)
 1261 format(1h*,10x,30h atmospheric model identity : ,t79,1h*,/,
     s       1h*,15x,a51,t79,1h*)
 1272 format(1h*,30h atmospheric model identity : ,t79,1h*,/,
     s       1h*,12x,33h user defined atmospheric model  ,t79,1h*,/,
     s       1h*,12x,11h*altitude  ,11h*pressure  ,
     s           11h*temp.     ,11h*h2o dens. ,11h*o3 dens.  ,t79,1h*)
 1271 format(1h*,12x,5e11.4,t79,1h*)
 1281 format(1h*,10x,31h atmospheric model identity :  ,t79,1h*,
     s     /,1h*,12x,35h user defined water content : uh2o=,f6.3,
     s                  7h g/cm2 ,t79,1h*,
     s     /,1h*,12x,35h user defined ozone content : uo3 =,f6.3,
     s                  7h cm-atm,t79,1h*)
 1301 format(1h*,10x,25h aerosols type identity :,t79,1h*,/,
     s     1h*,15x,24h no aerosols computed   ,t79,1h*)
  131 format(1h*,10x,25h aerosols type identity :,t79,1h*,/,
     s     1h*,15x,a15,15h aerosols model,t79,1h*)
  133 format(1h*,10x,25h aerosols type identity :,t79,1h*,/,
     s     1h*,15x,29h user defined aerosols model ,t79,1h*,/,
     s  1h*,26x,f6.3,15h % of dust-like,t79,1h*,/,
     s  1h*,26x,f6.3,19h % of water-soluble,t79,1h*,/,
     s  1h*,26x,f6.3,13h % of oceanic,t79,1h*,/,
     s  1h*,26x,f6.3,10h % of soot,t79,1h*)
  134 format(1h*,10x,25h aerosols type identity :,t79,1h*,/,
     s  1h*,15x,29h user defined aerosols model ,t79,1h*,/,
     s  1h*,15x,6husing ,i1,32h Log-normal size-distribution(s),t79,
     s  1h*,/,1h*,15x,42hMean radius  Stand. Dev.  Percent. dencity,
     s  t79,1h*)
  135 format(1h*,T41,f6.4,T55,f5.3,T69,e8.3,T79,1h*)
  136 format(1h*,10x,25h aerosols type identity :,t79,1h*,/,
     s  1h*,15x,29h user defined aerosols model ,t79,1h*,/,
     s  1h*,15x,40husing a Modified Gamma size-distribution,t79,1h*,/,
     s  1h*,19x,33hAlpha         b             Gamma,t79,1h*,/,
     s  1h*,T20,f6.3,T31,f6.3,T47,f6.3,T79,1h*)
  137 format(1h*,10x,25h aerosols type identity :,t79,1h*,/,
     s  1h*,15x,29h user defined aerosols model ,t79,1h*,/,
     s  1h*,15x,47husing a Power law size-distribution with alpha=,
     s  f3.1,T79,1h*)
  138 format(1h*,10x,25h aerosols type identity :,t79,1h*,/,
     s  1h*,15x,29h user defined aerosols model ,t79,1h*,/,
     s  1h*,15x,25husing data from the file:,T79,1h*,/,
     s  1h*,T25,A30,T79,1h*)
  139 format(1h*,15x,29h results saved into the file:,T79,1h*,/,
     s  1h*,T25,A30,T79,1h*)
  140 format(1h*,10x,29h optical condition identity :,t79,1h*,/,
     s       1h*,15x,31h user def. opt. thick. at 550nm :,f7.4,
     s       t79,1h*,/,1h*,t79,1h*)
  141 format(1h*,10x,29h optical condition identity :,t79,1h*,/,
     s       1h*,15x,13h visibility :,f6.2,4h km ,
     s                 20h opt. thick. 550nm :,f7.4,t79,1h*,/,
     s                  1h*,t79,1h*)
  148 format(1h*,22x,21h spectral condition  ,t79,1h*,/,1h*,
     s             22x,21h ------------------  ,t79,1h*)
  149 format(1h*,12x,34h monochromatic calculation at wl :,
     s                              f6.3,8h micron ,t79,1h*)
 1510 format(1h*,10x,a17,t79,1h*,/,
     s 1h*,15x,26hvalue of filter function :,t79,1h*,/,1h*,
     s 15x,8h wl inf=,f6.3,4h mic,2x,8h wl sup=,f6.3,4h mic,t79,1h*)
  168 format(1h*,t79,1h*,/,1h*,22x,14h target type  ,t79,1h*,/,1h*,
     s                         22x,14h -----------  ,t79,1h*,/,1h*,
     s                         10x,20h homogeneous ground ,t79,1h*)
  169 format(1h*,t79,1h*,/,1h*,22x,14h target type  ,t79,1h*,/,1h*,
     s                         22x,14h -----------  ,t79,1h*,/,1h*,
     s    10x,41h inhomogeneous ground , radius of target ,f6.3,
     s         5h km  ,t79,1h*)
  170 format(1h*,15x,22h target reflectance : ,t79,1h*)
  171 format(1h*,15x,29h environmental reflectance : ,t79,1h*)
  172 format(1h*,t79,1h*,/,79(1h*),///)
  173 format(1h*,t79,1h*,/,
     s       1h*,22x,30h target elevation description ,t79,1h*,/,
     s       1h*,22x,30h ---------------------------- ,t79,1h*)
  174 format(1h*,10x,22h ground pressure  [mb]    ,1x,f7.2,1x,t79,1h*)
  175 format(1h*,10x,22h ground altitude  [km]    ,f6.3,1x,t79,1h*)
  176 format(1h*,15x,34h gaseous content at target level: ,t79,1h*,
     s     /,1h*,15x,6h uh2o=,f6.3,7h g/cm2 ,
     s           5x,6h  uo3=,f6.3,7h cm-atm,t79,1h*)

c pressure at ground level (174) and altitude (175)
  178 format(1h*,t79,1h*,/,
     s       1h*,22x,30h plane simulation description ,t79,1h*,/,
     s       1h*,22x,30h ---------------------------- ,t79,1h*)
  179 format(1h*,10x,31h plane  pressure          [mb] ,f7.2,1x,t79,1h*)
  180 format(1h*,10x,31h plane  altitude absolute [km] ,f6.3,1x,t79,1h*)
  181 format(1h*,15x,37h atmosphere under plane description: ,t79,1h*)
  182 format(1h*,15x,26h ozone content            ,f6.3,1x,t79,1h*)
  183 format(1h*,15x,26h h2o   content            ,f6.3,1x,t79,1h*)
  184 format(1h*,15x,26haerosol opt. thick. 550nm ,f6.3,1x,t79,1h*)
 
 1401 format(1h*,t79,1h*)
 1500 format(1h*,1x,42hwave   total  total  total  total  atm.   ,
     s           33hswl    step   sbor   dsol   toar ,t79,1h*,/,
     s  1h*,1x,42h       gas    scat   scat   spheri intr   ,t79,1h*,/,
     s  1h*,1x,42h       trans  down   up     albedo refl   ,t79,1h*)
C--1501 format(1h*,6(F6.4,1X),F6.1,1X,4(F6.4,1X),t79,1h*)
 1501 format(1x,6(F6.4,1X),F6.4,1X,4(F6.4,1X),t79,1x)

      end

      subroutine aeroso (iaer,co,xmud,wldis,FILE)

      double precision cij(4),vi(4),nis,sumni,ni(4)
      real co(4),dd(4,10),ci(4),ex(4,10),sc(4,10),asy(4,10)
      real pha(5,10,83),sca(10),wldis(10)
      real ex2(1,10),sc2(1,10),asy2(1,10)
      real ex3(1,10),sc3(1,10),asy3(1,10)
      real ex4(1,10),sc4(1,10),asy4(1,10)
      real xmud,ext,ome,gasym,phase,ph,phasel,cgaus,pdgs
      real coef,sigm,pi
      integer i,j,k,l,j1,j2,iaer,icp
      character FILE*80
c sra basic components for aerosol model, extinction coefficients are 
c in km-1.
c     dust-like = 1
c     water-soluble = 2
c     oceanique = 3
c     soot = 4
 
      data vi /113.983516,113.983516d-06,5.1444150196,
     a          59.77353425d-06/
      data ni /54.734,1.86855d+06,276.05,1.80582d+06/
 
c     i: 1=dust-like 2=water-soluble 3=oceanic 4=soot
      data ((ex(i,j),sc(i,j),j=1,10),i=1,1) /
     a 0.1796674e-01,0.1126647e-01,0.1815135e-01,0.1168918e-01,
     a 0.1820247e-01,0.1180978e-01,0.1827016e-01,0.1196792e-01,
     a 0.1842182e-01,0.1232056e-01,0.1853081e-01,0.1256952e-01,
     a 0.1881427e-01,0.1319347e-01,0.1974608e-01,0.1520712e-01,
     a 0.1910712e-01,0.1531952e-01,0.1876025e-01,0.1546761e-01/
      data ((ex(i,j),sc(i,j),j=1,10),i=2,2) /
     a 0.7653460e-06,0.7377123e-06,0.6158538e-06,0.5939413e-06,
     a 0.5793444e-06,0.5587120e-06,0.5351736e-06,0.5125148e-06,
     a 0.4480091e-06,0.4289210e-06,0.3971033e-06,0.3772760e-06,
     a 0.2900993e-06,0.2648252e-06,0.1161433e-06,0.9331806e-07,
     a 0.3975192e-07,0.3345499e-07,0.1338443e-07,0.1201109e-07/
      data ((ex(i,j),sc(i,j),j=1,10),i=3,3) /
     a 0.3499458e-02,0.3499455e-02,0.3574996e-02,0.3574993e-02,
     a 0.3596592e-02,0.3596591e-02,0.3622467e-02,0.3622465e-02,
     a 0.3676341e-02,0.3676338e-02,0.3708866e-02,0.3708858e-02,
     a 0.3770822e-02,0.3770696e-02,0.3692255e-02,0.3677038e-02,
     a 0.3267943e-02,0.3233194e-02,0.2801670e-02,0.2728013e-02/
      data ((ex(i,j),sc(i,j),j=1,10),i=4,4) /
     a 0.8609083e-06,0.2299196e-06,0.6590103e-06,0.1519321e-06,
     a 0.6145787e-06,0.1350890e-06,0.5537643e-06,0.1155423e-06,
     a 0.4503008e-06,0.8200095e-07,0.3966041e-06,0.6469735e-07,
     a 0.2965532e-06,0.3610638e-07,0.1493927e-06,0.6227224e-08,
     a 0.1017134e-06,0.1779378e-08,0.6065031e-07,0.3050002e-09/
 
       data ((ex2(i,j),sc2(i,j),j=1,10),i=1,1) /
     A 0.4383631E+02,0.4028625E+02,0.4212415E+02,0.3904473E+02,        
     A 0.4157425E+02,0.3861470E+02,0.4085399E+02,0.3803645E+02,        
     A 0.3914040E+02,0.3661054E+02,0.3789763E+02,0.3554456E+02,        
     A 0.3467506E+02,0.3269951E+02,0.2459000E+02,0.2341019E+02,        
     A 0.1796726E+02,0.1715375E+02,0.1057569E+02,0.1009731E+02/        

       data ((ex3(i,j),sc3(i,j),j=1,10),i=1,1) /
     A 0.9539786E+05,0.9297790E+05,0.7530360E+05,0.7339717E+05,
     A 0.7021064E+05,0.6842549E+05,0.6421828E+05,0.6257180E+05,
     A 0.5243056E+05,0.5104987E+05,0.4557768E+05,0.4434877E+05,
     A 0.3193777E+05,0.3100621E+05,0.9637680E+04,0.9202678E+04,
     A 0.3610691E+04,0.3344476E+04,0.8105614E+03,0.6641915E+03/
 
       data ((ex4(i,j),sc4(i,j),j=1,10),i=1,1) /
     A .5427304E+08, .5427304E+08, .6198144E+08, .6198144E+08,
     A .6302432E+08, .6302432E+08, .6348947E+08, .6348947E+08,
     A .6146760E+08, .6146760E+08, .5817972E+08, .5817972E+08,
     A .4668909E+08, .4668909E+08, .1519062E+08, .1519062E+08,
     A .5133055E+07, .5133055E+07, .8998594E+06, .8998594E+06/
 
      data ((asy(i,j),j=1,10),i=1,4) /
     a 0.896,0.885,0.880,0.877,0.867,0.860,0.845,0.836,0.905,0.871,
     a 0.642,0.633,0.631,0.628,0.621,0.616,0.610,0.572,0.562,0.495,
     a 0.795,0.790,0.788,0.781,0.783,0.782,0.778,0.783,0.797,0.750,
     a 0.397,0.359,0.348,0.337,0.311,0.294,0.253,0.154,0.103,0.055/
 
      data ((asy2(i,j),j=1,10),i=1,1)/
     A 0.718,0.712,0.710,0.708,0.704,0.702,0.696,0.680,0.668,0.649/    
 
      data ((asy3(i,j),j=1,10),i=1,1)/
     A 0.704,0.690,0.686,0.680,0.667,0.659,0.637,0.541,0.437,0.241/    
 
      data ((asy4(i,j),j=1,10),i=1,1)/
     A .705, .744, .751, .757, .762, .759, .737, .586, .372, .139/
 
      common /sixs_aer/ ext(10),ome(10),gasym(10),phase(10)
      common /sixs_aerbas/ ph(10,83)
      common /sixs_sos/phasel(10,83),cgaus(83),pdgs(83)
c
c     optical properties of aerosol model computed from sra basic comp
      pi=4.*atan(1.)
      do 1 l=1,10
       ext(l)=0.
       sca(l)=0.
       if(l.eq.4.and.iaer.eq.0) ext(l)=1.
       ome(l)=0.
       gasym(l)=0.
       phase(l)=0.
       do 1 k=1,83
        phasel(l,k)=0.
    1 continue
 
      do 2 j=1,4
       ci(j)=co(j)
    2 continue
 
      if(iaer.eq.0) return

      do 7 k=1,82
      if((xmud.ge.cgaus(k)).and.(xmud.lt.cgaus(k+1))) go to 8
    7 continue
      return
    8 j1=k
      j2=j1+1
      coef=-(xmud-cgaus(j1))/(cgaus(j2)-cgaus(j1))

      if (iaer.eq.12) then
        open(10,file=FILE)
        read(10,*)
        do l=1,10
         read(10,'(8x,4(3x,f6.4,3x))')ext(l),sca(l),ome(l),gasym(l)
        enddo    
        read(10,'(///)')
        do k=1,83
         read(10,'(8x,10(1x,e10.4))')(phasel(l,k),l=1,10)
        enddo   
        close(10)
        do l=1,10
         phase(l)=phasel(l,j1)+coef*(phasel(l,j1)-phasel(l,j2))
        enddo
        return
      endif
c
      if (iaer.eq.5) then
        do k=1,10
        asy(1,k)=asy2(iaer-4,k)
        ex(1,k)=ex2(iaer-4,k)
        sc(1,k)=sc2(iaer-4,k)
        enddo
      endif
c
      if (iaer.eq.6) then
        do k=1,10
        asy(1,k)=asy3(iaer-5,k)
        ex(1,k)=ex3(iaer-5,k)
        sc(1,k)=sc3(iaer-5,k)
        enddo
      endif
c
      if (iaer.eq.7) then
        do k=1,10
        asy(1,k)=asy4(iaer-6,k)
        ex(1,k)=ex4(iaer-6,k)
        sc(1,k)=sc4(iaer-6,k)
        enddo
      endif
c
c
      if (iaer.ge.5.and.iaer.le.11) then
c calling a special aerosol model 
C     (background desert model...)
         if (iaer.eq.5) call bdm
C     (biomass burning model...)
         if (iaer.eq.6) call bbm
C     (stratospherique aerosol model...)
         if (iaer.eq.7) call stm
C     (user defined model from size distribution)
         if (iaer.ge.8.and.iaer.le.11) call mie(iaer,wldis,ex,sc,asy)

         do l=1,10
          dd(1,l)=ph(l,j1)+coef*(ph(l,j1)-ph(l,j2))
          do k=1,83
           pha(1,l,k)=ph(l,k)
          enddo
         enddo
         icp=1
         cij(1)=1.00
c for normalization of the extinction coefficient
         nis=1.d+00/ex(1,4)
      else
c calling each sra components
         icp=4
c  -dust
         call dust
         do l=1,10
         dd(1,l)=ph(l,j1)+coef*(ph(l,j1)-ph(l,j2))
         do k=1,83
         pha(1,l,k)=ph(l,k)
         enddo
         enddo
c  -water soluble
         call wate
         do l=1,10
         dd(2,l)=ph(l,j1)+coef*(ph(l,j1)-ph(l,j2))
         do k=1,83
         pha(2,l,k)=ph(l,k)
         enddo
         enddo
c  -oceanic type
         call ocea
         do l=1,10
         dd(3,l)=ph(l,j1)+coef*(ph(l,j1)-ph(l,j2))
         do  k=1,83
         pha(3,l,k)=ph(l,k)
         enddo
         enddo
c  - soot
         call soot
         do l=1,10
         dd(4,l)=ph(l,j1)+coef*(ph(l,j1)-ph(l,j2))
         do k=1,83
         pha(4,l,k)=ph(l,k)
         enddo
         enddo
c     summ of the ci/vi calculation
         sumni=0.
         sigm=0.
         do 3 i=1,4
    3    sigm=sigm+ci(i)/vi(i)
 
c     cij coefficients calculation
         do 4 j=1,4
         cij(j)=(ci(j)/vi(j)/sigm)
    4    sumni=sumni+cij(j)/ni(j)

c nis=1/Kext(550)
         nis=1.d+00/sumni
      endif
      
c     mixing parameters calculation
      do 5 l=1,10
      do 6 j=1,icp
      ext(l)=ex(j,l)*cij(j)+ext(l)
      sca(l)=sc(j,l)*cij(j)+sca(l)
      gasym(l)=sc(j,l)*cij(j)*asy(j,l)+gasym(l)
      phase(l)=sc(j,l)*cij(j)*dd(j,l)+phase(l)
      do 77 k=1,83
      phasel(l,k)=sc(j,l)*cij(j)*pha(j,l,k)+phasel(l,k)
   77 continue
    6 continue
      ome(l)=sca(l)/ext(l)
      gasym(l)=gasym(l)/sca(l)
      phase(l)=phase(l)/sca(l)
      do 78 k=1,83
      phasel(l,k)=phasel(l,k)/sca(l)
   78 continue
      ext(l)=ext(l)*nis
      sca(l)=sca(l)*nis
    5 continue
      if (iaer.ge.8.and.iaer.le.11) then
       open(10,file=FILE)
        write(10,'(3x,A5,1x,5(1x,A10,1x),1x,A10)')'Wlgth',
     s'Nor_Ext_Co','Nor_Sca_Co','Sg_Sca_Alb',
     s'Asymm_Para','Extinct_Co','Scatter_Co'
        do 79 l=1,10
         write(10,'(2x,f6.4,4(3x,f6.4,3x),2(2x,e10.4))')
     s wldis(l),ext(l),sca(l),ome(l),gasym(l),ext(l)/nis,sca(l)/nis
 79     continue
         write(10,'(//,T20,A16,/,3x,A4,1x,10(3x,f6.4,2x))')
     s   ' Phase Function ','TETA',(wldis(l),l=1,10)
        do 76 k=1,83
         write(10,'(2x,f6.2,10(1x,e10.4))')180.*acos(cgaus(k))/pi,
     s                 (phasel(l,k),l=1,10)
 76     continue
        close(10)
      endif
      return
      end
      subroutine msrm
c
c   MultiSpectral Reflectance Model 93         A.Kuusk   24.03.1993
c
      implicit double precision (a-h, o-z)
      save /count/, /soildata/, /aaa/, /ggg/, /ladak/
c
      dimension u1(10), u2(10), a1(10), a2(10)
      common /count/ jl, jj, lg, jg, lf, nnx, n1, n2, u1, u2, a1, a2
c
      double precision nnl, kk
      common /leafin/ nnl, vai, kk
      common /leafout/ refl, tran
c
      double precision ke, kab, kw
      dimension refr(200), ke(200), kab(200), kw(200)
      common /dat/ refr, ke, kab, kw
c
      dimension phis1(200), phis2(200), phis3(200), phis4(200)
      common /soildata/ phis1, phis2, phis3, phis4, rsl1, rsl2, rsl3,  
     & rsl4, th2, rsl, rsoil, rr1soil, rrsoil
c
      common /aaa/ rrl, ttl, ul, sl, clmp, clmp1, bi, bd, bqint
      common /ggg/ gr, gt, g, g1, th, sth, cth, th1, sth1, cth1, 
     & phi, sp, cp, th22, st, ct, st1, ct1, t10, t11, e1, e2,
     & s2, s3, ctg, ctg1, ctt1, stt1, calph, alp2, salp2, calp2, 
     & alph, salph, alpp, difmy, difsig
      common /cfresn/ rn, rk
      common /ladak/ ee, thm, sthm, cthm
      common /msrmdata/ th10, rncoef, cab, cw, bq
c
      data pi12/1.570796326794895d0/, pi/3.141592653589793d0/
      data eps4/.1d-3/
c
*           print *, 'msrm'
c
      sth10 = sin(th10)
      cth10 = cos(th10)
c
      sp    = sin(phi)
      cp    = cos(phi)
      th1   = th10
      sth1  = sth10
      cth1  = cth10
      sth   = sin(th)
      cth   = cos(th)
      rrls  = rrl
c
      call biz
c
      rrl  = refl
      rtp  = rrl + ttl
c
      call difr92
c
10    continue
c
      rrl = rrls
      bq  = bi + bd
c
      return
      end
*
******************************************************************
*
      subroutine akd
c  bdz   A.Kuusk    4.03.1988
c
      implicit double precision (a-h, o-z)
      save /count/, /aaa/, /ggg/
c
      dimension tt3(10), stt3(10), ctt3(10), tt2(10), stt2(10), ctt2(10)
c
      dimension u1(10), u2(10), a1(10), a2(10)
      common /count/ jl, jj, lg, jg, lf, nnx, n1, n2, u1, u2, a1, a2
c
      double precision nnl, kk
      common /leafin/ nnl, vai, kk
      common /leafout/ refl, tran
c
      common /aaa/ rrl, ttl, ul, sl, clmp, clmp1, bi, bd, bqint
      common /ggg/ gr, gt, g, g1, th, sth, cth, th1, sth1, cth1, 
     & phi, sp, cp, th22, st, ct, st1, ct1, t10, t11, e1, e2,
     & s2, s3, ctg, ctg1, ctt1, stt1, calph, alp2, salp2, calp2, 
     & alph, salph, alpp, difmy, difsig
c
      data pi/3.141592653589793d0/, pi1/1.5707963268d0/, eps/.005d0/
c
*                    print *, 'akd'
      bqint = 0.d0
      if (th .gt. eps) goto 4
      phi = 0.d0
      sp  = 0.d0
      cp  = 1.d0
c
      do 10 i2 = 1, n2
         th1  = (1.d0 - u2(i2))*pi1
         sth1 = sin(th1)
         cth1 = cos(th1)
         rrls = rrl
c
         call biz
c
         rrl = refl
         rtp = rrl + ttl
c
         call difr92
c
         rrl   = rrls
         bqint = bqint + a2(i2)*(bi + bd)*sth1*cth1
10    continue
c
      bqint = bqint*pi
      goto 1
c
4     continue
      do 14 i = 1, n1
         thi     = u1(i)*th
         tt3(i)  = thi
         stt3(i) = sin(thi)
         ctt3(i) = cos(thi)
14    continue
c
      do 15 i = 1, n2
         thi     = u2(i)*(th - pi1) + pi1
         tt2(i)  = thi
         stt2(i) = sin(thi)
         ctt2(i) = cos(thi)
15    continue
c
      do 11 j = 1, n1
         phi  = (1.d0 - u1(j))*pi
         sp   = sin(phi)
         cp   = cos(phi)
         bd1  = 0.d0
         bd2  = 0.d0
         do 12 i1 = 1, n1
            th1  = tt3(i1)
            sth1 = stt3(i1)
            cth1 = ctt3(i1)
c
         rrls = rrl
c
         call biz
c
         rrl = refl
         rtp = rrl + ttl
c
         call difr92
c
         rrl = rrls
c
            bd1 = bd1 + a1(i1)*(bi + bd)*sth1*cth1
12       continue
c
         do 13 i2 = 1, n2
            th1  = tt2(i2)
            sth1 = stt2(i2)
            cth1 = ctt2(i2)
c
         rrls = rrl
c
         call biz
c
         rrl = refl
         rtp = rrl + ttl
c
         call difr92
c
         rrl = rrls
c
            bd2 = bd2 + a2(i2)*(bi + bd)*sth1*cth1
13       continue
c
         bqint = bqint + ((pi1 - th)*bd2 + th*bd1)*a1(j)
11    continue
c
      bqint = bqint + bqint
c
1     return
      end
*
******************************************************************
*
      subroutine biz
c     canopy reflectance of single scattering for direct radiation
c     A. Kuusk   6.02.1992
c
      implicit double precision (a-h, o-z)
      double precision integr
      save /count/, /soildata/, /aaa/, /ggg/, /ladak/
c
*     dimension gj(2), g1j(2), grj(2), gtj(2), gfj(2)
c
      dimension u1(10), u2(10), a1(10), a2(10)
      common /count/ jl, jj, lg, jg, lf, nnx, n1, n2, u1, u2, a1, a2
c
      dimension phis1(200), phis2(200), phis3(200), phis4(200)
      common /soildata/ phis1, phis2, phis3, phis4, rsl1, rsl2, 
     & rsl3, rsl4, th2, rsl, rsoil, rr1soil, rrsoil
c
      common /aaa/ rrl, ttl, ul, sl, clmp, clmp1, bi, bd, bqint
      common /ggg/ gr, gt, g, g1, th, sth, cth, th1, sth1, cth1, 
     & phi, sp, cp, th22, st, ct, st1, ct1, t10, t11, e1, e2,
     & s2, s3, ctg, ctg1, ctt1, stt1, calph, alp2, salp2, calp2, 
     & alph, salph, alpp, difmy, difsig
      common /ladak/ ee, thm, sthm, cthm
c
      data pi/3.14159265358979d0/, eps/.1d-4/, eps3/.01d0/
c
      integr(xx) = (1.d0 - exp(-xx))/xx
*           print *, 'biz in'
      ths   = th
      sths  = sth
      cths  = cth
      th1s  = th1
      sth1s = sth1
      cth1s = cth1
*     thms  = thm
c
      call soil
c
      if (ul .gt. eps) goto 2
      bi  = rsoil
      goto 1
c
2     continue
      if (th1 .lt. th) goto 12
      t11  = th1
      st   = sth
      st1  = sth1
      ct   = cth
      ct1  = cth1
      t10  = th
      jj   = 0
      goto 7
c
12    t10  = th1
      st   = sth1
      st1  = sth
      ct   = cth1
      ct1  = cth
      t11  = th
      jj   = 1
c
7     continue
      ctt1  = ct*ct1
      stt1  = st*st1
      calph = stt1*cp + ctt1
      catmp = calph
      alph  = acos(catmp)
      alp2  = alph*.5d0
*     if (lf .ne. 2) then
*        if( jg .gt. 2) then
*           print *, ' ***  biz3:  jg > 2  ***'
*           stop
*        endif
         e1   = st*ct1
         e2   = ct*st1
         s2   = e1*cp + e2
         s3   = e1*sp
         ctg  = 1.d30
         ctg1 = 1.d30
         if (st .ne. 0.d0) ctg = ct/st
         if (st1 .ne. 0.d0) ctg1 = ct1/st1
         salph = sin(alph)
         alpp  = pi - alph
         salp2 = sin(alp2)
         calp2 = cos(alp2)
c
         call gmf(gf)
c
         if (ee .le. eps3) goto 95
         y4  = abs(cth + cth1)*.5d0/calp2
         if (y4.lt.1.d0) thp = acos(y4)
c
95       call glak(glthp, thp)
c
         x2 = glthp*.125d0
         gf = gf*x2
c
         call gmd92
c
      gammd = gr*rrl + gt*ttl
c
      t11 = th1
      st  = sth
      st1 = sth1
      ct  = cth
      ct1 = cth1
      t10 = th
      if (jj .eq. 1) then
         x = g1
         g1 = g
         g = x
      endif
c
*           print *, 'biz:2'
      gg   = g*g1
      g    = g*clmp
      g1   = g1*clmp1
      gg1  = g*ct1 + g1*ct
      sct  = sqrt(ctt1)
      alpd = alp2/sl
      bam  = alpd*sct/ul
c
      xx1  = 0.d0
      if (ctt1 .gt. eps) then
         gma  = alpd/sct
         ulg  = gg1/ctt1*ul
         ulg1 = ulg*.5d0
         xx1  = ulg + gma
      endif
      if ((xx1 .gt. 30.d0) .or. (ctt1 .le. eps)) then
         easte  = 0.d0
         easte2 = 0.d0
         easte4 = 0.d0
         bs1    = 0.d0
      else
         easte  = exp(-ulg)
         easte2 = exp(-ulg1 - gma)
         easte4 = exp(-ulg - gma)
         bs1    = (easte + easte2 - easte4)*rsoil
      endif
c
      xx1   = (1.d0 - easte)/gg1
      xx2   = (1.d0 - easte2)/(gg1*.5d0 + bam) - 
     & (1.d0 - easte4)/(gg1 + bam)
      bc1d  = xx1*gammd
      bc1hs = xx2*(gammd + gf)
      bcsp  = xx1*gf
      bc1   = bc1d + bcsp + bc1hs
      bi    = bc1 + bs1
c
1     continue
      th    = ths
      sth   = sths
      cth   = cths
      th1   = th1s
      sth1  = sth1s
      cth1s = cth1
*     thm   = thms
c
      return
      end
*
******************************************************************
*
      subroutine difr92
c   diffuse fluxes according to SAIL for an elliptical LAD
c   A. Kuusk 16.06.1992
c
      implicit double precision (a-h, o-z)
      double precision ks, ko, m, m11, m12, m21, m22, integr
      save /soildata/, /aaa/, /ggg/, /ladak/
c
      dimension phis1(200), phis2(200), phis3(200), phis4(200)
      common /soildata/ phis1, phis2, phis3, phis4, rsl1, rsl2, 
     & rsl3, rsl4, th2, rsl, rsoil, rr1soil, rrsoil
c
      common /aaa/ rrl, ttl, ul, sl, clmp, clmp1, bi, bd, bqint
      common /ggg/ gr, gt, g, g1, th, sth, cth, th1, sth1, cth1, 
     & phi, sp, cp, th22, st, ct, st1, ct1, t10, t11, e1, e2,
     & s2, s3, ctg, ctg1, ctt1, stt1, calph, alp2, salp2, calp2, 
     & alph, salph, alpp, difmy, difsig
      common /ladak/ ee, thm, sthm, cthm
c
      integr(x) = (1.d0 - exp(-x))/x
*           print *, 'difr92'
c
      tsun  = th1
      tview = th
      tants = sth1/cth1
      tanto = sth/cth
      rtp   = (rrl + ttl)/2.d0
c
      ks    = g1*ul/cth1
      ko    = g*ul/cth
      gg    = (1.289d0*difmy - 1.816d0*difsig)*(cthm**2 - 
     & .33333333333d0) + .31823d0
      bf    = (rrl - ttl)/2.d0*ul*gg
      att   = (1.d0 - rtp)*ul + bf
      sig   = rtp*ul + bf
      sb    = ks*rtp + bf
      sf    = ks*rtp - bf
      ub    = ko*rtp + bf
      uf    = ko*rtp - bf
      m     = sqrt(att**2 - sig**2)
      h1    = (att + m)/sig
      h2    = 1.d0/h1
      c     = (sf*sig - sb*(ks - att))/(m**2 - ks**2)
      d     = (sb*sig + sf*(ks + att))/(m**2 - ks**2)
*     epso  = skyl - d*sq
      epso  =  - d
*     epss  = (rrsoil*(d + 1.d0) - c)*sq*exp(-ks)
      epss  = (rrsoil*(d + 1.d0) - c)*exp(-ks)
      m11   = h1
      m12   = h2
      m21   = (1.d0 - rrsoil*h1)*exp(-m)
      m22   = (1.d0 - rrsoil*h2)*exp(m)
      det   = m11*m22 - m12*m21
      a     = (m22*epso - m12*epss)/det
      b     = (-m21*epso + m11*epss)/det
      ep    = integr(ko + m)
      em    = integr(ko - m)
      ek    = integr(ko + ks)
*     gp    = a*ep + b*em + c*ek*sq
      gp    = a*ep + b*em + c*ek
*     gm    = h1*a*ep + h2*b*em + d*ek*sq
      gm    = h1*a*ep + h2*b*em + d*ek
*     ems   = h1*a*exp(-m) + h2*b*exp(m) + d*sq*exp(-ks)
      ems   = h1*a*exp(-m) + h2*b*exp(m) + d*exp(-ks)
      rplants = uf*gp + ub*gm
      rdsoil  = rrsoil*ems*exp(-ko)
      bd    = rplants + rdsoil
c
      return
      end
*
**********************************************************************
*
      subroutine glak(glth, th)
c  elliptical distribution
c  A.Kuusk   1.03.1988
c
      implicit double precision (a-h, o-z)
      save /aaa/, /ladak/
      save bb, es, tms
c
      common /aaa/ rrl, ttl, ul, sl, clmp, clmp1, bi, bd, bqint
      common /ladak/ ee, thm, sthm, cthm
c
      data bb/1.d0/, es/0.d0/, tms/0.d0/, eps/.1d0/
c
*           print *, 'gl'
c
      if (ee .lt. eps) then
      glth = 1.d0
      return
      endif
c
      if (ee .eq. 1.d0) ee = .999999d0
      if ((ee .ne. es) .or. (thm .ne. tms)) then
        u1  = ee*cthm
        u3  = ee*sthm
        u2  = sqrt(1.d0 - u1*u1)
        u4  = sqrt(1.d0 - u3*u3)
        x   = log((u4 + u1)/(u2 - u3))
        x1  = atan2(u3, u4) - atan2(u1, u2)
        x2  = sthm*x - cthm*x1
        bb  = ee/x2
        es  = ee
        tms = thm
      endif
c
      glth = bb/sqrt(1.d0 - (ee*cos(thm - th))**2)
c
      return
      end
*
******************************************************************
*
      subroutine gmf(gf)
c  Fresnel' reflection                    A.Kuusk 02.01.1991
c  input parameters are ca = cos(th_incident),  rn=refract.ind., 
c  rk = leaf hair index
c
      implicit double precision (a-h, o-z)
      save /aaa/, /ggg/
c
      common /aaa/ rrl, ttl, ul, sl, clmp, clmp1, bi, bd, bqint
      common /ggg/ gr, gt, g, g1, th, sth, cth, th1, sth1, cth1, 
     & phi, sp, cp, th22, st, ct, st1, ct1, t10, t11, e1, e2,
     & s2, s3, ctg, ctg1, ctt1, stt1, calph, alp2, salp2, calp2, 
     & alph, salph, alpp, difmy, difsig
      common /cfresn/ rn, rk
c
      data pi12/1.570796326794895d0/
c
*           print *, 'gmf'
c
      ca=calp2
      x2  = ca*ca
      ag  = x2*2.d0 - 1.d0 + rn*rn
      bg  = 1.d0 + (ag - 2.d0)*x2
      xy  = ag - x2
      cg  = 2.d0*ca*sqrt(xy)
      sa2 = 1.d0 - x2
      y   = (bg + sa2*cg)*(ag + cg)
      y   = (ag - cg)*bg/y
      yy  = sqrt(sa2)/pi12/ca*rk
      gf  = exp(-yy)*y
c
      return
      end
*
******************************************************************
*
      subroutine soil
c   Soil directional reflectance and reflectance (albedo)
c   th, th1, th2 in radianes,  cp = cos(phi)
c   A.Kuusk     1.03.1988
c
      implicit double precision (a-h, o-z)
      save a, b, c, cts, ths1, ths2
      save /count/, /soildata/, /aaa/, /ggg/
c
      dimension phis1(200), phis2(200), phis3(200), phis4(200)
      common /soildata/ phis1, phis2, phis3, phis4, rsl1, rsl2, 
     & rsl3, rsl4, th2, rsl, rsoil, rr1soil, rrsoil
c
      common /count/ jl, jj, lg, jg, lf, nnx, n1, n2, u1, u2, a1, a2
      common /aaa/ rrl, ttl, ul, sl, clmp, clmp1, bi, bd, bqint
      common /ggg/ gr, gt, g, g1, th, sth, cth, th1, sth1, cth1, 
     & phi, sp, cp, th22, st, ct, st1, ct1, t10, t11, e1, e2,
     & s2, s3, ctg, ctg1, ctt1, stt1, calph, alp2, salp2, calp2, 
     & alph, salph, alpp, difmy, difsig
c
      data a/.45098d0/, b/5.7829d0/, c, cts/2*13.7575d0/
      data ths1, ths2/2*.785398163d0/
c
*           print *, 'soil'
      if (th2 .ne. ths2) then
         cts  = 16.41d0 - th2*th2*4.3d0
         ths2 = th2
      endif
      if (th1 .ne. ths1) then
         ths1 = th1
         x    = th1*th1
         a    = x*7.702d0 - 4.3d0
         b    = th1*7.363d0
         c    = 16.41d0 - x*4.3d0
      endif
      x2      = rsl/cts
      rsoil   = ((a*th + b*cp)*th + c)*x2
      rr1soil = (.7337d0*a + c)*x2
      rrsoil  = 14.25d0*x2
c
      return
      end
*
******************************************************************
*
      subroutine soilspec
c
c   Soil spectral reflectance,  Price,  RSE 33:113 - 121 (1990)
c
      implicit double precision (a-h, o-z)
      save /count/, /soildata/
c
      dimension u1(10), u2(10), a1(10), a2(10)
      common /count/ jl, jj, lg, jg, lf, nnx, n1, n2, u1, u2, a1, a2
c
      dimension phis1(200), phis2(200), phis3(200), phis4(200)
      common /soildata/phis1, phis2, phis3, phis4, rsl1, rsl2, 
     & rsl3, rsl4, th2, rsl, rsoil, rr1soil, rrsoil
c
      rsl = rsl1*phis1(jl) + rsl2*phis2(jl) + 
     &      rsl3*phis3(jl) + rsl4*phis4(jl)
c
      return
      end
*
**********************************************************************
*
      subroutine gmd92
c  phase function and G-funktion
c  A. Kuusk    22.03.1988 & 16.06.1992
c  0< = th,  th1, th2<=pi/2,  0<=phi<=pi
c
      implicit double precision (a-h, o-z)
      dimension f(5)
      save /aaa/, /ggg/, /ladak/
c
      common /aaa/ rrl, ttl, ul, sl, clmp, clmp1, bi, bd, bqint
      common /ggg/ gr, gt, g, g1, th, sth, cth, th1, sth1, cth1, 
     & phi, sp, cp, th22, st, ct, st1, ct1, t10, t11, e1, e2,
     & s2, s3, ctg, ctg1, ctt1, stt1, calph, alp2, salp2, calp2, 
     & alph, salph, alpp, difmy, difsig
      common /ladak/ ee, thm, sthm, cthm
c
      data pi/3.14159265358979d0/, pi4/6.28318531717958d0/, 
     & pi12/.159154943d0/, pi14/.636619773d0/, eps5/.1d-2/
     & , pi13/.1061032953d0/
c
*           print *, 'gmd92'
c
c                            ***  gammad,  e = 0.  ***
      gr0 = (salph + alpp*calph)*pi13
      gt0 = (salph - alph*calph)*pi13
      if (ee .lt. .4d0) then
         gr = gr0
         gt = gt0
         g  = .5d0
         g1 = .5d0
         return
      endif
c                            ***  gammad,  e = 1.  ***
      sg   = 0.d0
      sg1  = 0.d0
      sgmr = 0.d0
      sgmt = 0.d0
      if (th22 .lt. t11) goto 47
      assign 46 to l4
      goto 61
c
46    continue
      assign 48 to l4
      goto 64
c
47    continue
      if (th22 .lt. t10) goto 50
      assign 51 to l4
      goto 62
c
51    continue
      assign 46 to l4
      goto 65
c
50    continue
      assign 52 to l4
      goto 63
c
52    continue
      assign 48 to l4
      goto 65
c
48    continue
c
      gr1 = sgmr*pi12
      gt1 = sgmt*pi12
      gr  = gr0 - .0102d0 + 
     &    (1.742d0*difmy - .4557d0*difsig)*(gr1 - gr0)
      gt  = gt0 + .00653d0 + 
     &    (.2693d0*difmy + 5.821d0*difsig)*(gt1 - gt0)
      g   = (2.653d0*difmy + 1.432d0*difsig)*(sg  - .5d0) + .50072d0
      g1  = (2.653d0*difmy + 1.432d0*difsig)*(sg1 - .5d0) + .50072d0
c
49    continue
      return
c
c  ******************************    tl1 = 0.,  tl2=pi/2 - th1
c
61    assign 71 to l2
      goto 130
71    y = pp
      if (y .gt. 0.d0) sgmr = sgmr + y
      if (y .lt. 0.d0) sgmt = sgmt - y
      y1  = ct1*cthm
      sg1 = sg1 + abs(y1)
      goto l4
c
c  ******************************    tl1 = pi/2 - th1,  tl2=pi/2 - th
c
62    continue
      x2 = cthm/sthm
      x  = -ctg1*x2
      x1 = sqrt(1.d0 - x*x)
      fa = atan2(x1, x)
      fb = pi4 - fa
      assign 72 to l2
      goto 30
c
72    continue
      y = pp
      if (y .gt. 0.d0) sgmr = sgmr + y
      if (y .lt. 0.d0) sgmt = sgmt - y
      assign 73 to l2
      goto 130
c
73    y = pp - y
      if (y .gt. 0.d0) sgmr = sgmr + y
      if (y .lt. 0.d0) sgmt = sgmt - y
      goto l4
c
c  ******************************   tl1 = pi/2 - th,  tl2=pi/2
c
63    continue
      x2 = cthm/sthm
      x  = -ctg1*x2
      x1 = sqrt(1.d0 - x*x)
      fa = atan2(x1, x)
      f(2) = fa
      f(3) = pi4 - fa
      x  = -ctg*x2
      x1 = sqrt(1.d0 - x*x)
      fa = atan2(x1, x)
      fb = phi - fa
      if (fb .lt. 0.d0) fb = fb + pi4
      f(4) = fb
      f(5) = phi + fa
      do 75 ii = 2, 4
         i1 = ii + 1
         do 75 j = i1, 5
            fa = f(ii)
            fb = f(j)
            if (fb .gt. fa) goto 75
            f(ii) = fb
            f(j)  = fa
75    continue
      f(1) = f(5) - pi4
      i1   = 1
76    ii   = i1
      i1   = ii + 1
      fa   = f(ii)
      fb   = f(i1)
      assign 74 to l2
      goto 30
c
c  ******************************   tl1 = pi/2 - th,  tl2=pi/2
c
74    continue
      y = pp
      if (y .gt. 0.d0) sgmr = sgmr + y
      if (y .lt. 0.d0) sgmt = sgmt - y
      if (i1 .le. 4) goto 76
c
      x2 = ct*cthm
      x1 = st*sthm/x2
      x1 = sqrt(x1*x1 - 1.d0)
      x  = atan2(1.d0, x1)
      x  = (x + x1)*x2
      y  = x*pi14
      sg = sg + abs(y)
      goto l4
c
c  ******************************    tl1 = 0,  tl2=pi/2 - th
c
64    y1 = ct*cthm
      sg = sg + abs(y1)
      goto l4
c
c  ******************************    tl1 = pi/2 - th1,  tl2=pi/2
c
65    continue
      x2  = ct1*cthm
      x1  = st1*sthm/x2
      x1  = sqrt(x1*x1 - 1.d0)
      x   = atan2(1.d0, x1)
      x   = (x + x1)*x2
      y   = x*pi14
      sg1 = sg1 + abs(y)
      goto l4
c
c  ******************************    p(fa, fb)
c
30    x  = fb - fa
      if (x .gt. eps5) goto 31
      pp = 0.d0
      goto l2
31    if ((pi4 - x) .lt. eps5) goto 130
      sfa = sin(fa)
      sfb = sin(fb)
      cfa = cos(fa)
      cfb = cos(fb)
      pp  = x*ctt1*cthm*cthm
      y1  = x + sfb*cfb - sfa*cfa
      x   = cfa - cfb
      y1  = y1*cp + sp*x*(cfa + cfb)
      pp  = pp + stt1*.5d0*y1*sthm*sthm
      y1  = s2*(sfb - sfa) + s3*x
      pp  = pp + y1*sthm*cthm
      goto l2
c
130   x  = sthm*sthm
      pp = calph*x + ctt1*(2.d0 - 3.d0*x)
      pp = pp*pi
      goto l2
c
      end
*
******************************************************************
*
*
c     ******************************************************************
c     leaf reflectance and transmittance. 
c     Input data are refractive index n,  a structure parameter N 
c     and an absorption coefficient k:
c     the PROSPECT model,  Jacquemoud & Baret,  RSE 34:75-91 (1990)
c     ******************************************************************

      subroutine leaf
c
      implicit double precision (a-h, o-z)
c
      double precision nn, k, inex
      common /leafin/ nn, vai, k
      common /leafout/ refl, tran
      common /nagout/ inex
      common /tauin/ teta, ref
      common /tauout/ tau

c     ******************************************************************
c     determination of elementary reflectances et transmittances
c     ******************************************************************
c     ALLEN et al.,  1969,  Interaction of isotropic ligth with a compact
c     plant leaf,  J. Opt. Soc. Am.,  Vol.59,  10:1376-1379
c     JACQUEMOUD S. and BARET F.,  1990,  Prospect : a model of leaf
c     optical properties spectra,  Remote Sens. Environ.,  34:75-91
c     ******************************************************************

*                     print *, 'leaf'
      if (k .le. 0.d0) then
         k = 1.d0
      else
         call s13aaf
         k = (1.d0 - k)*exp(-k) + k**2*inex
      endif

      teta = 90.d0
      ref  = nn
c
      call tav
c
      t1   = tau
      teta = 59.d0
c
      call tav
c
      t2 = tau
      x1 = 1.d0 - t1
      x2 = t1**2*k**2*(nn**2 - t1)
      x3 = t1**2*k*nn**2
      x4 = nn**4 - k**2*(nn**2 - t1)**2
      x5 = t2/t1
      x6 = x5*(t1 - 1.d0) + 1.d0 - t2
      r  = x1 + x2/x4
      t  = x3/x4
      ra = x5*r + x6
      ta = x5*t

c     ******************************************************************
c     reflectances et transmittances corresponding to N elementary
c     layers
c     ******************************************************************
c     STOKES G.G.,  1862,  On the intensity of the light reflected from or
c     transmitted through a pile of plates,  Proceedings of the Royal
c     Society of London,  Vol.11,  545-556
c     ******************************************************************

      delta = (t**2 - r**2 - 1.d0)**2 - 4.d0*r**2
      alfa  = (1.d0 + r**2 - t**2 + sqrt(delta))/(2.d0*r)
      beta  = (1.d0 + r**2 - t**2 - sqrt(delta))/(2.d0*r)
      va    = (1.d0 + r**2 - t**2 + sqrt(delta))/(2.d0*r)
      vb    = sqrt(beta*(alfa - r)/(alfa*(beta - r)))
      s1    = ra*(va*vb**(vai - 1.d0) - 
     & va**(-1.d0)*vb**(-(vai - 1.d0))) + 
     & (ta*t - ra*r)*(vb**(vai - 1.d0) - vb**(-(vai - 1.d0)))
      s2    = ta*(va - va**(-1.d0))
      s3    = va*vb**(vai - 1.d0) - va**(-1.d0)*vb**(-(vai - 1.d0))
     &        - r*(vb**(vai - 1.d0) - vb**(-(vai - 1.d0)))
      refl  = s1/s3
      tran  = s2/s3
c
      return
      end


c     ******************************************************************
c     exponential integral: int(exp(-t)/t, t = x..inf)
c     ******************************************************************

      subroutine s13aaf
c
      implicit double precision (a-h, o-z)
c
      double precision nn, k, inex
      common /leafin/ nn, vai, k
      common /nagout/ inex
*                     print *, 's13aafin'

      if (k .gt. 4.d0) goto 10

      x  =  0.5d0 * k  -  1.d0
      y  =  (((((((((((((((-3.60311230482612224d-13
     &    *x + 3.46348526554087424d-12)*x - 2.99627399604128973d-11)
     &    *x + 2.57747807106988589d-10)*x - 2.09330568435488303d-9)
     &    *x + 1.59501329936987818d-8)*x - 1.13717900285428895d-7)
     &    *x + 7.55292885309152956d-7)*x - 4.64980751480619431d-6)
     &    *x + 2.63830365675408129d-5)*x - 1.37089870978830576d-4)
     &    *x + 6.47686503728103400d-4)*x - 2.76060141343627983d-3)
     &    *x + 1.05306034687449505d-2)*x - 3.57191348753631956d-2)
     &    *x + 1.07774527938978692d-1)*x - 2.96997075145080963d-1
      y  =  (y*x + 8.64664716763387311d-1)*x  +  7.42047691268006429d-1
      inex  =  y  -  log(k)
      goto 30

10    if (k .ge. 85.d0) go to 20
      x  =  14.5d0 / (k + 3.25d0)  -  1.d0
      y  =  (((((((((((((((-1.62806570868460749d-12
     &    *x - 8.95400579318284288d-13)*x - 4.08352702838151578d-12)
     &    *x - 1.45132988248537498d-11)*x - 8.35086918940757852d-11)
     &    *x - 2.13638678953766289d-10)*x - 1.10302431467069770d-9)
     &    *x - 3.67128915633455484d-9)*x - 1.66980544304104726d-8)
     &    *x - 6.11774386401295125d-8)*x - 2.70306163610271497d-7)
     &    *x - 1.05565006992891261d-6)*x - 4.72090467203711484d-6)
     &    *x - 1.95076375089955937d-5)*x - 9.16450482931221453d-5)
     &    *x - 4.05892130452128677d-4)*x - 2.14213055000334718d-3
      y  =  ((y*x - 1.06374875116569657d-2)*x - 
     &     8.50699154984571871d-2)*x  + 
     &     9.23755307807784058d-1
      inex  =  exp(-k) * y / k
      goto 30

20    inex  =  0.d0
      goto 30

30    continue
*                     print *, 's13aafout'
      return
      end

c     ******************************************************************
c     determination of tav for any solid angle
c     ******************************************************************
c     STERN F.,  1964,  Transmission of isotropic radiation across an
c     interface between two dielectrics,  Appl.Opt.,  Vol.3,  1:111-113
c     ALLEN W.A.,  1973,  Transmission of isotropic light across a
c     dielectric surface in two and three dimensions,  J.Opt.Soc.Am., 
c     Vol.63,  6:664-666
c     ******************************************************************

      subroutine tav
c
      implicit double precision (a-h, o-z)
      double precision k
c
      common /tauin/ teta, ref
      common /tauout/ tau
c
      data dr/1.745329251994330d-2/, eps/.1d-6/, 
     &     pi12/1.570796326794895d0/

*                     print *, 'tavin'
      teta = teta*dr
      r2   = ref**2
      rp   = r2 + 1.d0
      rm   = r2 - 1.d0
      a    = (ref + 1.d0)**2/2.d0
      k    = -(r2 - 1.d0)**2/4.d0
      ds   = sin(teta)

      if (abs(teta) .le. eps) then
         tau = 4.d0*ref/(ref + 1.d0)**2
      else

         if (abs(teta - pi12) .le. eps) then
            b1 = 0.d0
         else
            xxx = (ds**2 - rp/2.d0)**2 + k
            b1 = sqrt(xxx)
         endif

         b2 = ds**2 - rp/2.d0
         b  = b1 - b2
         ts = (k**2/(6.d0*b**3) + k/b - b/2.d0) - 
     &        (k**2/(6.d0*a**3) + k/a - a/2.d0)
         tp1 = -2.d0*r2*(b - a)/rp**2
         tp2 = -2.d0*r2*rp*log(b/a)/rm**2
         tp3 = r2*(1.d0/b - 1.d0/a)/2.d0
         tp4 = 16.d0*r2**2*(r2**2 + 1.d0)*dlog((2.d0*rp*b - rm**2)/
     &    (2.d0*rp*a - rm**2))/(rp**3*rm**2)
         tp5 = 16.d0*r2**3*(1.d0/(2.d0*rp*b - rm**2) - 1.d0/
     &    (2.d0*rp*a - rm**2))/rp**3
         tp  = tp1 + tp2 + tp3 + tp4 + tp5
         tau = (ts + tp)/(2.d0*ds**2)
      endif
*                     print *, 'tavout'
      return
      end
*
******************************************************************
*
c     constant values: refractive index (ref), albino and dry leaf
c     absorption (ke), chlorophyll a+b specific absorption coefficient
c     (kab), water specific absorption coefficient (kw), 
*     and basis functions for soil spectral reflectance phis1, phis2,
*     phis3 and phis4 (Price, 1990)
c     ******************************************************************
c     JACQUEMOUD S. AND BARET F., 1990, Prospect : a model of leaf
c     optical properties spectra, Remote Sens. Environ., 34:75-91
c     JACQUEMOUD S. et al., 1991, Validation d'un modele de reflectance
c     spectrale et directionnnelle de sol, 5ieme Colloque International
c     Mesures Physiques et Signatures en Teledetection, Courchevel
c     (France), 14-18 Janvier 1991
c     ******************************************************************

      block data valeur
c
      implicit double precision (a-h, o-z)
c
      double precision ke, kab, kw
      dimension ref(200), ke(200), kab(200), kw(200)
      common /dat/ ref, ke, kab, kw
c
      dimension phis1(200), phis2(200), phis3(200), phis4(200)
      common /soildata/ phis1, phis2, phis3, phis4, rsl1, rsl2,
     & rsl3, rsl4, th2, rsl, rsoil, rr1soil, rrsoil
c
      data (ref(i), i = 1, 100)/
     & 1.5123,1.5094,1.5070,1.5050,1.5032,1.5019,1.5007,1.4997,1.4988,
     & 1.4980,1.4969,
     & 1.4959,1.4951,1.4943,1.4937,1.4930,1.4925,1.4920,1.4915,1.4910,
     & 1.4904,1.4899,1.4893,1.4887,1.4880,1.4873,1.4865,1.4856,1.4846,
     & 1.4836,1.4825,1.4813,1.4801,1.4788,1.4774,1.4761,1.4746,1.4732,
     & 1.4717,1.4701,1.4685,1.4670,1.4654,1.4639,1.4624,1.4609,1.4595,
     & 1.4582,1.4570,1.4559,1.4548,1.4538,1.4528,1.4519,1.4510,1.4502,
     & 1.4495,1.4489,1.4484,1.4480,1.4477,1.4474,1.4472,1.4470,1.4468,
     & 1.4467,1.4465,1.4463,1.4461,1.4458,1.4456,1.4453,1.4450,1.4447,
     & 1.4444,1.4440,1.4435,1.4430,1.4423,1.4417,1.4409,1.4402,1.4394,
     & 1.4387,1.4380,1.4374,1.4368,1.4363,1.4357,1.4352,1.4348,1.4345,
     & 1.4342,1.4341,1.4340,1.4340,1.4341,1.4342,1.4343,1.4345/

      data (ref(i), i = 101, 200)/
     & 1.4347,1.4348,1.4347,1.4345,1.4341,1.4336,1.4331,1.4324,1.4317,
     & 1.4308,1.4297,1.4284,1.4269,1.4253,1.4235,1.4216,1.4196,1.4176,
     & 1.4156,1.4137,1.4118,1.4100,1.4082,1.4065,1.4047,1.4029,1.4011,
     & 1.3993,1.3975,1.3958,1.3940,1.3923,1.3906,1.3888,1.3870,1.3851,
     & 1.3830,1.3808,1.3784,1.3758,1.3731,1.3703,1.3676,1.3648,1.3620,
     & 1.3592,1.3565,1.3537,1.3510,1.3484,1.3458,1.3433,1.3410,1.3388,
     & 1.3368,1.3350,1.3333,1.3317,1.3303,1.3289,1.3275,1.3263,1.3251,
     & 1.3239,1.3228,1.3217,1.3205,1.3194,1.3182,1.3169,1.3155,1.3140,
     & 1.3123,1.3105,1.3086,1.3066,1.3046,1.3026,1.3005,1.2985,1.2964,
     & 1.2944,1.2923,1.2902,1.2882,1.2863,1.2844,1.2826,1.2808,1.2793,
     & 1.2781,1.2765,1.2750,1.2738,1.2728,1.2719,1.2712,1.2708,1.2712,
     & 1.2736/

      data (ke(i), i = 1, 100)/
     &.1104,.0893,.0714,.0567,.0442,.0348,.0279,.0232,.0197,.0173,.0154,
     &.0142,.0120,.0108,.0093,.0092,.0092,.0092,.0092,.0092,.0091,.0091,
     &.0091,.0091,.0091,.0090,.0090,.0090,.0090,.0090,.0089,.0089,.0089,
     &.0089,.0088,.0088,.0088,.0088,.0088,.0087,.0087,.0087,.0087,.0087,
     &.0086,.0086,.0086,.0086,.0086,.0085,.0085,.0085,.0085,.0085,.0084,
     &.0084,.0084,.0084,.0084,.0083,.0083,.0083,.0082,.0082,.0082,.0082,
     &.0082,.0081,.0081,.0081,.0081,.0081,.0080,.0080,.0080,.0080,.0080,
     &.0079,.0079,.0079,.0079,.0079,.0078,.0078,.0078,.0078,.0078,.0077,
     &.0077,.0077,.0077,.0077,.0076,.0076,.0076,.0076,.0076,.0075,.0075,
     &.0075/

      data (ke(i), i = 101, 200)/
     &.0074,.0073,.0072,.0071,.0070,.0069,.0068,.0068,.0067,.0066,.0065,
     &.0064,.0063,.0062,.0062,.0061,.0060,.0059,.0058,.0057,.0056,.0056,
     &.0054,.0053,.0053,.0052,.0051,.0050,.0049,.0048,.0047,.0047,.0046,
     &.0045,.0044,.0043,.0042,.0041,.0040,.0039,.0039,.0037,.0037,.0036,
     &.0035,.0034,.0033,.0032,.0031,.0031,.0030,.0029,.0028,.0027,.0026,
     &.0025,.0025,.0024,.0023,.0022,.0021,.0020,.0019,.0019,.0018,.0017,
     &.0016,.0015,.0014,.0014,.0013,.0012,.0010,.0010,.0009,.0008,.0007,
     &.0006,.0006,.0005,.0004,.0003,.0002,.0002,.0001,15*.0000/

      data kab/
     & .04664,.04684,.04568,.04482,.04344,.04257,.04287,.04189,.04116,
     & .03847,.03409,
     & .03213,.03096,.03116,.03051,.03061,.02998,.02965,.02913,.02902,
     & .02769,.02707,.02539,.02409,.02150,.01807,.01566,.01317,.01095,
     & .00929,.00849,.00803,.00788,.00757,.00734,.00713,.00692,.00693,
     & .00716,.00758,.00815,.00877,.00938,.00976,.01041,.01089,.01105,
     & .01127,.01170,.01222,.01280,.01374,.01441,.01462,.01495,.01499,
     & .01506,.01580,.01686,.01810,.01961,.02112,.02336,.02702,.02880,
     & .02992,.03142,.03171,.02961,.02621,.02078,.01518,.01020,.00718,
     & .00519,.00390,.00298,.00218,.00163,.00116,.00083,.00057,.00039,
     & .00027,.00014,.00011,.00009,.00005,112*.00000/

      data kw/
     & 111*0.,00.100,00.200,00.278,00.206,00.253,00.260,00.313,00.285,
     & 00.653,00.614,00.769,00.901,00.872,00.812,00.733,00.724,00.855,
     & 00.900,01.028,01.500,02.026,02.334,03.636,08.942,14.880,17.838,
     & 19.497,19.419,17.999,12.024,10.709,08.384,07.081,06.155,05.619,
     & 05.112,04.512,04.313,04.064,03.804,03.709,03.877,04.348,04.574,
     & 05.029,05.804,06.345,05.823,05.886,06.315,08.432,15.588,32.247,
     & 51.050,58.694,55.135,50.454,42.433,40.670,36.030,29.771,25.153,
     & 24.378,22.008,20.608,18.576,17.257,15.921,14.864,12.861,12.773,
     & 12.426,13.090,14.013,15.066,15.857,16.776,19.113,21.066,22.125,
     & 26.438,28.391,28.920,31.754,36.375,40.056,41.019,45.471,43.126/

      data (phis1(i), i = 1, 100)/
     &  .088, .095, .102, .109, .116, .123, .130, .136, .143, .150,
     &  .157, .164, .171, .178, .185, .192, .199, .206, .213, .220,
     &  .227, .233, .240, .247, .254, .261, .268, .275, .282, .289,
     &  .295, .302, .309, .316, .326, .335, .345, .356, .366, .376,
     &  .386, .395, .404, .412, .421, .429, .436, .443, .450, .457,
     &  .464, .470, .476, .483, .489, .495, .502, .508, .514, .520,
     &  .526, .532, .538, .543, .549, .555, .561, .568, .574, .580,
     &  .587, .594, .601, .608, .615, .622, .629, .637, .644, .652,
     &  .659, .667, .674, .681, .689, .696, .702, .709, .716, .723,
     &  .729, .735, .742, .748, .754, .760, .766, .771, .777, .782/

      data (phis1(i), i = 101, 200)/
     &  .802, .819, .832, .842, .854, .868, .883, .899, .917, .935,
     &  .954, .974, .993,1.012,1.030,1.047,1.063,1.078,1.091,1.102,
     & 1.111,1.118,1.126,1.137,1.150,1.163,1.176,1.187,1.192,1.188,
     & 1.177,1.159,1.134,1.090, .979, .830, .764, .744, .748, .777,
     &  .823, .878, .932, .983,1.026,1.062,1.091,1.115,1.133,1.147,
     & 1.156,1.161,1.162,1.158,1.149,1.132,1.109,1.087,1.072,1.056,
     & 1.035, .989, .886, .659, .456, .350, .323, .335, .361, .396,
     &  .438, .484, .530, .576, .622, .664, .705, .740, .768, .788,
     &  .800, .802, .796, .794, .797, .789, .779, .756, .725, .715,
     &  .675, .635, .585, .535, .485, .435, .385, .335, .285, .235/

      data (phis2(i), i = 1, 100)/
     &  .249, .245, .241, .237, .232, .228, .222, .217, .211, .205,
     &  .199, .193, .186, .179, .171, .163, .155, .147, .139, .130,
     &  .121, .111, .102, .092, .081, .071, .060, .049, .038, .026,
     &  .014, .002,-.011,-.024,-.037,-.050,-.064,-.078,-.092,-.107,
     & -.121,-.137,-.152,-.168,-.184,-.200,-.216,-.232,-.246,-.259,
     & -.270,-.280,-.289,-.297,-.303,-.308,-.313,-.317,-.322,-.325,
     & -.329,-.332,-.335,-.338,-.340,-.342,-.345,-.347,-.350,-.352,
     & -.355,-.358,-.360,-.363,-.366,-.369,-.372,-.374,-.377,-.378,
     & -.380,-.381,-.382,-.382,-.383,-.382,-.382,-.381,-.380,-.378,
     & -.376,-.373,-.370,-.367,-.363,-.359,-.354,-.349,-.344,-.338/

      data (phis2(i), i = 101, 200)/
     & -.310,-.283,-.258,-.234,-.212,-.190,-.167,-.143,-.118,-.092,
     & -.066,-.039,-.014, .011, .034, .057, .083, .114, .151, .192,
     &  .233, .272, .311, .348, .380, .407, .438, .476, .521, .570,
     &  .624, .674, .708, .766, .824, .853, .854, .852, .858, .881,
     &  .916, .947, .973, .997,1.017,1.036,1.052,1.067,1.082,1.095,
     & 1.107,1.119,1.131,1.142,1.154,1.166,1.175,1.179,1.178,1.172,
     & 1.162,1.148,1.083, .900, .678, .538, .499, .515, .552, .598,
     &  .653, .716, .777, .834, .886, .932, .973,1.007,1.036,1.058,
     & 1.075,1.086,1.091,1.091,1.086,1.076,1.060,1.039,1.012, .980,
     &  .943, .900, .852, .799, .740, .676, .606, .532, .451, .366/

      data (phis3(i), i = 1, 100)/
     & -.417,-.384,-.351,-.318,-.285,-.253,-.221,-.189,-.157,-.126,
     & -.095,-.064,-.033,-.003, .027, .057, .087, .117, .146, .175,
     &  .204, .232, .260, .289, .316, .344, .371, .399, .425, .452,
     &  .478, .505, .525, .545, .566, .587, .606, .626, .652, .676,
     &  .699, .722, .744, .764, .784, .804, .822, .839, .856, .872,
     &  .886, .900, .913, .926, .937, .948, .957, .966, .974, .981,
     &  .988, .993, .998,1.002,1.006,1.009,1.012,1.014,1.016,1.017,
     & 1.018,1.018,1.018,1.017,1.016,1.014,1.012,1.010,1.007,1.003,
     &  .999, .995, .990, .984, .978, .972, .965, .957, .949, .941,
     &  .932, .923, .913, .902, .891, .880, .868, .855, .842, .829/

      data (phis3(i), i = 101, 200)/
     &  .766, .694, .620, .550, .484, .421, .361, .303, .247, .190,
     &  .134, .079, .023,-.031,-.086,-.140,-.190,-.235,-.275,-.310,
     & -.340,-.367,-.394,-.422,-.452,-.484,-.513,-.541,-.565,-.578,
     & -.575,-.556,-.525,-.468,-.323,-.115,-.018, .002,-.003,-.029,
     & -.076,-.142,-.211,-.274,-.333,-.386,-.432,-.471,-.503,-.528,
     & -.544,-.551,-.549,-.538,-.517,-.491,-.463,-.436,-.419,-.417,
     & -.401,-.348,-.216, .014, .160, .203, .209, .210, .207, .200,
     &  .189, .174, .155, .132, .105, .075, .043, .013,-.012,-.035,
     & -.053,-.068,-.078,-.082,-.080,-.073,-.060,-.041,-.017, .006,
     &  .035, .065, .097, .125, .168, .180, .168, .125, .097, .065/

      data (phis4(i), i = 1, 100)/
     &  .067, .077, .086, .094, .102, .111, .118, .126, .133, .140,
     &  .146, .152, .158, .164, .169, .174, .179, .184, .188, .192,
     &  .195, .198, .201, .204, .206, .208, .210, .212, .213, .214,
     &  .214, .214, .214, .214, .213, .212, .211, .210, .210, .209,
     &  .207, .205, .202, .198, .194, .189, .184, .179, .173, .167,
     &  .161, .155, .149, .143, .136, .130, .123, .116, .108, .101,
     &  .093, .085, .077, .068, .060, .051, .043, .034, .026, .018,
     &  .010, .002,-.006,-.014,-.022,-.030,-.037,-.045,-.052,-.060,
     & -.067,-.074,-.081,-.087,-.093,-.098,-.103,-.108,-.112,-.116,
     & -.120,-.123,-.126,-.129,-.132,-.134,-.136,-.138,-.140,-.141/

      data (phis4(i), i = 101, 200)/
     & -.147,-.152,-.158,-.166,-.170,-.165,-.157,-.151,-.144,-.128,
     & -.104,-.078,-.049,-.009, .038, .082, .122, .169, .222, .272,
     &  .317, .364, .413, .469, .532, .591, .642, .694, .748, .790,
     &  .810, .817, .819, .740, .494, .215, .110, .125, .155, .204,
     &  .291, .408, .521, .627, .724, .811, .884, .940, .987,1.025,
     & 1.053,1.071,1.077,1.072,1.046, .996, .941, .892, .857, .842,
     &  .809, .713, .509, .055,-.236,-.324,-.336,-.320,-.308,-.294,
     & -.275,-.248,-.205,-.144,-.094,-.048, .005, .058, .105, .132,
     &  .123, .079, .045, .024, .014, .018, .022,-.010,-.042,-.054,
     & -.055,-.060,-.060,-.055,-.050,-.046,-.042,-.038,-.034,-.030/

       end
*
******************************************************************
*
      subroutine dakg(u, a, nq)
c Gaussi kvadratuuri sqlmed ja kordajad, nq = 2*n, u=(-1., 1.)
      implicit double precision (a-h, o-z)
      dimension u(48), a(48)
c
*              print *,'dakg'
      n = nq/2
      goto (1, 2, 1, 4, 1, 6, 1, 8, 1, 10, 1, 12, 1, 14, 1, 16, 1, 
     & 1, 1, 20, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
     & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 48), nq
1     continue
      print *,  ' ***   dakg - inacceptable nq'
      stop 'dakg'
c
2     continue
      u(2) = .577350269189626d0
      a(2) = 1.d0
      goto 13
c
4     continue
      u(3) = .339981043584856d0
      u(4) = .861136311594053d0
      a(3) = .652145154862546d0
      a(4) = .347854845137454d0
      goto 13
c
6     continue
      u(4) = .238619186083197d0
      u(5) = .661209386466265d0
      u(6) = .932469514203152d0
      a(4) = .467913934572691d0
      a(5) = .360761573048139d0
      a(6) = .171324492379170d0
      goto 13
c
8     continue
      u(5) = .183434642495650d0
      u(6) = .525532409916329d0
      u(7) = .796666477413627d0
      u(8) = .960289856497536d0
      a(5) = .362683783378362d0
      a(6) = .313706645877887d0
      a(7) = .222381034453374d0
      a(8) = .101228536290376d0
      goto 13
c
10    continue
      u(6)  = .148874338981631d0
      u(7)  = .433395394129247d0
      u(8)  = .679409568299024d0
      u(9)  = .865063366688985d0
      u(10) = .973906528517172d0
      a(6)  = .295524224714753d0
      a(7)  = .269266719309996d0
      a(8)  = .219086362515982d0
      a(9)  = .149451349150580d0
      a(10) = .666713443086881d-1
      goto 13
c
12    continue
      u(7)  = .125233408511469d0
      u(8)  = .367831498998180d0
      u(9)  = .587317954286617d0
      u(10) = .769902674194305d0
      u(11) = .904117256370475d0
      u(12) = .981560634246719d0
      a(7)  = .249147045813402d0
      a(8)  = .233492536538355d0
      a(9)  = .203167426723066d0
      a(10) = .160078328543346d0
      a(11) = .106939325995318d0
      a(12) = .471753363865118d-1
      goto 13
c
14    continue
      u( 8) = .108054948707344d0
      u( 9) = .319112368927890d0
      u(10) = .515248636358154d0
      u(11) = .687292904811685d0
      u(12) = .827201315069765d0
      u(13) = .928434883663574d0
      u(14) = .986283808696812d0
      a( 8) = .215263853463158d0
      a( 9) = .205198463721296d0
      a(10) = .185538397477938d0
      a(11) = .157203167158194d0
      a(12) = .121518570687903d0
      a(13) = .801580871597602d-1
      a(14) = .351194603317519d-1
      goto 13
c
16    continue
      u( 9) = .950125098376374d-1
      u(10) = .281603550779259d0
      u(11) = .458016777657227d0
      u(12) = .617876244402643d0
      u(13) = .755404408355003d0
      u(14) = .865631202387832d0
      u(15) = .944575023073233d0
      u(16) = .989400934991650d0
      a( 9) = .189450610455068d0
      a(10) = .182603415044924d0
      a(11) = .169156519395003d0
      a(12) = .149595988816577d0
      a(13) = .124628971255534d0
      a(14) = .951585116824928d-1
      a(15) = .622535239386479d-1
      a(16) = .271524594117541d-1
      goto 13
c
20    continue
      u(11) = .765265211334973d-1
      u(12) = .227785851141645d0
      u(13) = .373706088715420d0
      u(14) = .510867001950827d0
      u(15) = .636053680726515d0
      u(16) = .746331906460151d0
      u(17) = .839116971822219d0
      u(18) = .912234428251326d0
      u(19) = .963971927277914d0
      u(20) = .993128599185095d0
      a(11) = .152753387130726d0
      a(12) = .149172986472604d0
      a(13) = .142096109318382d0
      a(14) = .131688638449177d0
      a(15) = .118194531961518d0
      a(16) = .101930119817240d0
      a(17) = .832767415767047d-1
      a(18) = .626720483341091d-1
      a(19) = .406014298003869d-1
      a(20) = .176140071391521d-1
      goto 13
c
48    continue
      u(25) = .323801709628694d-1
      u(26) = .970046992094627d-1
      u(27) = .161222356068892d0
      u(28) = .224763790394689d0
      u(29) = .287362487355456d0
      u(30) = .348755886292161d0
      u(31) = .408686481990717d0
      u(32) = .466902904750958d0
      u(33) = .523160974722233d0
      u(34) = .577224726083973d0
      u(35) = .628867396776514d0
      u(36) = .677872379632664d0
      u(37) = .724034130923815d0
      u(38) = .767159032515740d0
      u(39) = .807066204029443d0
      u(40) = .843588261624394d0
      u(41) = .876572020274247d0
      u(42) = .905879136715570d0
      u(43) = .931386690706554d0
      u(44) = .952987703160431d0
      u(45) = .970591592546247d0
      u(46) = .984124583722827d0
      u(47) = .993530172266351d0
      u(48) = .998771007252426d0
      a(25) = .647376968126839d-1
      a(26) = .644661644359501d-1
      a(27) = .639242385846482d-1
      a(28) = .631141922862540d-1
      a(29) = .620394231598927d-1
      a(30) = .607044391658939d-1
      a(31) = .591148396983956d-1
      a(32) = .572772921004032d-1
      a(33) = .551995036999842d-1
      a(34) = .528901894851937d-1
      a(35) = .503590355538545d-1
      a(36) = .476166584924905d-1
      a(37) = .446745608566943d-1
      a(38) = .415450829434647d-1
      a(39) = .382413510658307d-1
      a(40) = .347772225647704d-1
      a(41) = .311672278327981d-1
      a(42) = .274265097083569d-1
      a(43) = .235707608393244d-1
      a(44) = .196161604573555d-1
      a(45) = .155793157229438d-1
      a(46) = .114772345792345d-1
      a(47) = .732755390127626d-2
      a(48) = .315334605230584d-2
13    continue
c
      nq1 = nq+1
      do 15 i = 1,n
         ii   = nq1-i
         u(i) = -u(ii)
         a(i) = a(ii)
15    continue
c
      return
      end
*
******************************************************************
c  akbrdf  -  an interface between 6s and msrm
c  MSRM93  -  MultiSpectral Reflectance Model   A. Kuusk   24.03.1993
c                                               Internet:  andres@aai.ee
c
c             A. Kuusk, A multispectral canopy reflectance model, 
c                       Remote Sens. Environ.,  1994,  50(2):75-82.
c
      subroutine akbrdf(eei, thmi, uli, sli, rsl1i, wlmoy, rnci, 
     & cabi, cwi, vaii, mu, np, rm, rp, brdfint)
c  See on tegelikult juba mcrm,  aga clx ja clz on fikseeritud
c
      implicit double precision (a-h, o-z)
      integer np, mu
      integer k, j
      real eei, thmi, uli, sli, rsl1i, wlmoy, rnci, cabi, cwi, 
     & vaii, pir
      real mu1, mu2, fi
      real rm(-mu:mu), rp(np), brdfint(-mu:mu, np)
      save /count/, /soildata/, /aaa/, /ggg/, /ladak/
c
      dimension u1(10), u2(10), a1(10), a2(10)
      common /count/ jl, jj, lg, jg, lf, nnx, n1, n2, u1, u2, a1, a2
c
      double precision nnl, kk, integr
      common /leafin/ nnl, vai, kk
      common /leafout/ refl, tran
c
      double precision ke, kab, kw
      dimension refr(200), ke(200), kab(200), kw(200)
      common /dat/ refr, ke, kab, kw
c
      dimension phis1(200), phis2(200), phis3(200), phis4(200)
      common /soildata/ phis1, phis2, phis3, phis4, rsl1, rsl2, 
     & rsl3, rsl4, th2, rsl, rsoil, rr1soil, rrsoil
c
      common /aaa/ rrl, ttl, ul, sl, clmp, clmp1, bi, bd, bqint
      common /ggg/gr, gt, g, g1, th, sth, cth, th1, sth1, cth1, 
     & phi, sp, cp, th22, st, ct, st1, ct1, t10, t11, e1, e2, 
     & s2, s3, ctg, ctg1, ctt1, stt1, calph, alp2, salp2, calp2,
     & alph, salph, alpp, difmy, difsig
      common /cfresn/ rn, rk
      common /ladak/ ee, thm, sthm, cthm
      common /msrmdata/ th10, rncoef, cab, cw, bq
c
c
      data pi/3.141592653589793d0/, pir/3.14159265/
      data pi12/1.570796326794895d0/, dr/1.745329251994330d-2/
      data eps/.1d-5/, eps4/.1d-3/
c---      data lf/1/
c
*           print *, 'msrm93'
c
      integr(xx) = (1.d0 - exp(-xx))/xx
      jg = 1
*     if (lf .eq. 1) then 
c
        ee    = eei
        thm   = thmi*dr
        ul    = uli
        sl    = sli
        clz   = .9d0
        clx   = .1d0
        th2   = 45.d0*dr
        rsl1  = rsl1i
        rsl2  = -.48d0*rsl1 + .0862d0
        rsl3  = 0.d0
        rsl4  = 0.d0
        rlambda = wlmoy*1000.d0
c
        if ((rlambda .gt. 2500.d0) .or. (rlambda .lt. 404.d0)) then
           print *, 'AKBRDF: wavelength out of range'
           stop
        endif
c
        if (rlambda .le. 800.d0) then
           jl = nint((rlambda - 400.d0)/4.d0)
        else
           jl = nint((rlambda - 800.d0)/17.d0) + 100
        endif
c
        rncoef = rnci
        cab    = cabi
        cw     = cwi
        vai    = vaii
        nnl    = refr(jl)
        kk     = ke(jl) + cab*kab(jl) + cw*kw(jl)
        call leaf
c
        rn   = rncoef*nnl
        rrl  = refl - ((1.d0 - rn)/(1.d0 + rn))**2
        rrls = rrl
        ttl  = tran
c
        call soilspec
c
        cthm = cos(thm)
        sthm = sin(thm)
c
        th22 = pi12 - thm
        if (abs(th22) .lt. eps4) th22 = 0.d0
        eln = -log(1.d0 - ee)
        difmy  = abs(.059d0*eln*(thm - 1.02d0) + .02d0)
        difsig = abs(.01771d0 - .0216d0*eln*(thm - .256d0))
c
*       lf = 2
*     endif
c
      sth10 = sin(th10)
      cth10 = cos(th10)
c
      mu1 = rm(0)
      do 1 k = 1, np
      do 2 j = 1, mu
        mu2 = rm(j)
        if (j .eq. mu) then
           fi = rm(-mu)
        else
           fi = rp(k) + rm(-mu)
        endif
        th10 = acos(mu1)
        if (fi .lt. 0.) fi = fi + 2.*pir
        if (fi .gt. (2.*pir)) fi = fi - 2.*pir
        if (fi .gt. pir) fi = 2.*pir - fi
        tgt1 = tan(th10)
        xx   = tgt1*clx/sl
c
        if (xx .lt. eps) then
            clmp1 = clz
        else
            clmp1 = 1.d0 - (1.d0 - clz)*integr(xx)
        endif
c
        phi = fi
        th1 = th10
        th  = acos(mu2)
        tgt = tan(th)
        xx  = tgt*clx/sl
c
        if (xx .lt. eps) then
            clmp = clz
        else
            clmp = 1.d0 - (1.d0 - clz)*integr(xx)
        endif
c
        call msrm
        brdfint(j, k) = bq
c
  2   continue
  1   continue
c
      return
      end
*
*
******************************************************************
*

      subroutine akalbe
*    & (eei, thmi, uli, sli, rsl1i, wlmoy, rnci, cabi, cwi, vaii, albbrdf)
     & (albbrdf)
c
c   aa94.f   -  albeedo integrating msrm93 over the hemisphere 
c   A. Kuusk    23.09.1994
c
      implicit double precision (a-h, o-z)
c
*     real eei, thmi, uli, sli, rsl1i, wlmoy, rnci, cabi, cwi, vaii, albbrdf
      real albbrdf
      save /count/, /soildata/, /aaa/, /ggg/, /ladak/
c
      dimension uu(20), aa(20)
c
      dimension u1(10), u2(10), a1(10), a2(10)
      common /count/ jl, jj, lg, jg, lf, nnx, n1, n2, u1, u2, a1, a2
c
      dimension phis1(200), phis2(200), phis3(200), phis4(200)
      common /soildata/ phis1, phis2, phis3, phis4, rsl1, rsl2, 
     & rsl3, rsl4, th2, rsl, rsoil, rr1soil, rrsoil
c
      common /aaa/ rrl, ttl, ul, sl, clmp, clmp1, bi, bd, bqint
      common /ggg/gr, gt, g, g1, th, sth, cth, th1, sth1, cth1,
     & phi, sp, cp, th22, st, ct, st1, ct1, t10, t11, e1, e2,  
     & s2, s3, ctg, ctg1, ctt1, stt1, calph, alp2, salp2, calp2, 
     & alph, salph, alpp, difmy, difsig
      common /ladak/ ee, thm, sthm, cthm
c
      data pi/3.141592653589793d0/, pi1/1.5707963268d0/
c
*           print *, 'aa94'
c
      n1 = 6
      n2 = 8
c
      n  = n2 + n2
      ng = n + 1
      call dakg(uu, aa, n)
c
      do 20 i = 1, n2
         i1    = ng - i
         a2(i) = aa(i)
20       u2(i) = uu(i1)
c
      n  = n1 + n1
      ng = n + 1
      call dakg(uu, aa, n)
c
      do 21 i = 1, n1
         i1    = ng - i
         a1(i) = aa(i)
21       u1(i) = uu(i1)
c
      bdd = 0.d0
      do 10 i2 = 1, n2
         th  = (1.d0 - u2(i2))*pi1
         sth = sin(th)
         cth = cos(th)
c
         call akd
c
         bdd = bdd + a2(i2)*bqint*sth*cth
10    continue
c
      albbrdf = bdd*pi
c
      return
      end
*
******************************************************************
*
      subroutine atmref (iaer,tamoy,trmoy,pizmoy,tamoyp,trmoyp,palt,
     s               phi,xmus,xmuv,
     s               phirad,nt,mu,np,rm,gb,rp,
     a                   rorayl,roaero,romix,xlm1,xlm2)
      integer mu,np
      real rm(-mu:mu),rp(np),gb(-mu:mu),xlm1(-mu:mu,np)
      real xlm2(-mu:mu,np)
      real tamoy,trmoy,pizmoy,tamoyp,trmoyp,palt
      real phi,xmus,xmuv,phirad
      real rorayl,roaero,romix,delta,sigma,tamol,tamolp
      integer iaer,nt
      common /sixs_del/ delta,sigma
c     atmospheric reflectances
c
      rorayl=0.
      roaero=0.
c     rayleigh reflectance 3 cases (satellite,plane,ground)
      if(palt.lt.900..and.palt.gt.0.0)then
        rm(-mu)=-xmuv
        rm(mu)=xmuv
        rm(0)=-xmus
        tamol=0.
        tamolp=0.
      call os(tamol,trmoy,pizmoy,tamolp,trmoyp,palt,
     s               phirad,nt,mu,np,rm,gb,rp,
     s                     xlm1)
        rorayl=xlm1(-mu,1)/xmus
        else
        if (palt.le.0.0) then
           rorayl=0.
           else
           call chand(phi,xmuv,xmus,trmoy,rorayl)
           endif
        endif
c
      if (iaer.eq.0) then
         romix=rorayl
         return
         endif
c     rayleigh+aerosol=romix,aerosol=roaero reflectance computed
c     using sucessive order of scattering method
c     3 cases: satellite,plane,ground
      if(palt.gt.0.0) then
        rm(-mu)=-xmuv
        rm(mu)=xmuv
        rm(0)=-xmus
c      write(6,*) "input OS",tamoy,trmoy,pizmoy,tamoyp,trmoyp,palt,
c    s phirad
        call os(tamoy,trmoy,pizmoy,tamoyp,trmoyp,palt,
     s               phirad,nt,mu,np,rm,gb,rp,
     s                     xlm2)
        romix=(xlm2(-mu,1)/xmus)
        tamol=0.
        tamolp=0.
        call os(tamoy,tamol,pizmoy,tamoyp,tamolp,palt,
     s               phirad,nt,mu,np,rm,gb,rp,
     s                     xlm2)
        roaero=(xlm2(-mu,1)/xmus)
        else
        roaero=0.
        romix=0.
        endif
c     write(6,*) " OS: ",rorayl,roaero,romix
      return
      end
      subroutine  bbm
      common/sixs_aerbas/ph(10,83)
      real phr(10,83)
      real ph
      integer i,j
c
c  model: biomass burning
      data ((phr(i,j),j=1,83),i=1,1)/
     &0.2150,0.2122,0.2027,0.1928,0.1884,0.1905,0.1952,0.1983,
     &0.1980,0.1954,0.1918,0.1874,0.1819,0.1752,0.1680,0.1612,
     &0.1553,0.1501,0.1457,0.1417,0.1382,0.1351,0.1326,0.1308,
     &0.1296,0.1292,0.1293,0.1299,0.1310,0.1328,0.1353,0.1387,
     &0.1429,0.1480,0.1539,0.1606,0.1682,0.1770,0.1870,0.1984,
     &0.2115,0.2186,0.2263,0.2432,0.2622,0.2838,0.3082,0.3358,
     &0.3671,0.4024,0.4423,0.4875,0.5386,0.5968,0.6630,0.7387,
     &0.8253,0.9247,1.0387,1.1695,1.3192,1.4909,1.6883,1.9162,
     &2.1797,2.4841,2.8350,3.2382,3.7008,4.2315,4.8393,5.5328,
     &6.3184,7.2028,8.1966,9.3190,10.591,12.016,13.541,15.036,
     &16.295,17.092,17.290/
      data ((phr(i,j),j=1,83),i=2,2)/
     &0.2180,0.2160,0.2091,0.2007,0.1951,0.1943,0.1972,0.2005,
     &0.2013,0.1986,0.1934,0.1874,0.1819,0.1771,0.1724,0.1673,
     &0.1619,0.1565,0.1518,0.1480,0.1449,0.1426,0.1408,0.1395,
     &0.1387,0.1383,0.1385,0.1392,0.1406,0.1427,0.1456,0.1492,
     &0.1535,0.1585,0.1644,0.1713,0.1793,0.1887,0.1995,0.2119,
     &0.2261,0.2339,0.2421,0.2601,0.2803,0.3029,0.3284,0.3571,
     &0.3896,0.4266,0.4687,0.5166,0.5710,0.6328,0.7029,0.7826,
     &0.8733,0.9769,1.0955,1.2314,1.3869,1.5649,1.7685,2.0010,
     &2.2665,2.5691,2.9134,3.3049,3.7496,4.2547,4.8268,5.4728,
     &6.1989,7.0122,7.9194,8.9236,10.016,11.166,12.309,13.350,
     &14.175,14.677,14.799/
      data ((phr(i,j),j=1,83),i=3,3)/
     &0.2171,0.2154,0.2091,0.2012,0.1955,0.1939,0.1960,0.1992,
     &0.2006,0.1987,0.1940,0.1879,0.1820,0.1770,0.1727,0.1684,
     &0.1638,0.1590,0.1544,0.1504,0.1473,0.1450,0.1433,0.1422,
     &0.1416,0.1414,0.1418,0.1426,0.1439,0.1459,0.1486,0.1522,
     &0.1566,0.1619,0.1681,0.1752,0.1835,0.1930,0.2039,0.2163,
     &0.2305,0.2383,0.2466,0.2650,0.2857,0.3090,0.3352,0.3647,
     &0.3981,0.4358,0.4785,0.5269,0.5816,0.6437,0.7143,0.7949,
     &0.8870,0.9923,1.1126,1.2501,1.4072,1.5866,1.7913,2.0244,
     &2.2894,2.5902,2.9315,3.3191,3.7594,4.2591,4.8242,5.4598,
     &6.1705,6.9612,7.8359,8.7939,9.8227,10.889,11.934,12.873,
     &13.609,14.053,14.161/
      data ((phr(i,j),j=1,83),i=4,4)/	
     &0.2183,0.2168,0.2113,0.2040,0.1981,0.1956,0.1966,0.1992,
     &0.2011,0.2003,0.1965,0.1907,0.1843,0.1786,0.1740,0.1701,
     &0.1664,0.1624,0.1583,0.1543,0.1510,0.1484,0.1466,0.1454,
     &0.1448,0.1447,0.1451,0.1461,0.1476,0.1497,0.1525,0.1560,
     &0.1605,0.1660,0.1725,0.1800,0.1886,0.1984,0.2095,0.2221,
     &0.2364,0.2442,0.2526,0.2710,0.2920,0.3158,0.3429,0.3735,
     &0.4081,0.4469,0.4906,0.5399,0.5957,0.6591,0.7311,0.8132,
     &0.9068,1.0134,1.1350,1.2737,1.4317,1.6116,1.8158,2.0475,
     &2.3101,2.6080,2.9466,3.3317,3.7694,4.2645,4.8209,5.4417,
     &6.1298,6.8888,7.7208,8.6227,9.5794,10.558,11.502,12.339,
     &12.989,13.377,13.471/     
      data ((phr(i,j),j=1,83),i=5,5)/
     &0.2249,0.2239,0.2197,0.2137,0.2078,0.2036,0.2019,0.2022,
     &0.2034,0.2038,0.2022,0.1985,0.1929,0.1865,0.1803,0.1751,
     &0.1711,0.1679,0.1651,0.1624,0.1597,0.1571,0.1549,0.1533,
     &0.1525,0.1525,0.1532,0.1545,0.1564,0.1589,0.1622,0.1662,
     &0.1710,0.1767,0.1832,0.1907,0.1995,0.2097,0.2215,0.2351,
     &0.2505,0.2589,0.2679,0.2874,0.3094,0.3341,0.3621,0.3936,
     &0.4294,0.4698,0.5157,0.5676,0.6266,0.6934,0.7691,0.8549,
     &0.9521,1.0619,1.1858,1.3254,1.4828,1.6604,1.8614,2.0895,
     &2.3492,2.6448,2.9802,3.3587,3.7823,4.2529,4.7729,5.3459,
     &5.9763,6.6683,7.4225,8.2312,9.0749,9.9191,10.715,11.405,
     &11.931,12.241,12.316/
      data ((phr(i,j),j=1,83),i=6,6)/
     &0.2268,0.2259,0.2225,0.2173,0.2117,0.2070,0.2041,0.2029,
     &0.2031,0.2034,0.2029,0.2008,0.1970,0.1919,0.1861,0.1806,
     &0.1758,0.1721,0.1693,0.1673,0.1655,0.1638,0.1622,0.1607,
     &0.1597,0.1593,0.1597,0.1609,0.1628,0.1654,0.1687,0.1728,
     &0.1778,0.1838,0.1908,0.1990,0.2083,0.2189,0.2309,0.2446,
     &0.2602,0.2687,0.2779,0.2981,0.3210,0.3470,0.3763,0.4093,
     &0.4464,0.4883,0.5354,0.5887,0.6488,0.7166,0.7931,0.8793,
     &0.9763,1.0854,1.2083,1.3469,1.5037,1.6817,1.8840,2.1140,
     &2.3747,2.6685,2.9974,3.3632,3.7679,4.2147,4.7080,5.2537,
     &5.8571,6.5215,7.2440,8.0132,8.8062,9.5882,10.315,10.936,
     &11.403,11.678,11.743/
      data ((phr(i,j),j=1,83),i=7,7)/
     &0.2427,0.2421,0.2399,0.2362,0.2317,0.2269,0.2224,0.2187,
     &0.2159,0.2139,0.2124,0.2110,0.2094,0.2072,0.2041,0.2004,
     &0.1962,0.1917,0.1875,0.1839,0.1810,0.1790,0.1779,0.1775,
     &0.1776,0.1782,0.1792,0.1805,0.1822,0.1846,0.1877,0.1917,
     &0.1968,0.2031,0.2106,0.2194,0.2295,0.2412,0.2545,0.2697,
     &0.2869,0.2963,0.3063,0.3284,0.3532,0.3812,0.4125,0.4476,
     &0.4868,0.5305,0.5793,0.6339,0.6951,0.7639,0.8414,0.9288,
     &1.0277,1.1395,1.2655,1.4074,1.5662,1.7433,1.9398,2.1574,
     &2.3979,2.6641,2.9597,3.2895,3.6589,4.0736,4.5385,5.0562,
     &5.6261,6.2424,6.8937,7.5622,8.2243,8.8516,9.4132,9.8788,
     &10.221,10.418,10.465/
      data ((phr(i,j),j=1,83),i=8,8)/
     &0.3408,0.3406,0.3396,0.3380,0.3356,0.3327,0.3292,0.3253,
     &0.3210,0.3165,0.3119,0.3072,0.3026,0.2981,0.2939,0.2898,
     &0.2861,0.2827,0.2797,0.2770,0.2747,0.2728,0.2712,0.2701,
     &0.2693,0.2690,0.2693,0.2700,0.2715,0.2737,0.2768,0.2808,
     &0.2861,0.2926,0.3005,0.3101,0.3214,0.3346,0.3499,0.3675,
     &0.3875,0.3984,0.4100,0.4354,0.4636,0.4951,0.5300,0.5686,
     &0.6114,0.6588,0.7114,0.7697,0.8346,0.9068,0.9874,1.0773,
     &1.1778,1.2898,1.4147,1.5535,1.7072,1.8766,2.0625,2.2649,
     &2.4840,2.7191,2.9691,3.2325,3.5070,3.7899,4.0777,4.3667,
     &4.6524,4.9302,5.1951,5.4422,5.6667,5.8638,6.0293,6.1597,
     &6.2518,6.3038,6.3160/
      data ((phr(i,j),j=1,83),i=9,9)/
     &0.4735,0.4733,0.4725,0.4711,0.4690,0.4664,0.4632,0.4596,
     &0.4554,0.4507,0.4457,0.4404,0.4347,0.4289,0.4229,0.4168,
     &0.4106,0.4046,0.3987,0.3930,0.3876,0.3825,0.3779,0.3738,
     &0.3704,0.3676,0.3656,0.3646,0.3645,0.3655,0.3677,0.3712,
     &0.3762,0.3827,0.3910,0.4011,0.4134,0.4278,0.4447,0.4643,
     &0.4868,0.4992,0.5124,0.5414,0.5742,0.6109,0.6519,0.6974,
     &0.7479,0.8036,0.8648,0.9317,1.0047,1.0838,1.1694,1.2615,
     &1.3601,1.4653,1.5768,1.6946,1.8182,1.9474,2.0814,2.2197,
     &2.3614,2.5057,2.6516,2.7980,2.9435,3.0870,3.2271,3.3623,
     &3.4914,3.6129,3.7254,3.8277,3.9184,3.9966,4.0611,4.1113,
     &4.1465,4.1662,4.1709/
      data ((phr(i,j),j=1,83),i=10,10)/
     &0.7907,0.7905,0.7895,0.7878,0.7852,0.7820,0.7780,0.7733,
     &0.7679,0.7619,0.7553,0.7481,0.7405,0.7324,0.7239,0.7151,
     &0.7061,0.6968,0.6875,0.6782,0.6690,0.6599,0.6512,0.6428,
     &0.6349,0.6276,0.6211,0.6154,0.6107,0.6071,0.6047,0.6037,
     &0.6042,0.6063,0.6102,0.6160,0.6239,0.6339,0.6462,0.6609,
     &0.6782,0.6878,0.6981,0.7207,0.7461,0.7743,0.8055,0.8396,
     &0.8768,0.9168,0.9599,1.0058,1.0545,1.1060,1.1601,1.2166,
     &1.2753,1.3362,1.3988,1.4630,1.5284,1.5948,1.6618,1.7290,
     &1.7962,1.8627,1.9284,1.9927,2.0553,2.1156,2.1734,2.2281,
     &2.2795,2.3271,2.3705,2.4095,2.4438,2.4729,2.4969,2.5154,
     &2.5283,2.5355,2.5372/
      do 1 i=1,10
      do 1 j=1,83
      ph(i,j)=phr(i,j)
   1  continue
      return
      end
      subroutine  bdm
      common/sixs_aerbas/ph(10,83)
      real phr(10,83)
      real ph
      integer i,j
c
c  model:background desert
      data ((phr(i,j),j=1,83),i=1,1)/
     &0.8352,0.8057,0.7377,0.6569,0.5760,0.5032,0.4427,0.3969,
     &0.3646,0.3385,0.3125,0.2863,0.2611,0.2380,0.2175,0.1998,
     &0.1848,0.1722,0.1619,0.1536,0.1469,0.1416,0.1376,0.1347,
     &0.1329,0.1319,0.1319,0.1327,0.1343,0.1366,0.1397,0.1437,
     &0.1485,0.1541,0.1607,0.1682,0.1768,0.1865,0.1974,0.2097,
     &0.2235,0.2309,0.2388,0.2559,0.2749,0.2960,0.3196,0.3459,
     &0.3750,0.4073,0.4432,0.4831,0.5276,0.5774,0.6331,0.6954,
     &0.7652,0.8433,0.9308,1.0291,1.1399,1.2648,1.4060,1.5661,
     &1.7479,1.9552,2.1925,2.4657,2.7822,3.1529,3.5921,4.1201,
     &4.7671,5.5787,6.6249,8.0218,9.9742,12.864,17.461,25.540,
     &42.106,87.294,183.39/
      data ((phr(i,j),j=1,83),i=2,2)/
     &0.8002,0.7733,0.7063,0.6273,0.5489,0.4793,0.4227,0.3810,
     &0.3524,0.3297,0.3071,0.2839,0.2611,0.2399,0.2207,0.2039,
     &0.1895,0.1773,0.1672,0.1590,0.1523,0.1471,0.1431,0.1401,
     &0.1382,0.1373,0.1372,0.1380,0.1396,0.1419,0.1451,0.1490,
     &0.1539,0.1596,0.1663,0.1739,0.1826,0.1925,0.2036,0.2161,
     &0.2301,0.2377,0.2458,0.2632,0.2826,0.3042,0.3284,0.3553,
     &0.3852,0.4184,0.4553,0.4964,0.5424,0.5938,0.6514,0.7159,
     &0.7882,0.8693,0.9603,1.0626,1.1779,1.3081,1.4554,1.6225,
     &1.8124,2.0287,2.2763,2.5611,2.8904,3.2747,3.7281,4.2698,
     &4.9281,5.7448,6.7836,8.1481,10.017,12.714,16.878,23.920,
     &37.639,72.434,130.26/
      data ((phr(i,j),j=1,83),i=3,3)/
     &0.7899,0.7637,0.6974,0.6190,0.5414,0.4728,0.4173,0.3766,
     &0.3489,0.3271,0.3054,0.2830,0.2609,0.2402,0.2214,0.2048,
     &0.1906,0.1786,0.1685,0.1603,0.1537,0.1485,0.1445,0.1415,
     &0.1396,0.1387,0.1386,0.1394,0.1409,0.1433,0.1464,0.1504,
     &0.1553,0.1610,0.1677,0.1754,0.1841,0.1940,0.2052,0.2178,
     &0.2319,0.2395,0.2476,0.2651,0.2846,0.3064,0.3307,0.3578,
     &0.3879,0.4214,0.4585,0.5000,0.5464,0.5982,0.6563,0.7215,
     &0.7944,0.8763,0.9682,1.0715,1.1881,1.3197,1.4685,1.6375,
     &1.8294,2.0481,2.2983,2.5859,2.9183,3.3060,3.7626,4.3073,
     &4.9677,5.7846,6.8200,8.1736,10.017,12.661,16.709,23.484,
     &36.489,68.980,119.09/
      data ((phr(i,j),j=1,83),i=4,4)/	
     &0.7770,0.7516,0.6862,0.6087,0.5323,0.4648,0.4106,0.3713,
     &0.3447,0.3239,0.3032,0.2817,0.2604,0.2403,0.2220,0.2058,
     &0.1918,0.1800,0.1700,0.1619,0.1553,0.1501,0.1461,0.1432,
     &0.1413,0.1403,0.1402,0.1410,0.1426,0.1449,0.1481,0.1521,
     &0.1570,0.1628,0.1695,0.1772,0.1860,0.1959,0.2072,0.2199,
     &0.2340,0.2417,0.2498,0.2675,0.2871,0.3091,0.3336,0.3608,
     &0.3912,0.4250,0.4625,0.5044,0.5512,0.6036,0.6623,0.7282,
     &0.8019,0.8848,0.9777,1.0823,1.2003,1.3336,1.4844,1.6555,
     &1.8498,2.0712,2.3245,2.6155,2.9515,3.3429,3.8033,4.3511,
     &5.0134,5.8298,6.8600,8.1995,10.012,12.589,16.495,22.954,
     &35.131,65.032,107.60/
      data ((phr(i,j),j=1,83),i=5,5)/
     &0.7483,0.7247,0.6618,0.5867,0.5127,0.4480,0.3967,0.3601,
     &0.3357,0.3169,0.2982,0.2787,0.2590,0.2403,0.2230,0.2076,
     &0.1942,0.1827,0.1730,0.1651,0.1586,0.1535,0.1495,0.1466,
     &0.1447,0.1437,0.1437,0.1445,0.1460,0.1484,0.1516,0.1557,
     &0.1606,0.1664,0.1732,0.1810,0.1899,0.2000,0.2114,0.2242,
     &0.2386,0.2464,0.2546,0.2725,0.2925,0.3147,0.3396,0.3674,
     &0.3983,0.4327,0.4710,0.5138,0.5616,0.6151,0.6752,0.7425,
     &0.8180,0.9029,0.9982,1.1054,1.2264,1.3631,1.5178,1.6933,
     &1.8926,2.1196,2.3789,2.6765,3.0196,3.4181,3.8851,4.4383,
     &5.1031,5.9162,6.9324,8.2384,9.9800,12.414,16.024,21.8313,
     &32.417,57.406,86.131/
      data ((phr(i,j),j=1,83),i=6,6)/
     &0.7290,0.7065,0.6456,0.5721,0.5001,0.4373,0.3879,0.3528,
     &0.3298,0.3122,0.2948,0.2764,0.2579,0.2400,0.2234,0.2085,
     &0.1955,0.1843,0.1748,0.1670,0.1606,0.1555,0.1516,0.1488,
     &0.1469,0.1459,0.1459,0.1467,0.1483,0.1507,0.1539,0.1579,
     &0.1629,0.1688,0.1756,0.1835,0.1924,0.2026,0.2141,0.2270,
     &0.2415,0.2494,0.2577,0.2758,0.2959,0.3185,0.3436,0.3717,
     &0.4030,0.4378,0.4766,0.5199,0.5684,0.6227,0.6836,0.7519,
     &0.8285,0.9146,1.0114,1.1203,1.2432,1.3821,1.5392,1.7174,
     &1.9197,2.1501,2.4131,2.7146,3.0617,3.4642,3.9347,4.4902,
     &5.1552,5.9644,6.9694,8.2514,9.9446,12.284,15.706,21.119,
     &30.763,52.922,75.133/
      data ((phr(i,j),j=1,83),i=7,7)/
     &0.6834,0.6633,0.6079,0.5390,0.4716,0.4134,0.3682,0.3368,
     &0.3165,0.3014,0.2865,0.2708,0.2546,0.2387,0.2237,0.2101,
     &0.1980,0.1875,0.1786,0.1711,0.1650,0.1601,0.1564,0.1536,
     &0.1518,0.1509,0.1509,0.1517,0.1534,0.1558,0.1591,0.1632,
     &0.1683,0.1742,0.1812,0.1892,0.1983,0.2087,0.2204,0.2337,
     &0.2485,0.2565,0.2650,0.2835,0.3042,0.3273,0.3531,0.3819,
     &0.4140,0.4499,0.4898,0.5345,0.5844,0.6405,0.7033,0.7739,
     &0.8531,0.9420,1.0421,1.1547,1.2818,1.4253,1.5877,1.7716,
     &1.9804,2.2176,2.4880,2.7971,3.1520,3.5615,4.0374,4.5953,
     &5.2568,6.0522,7.0261,8.2465,9.8238,11.947,14.951,19.512,
     &27.207,43.691,55.647/
      data ((phr(i,j),j=1,83),i=8,8)/
     &0.5664,0.5524,0.5105,0.4593,0.4056,0.3604,0.3252,0.3017,
     &0.2868,0.2764,0.2666,0.2562,0.2452,0.2340,0.2231,0.2127,
     &0.2033,0.1949,0.1876,0.1814,0.1763,0.1721,0.1689,0.1665,
     &0.1651,0.1644,0.1647,0.1657,0.1675,0.1702,0.1737,0.1781,
     &0.1835,0.1898,0.1972,0.2056,0.2153,0.2264,0.2388,0.2529,
     &0.2687,0.2773,0.2864,0.3062,0.3284,0.3531,0.3808,0.4118,
     &0.4464,0.4850,0.5281,0.5763,0.6303,0.6908,0.7586,0.8347,
     &0.9201,1.0159,1.1236,1.2446,1.3810,1.5344,1.7076,1.9029,
     &2.1236,2.3730,2.6551,2.9753,3.3383,3.7523,4.2252,4.7683,
     &5.3971,6.1289,6.9937,8.0270,9.2901,10.873,12.948,15.760,
     &20.227,26.155,28.327/
      data ((phr(i,j),j=1,83),i=9,9)/
     &0.5017,0.4916,0.4574,0.4166,0.3755,0.3366,0.3067,0.2874,
     &0.2748,0.2660,0.2585,0.2504,0.2418,0.2329,0.2241,0.2156,
     &0.2078,0.2007,0.1945,0.1891,0.1846,0.1810,0.1781,0.1761,
     &0.1750,0.1746,0.1750,0.1762,0.1782,0.1810,0.1848,0.1894,
     &0.1950,0.2016,0.2093,0.2181,0.2283,0.2398,0.2528,0.2676,
     &0.2841,0.2931,0.3026,0.3234,0.3466,0.3726,0.4016,0.4341,
     &0.4704,0.5108,0.5560,0.6065,0.6630,0.7261,0.7970,0.8761,
     &0.9649,1.0644,1.1759,1.3009,1.4415,1.5992,1.7762,1.9755,
     &2.1997,2.4512,2.7345,3.0541,3.4136,3.8188,4.2785,4.7998,
     &5.3909,6.0685,6.8504,7.7572,8.8313,10.118,11.724,13.933,
     &16.806,19.370,20.119/
      data ((phr(i,j),j=1,83),i=10,10)/
     &0.4481,0.4411,0.4148,0.3788,0.3444,0.3172,0.2972,0.2822,
     &0.2711,0.2632,0.2572,0.2514,0.2450,0.2379,0.2310,0.2245,
     &0.2183,0.2126,0.2074,0.2030,0.1993,0.1963,0.1939,0.1923,
     &0.1915,0.1914,0.1920,0.1934,0.1957,0.1988,0.2027,0.2076,
     &0.2135,0.2206,0.2287,0.2381,0.2488,0.2611,0.2750,0.2906,
     &0.3082,0.3178,0.3279,0.3501,0.3748,0.4024,0.4332,0.4677,
     &0.5062,0.5491,0.5968,0.6500,0.7094,0.7758,0.8499,0.9325,
     &1.0245,1.1273,1.2424,1.3710,1.5144,1.6743,1.8527,2.0524,
     &2.2759,2.5253,2.8026,3.1112,3.4553,3.8394,4.2681,4.7465,
     &5.2801,5.8742,6.5358,7.2843,8.1602,9.2141,10.458,11.804,
     &13.032,13.853,14.061/
      do 1 i=1,10
      do 1 j=1,83
      ph(i,j)=phr(i,j)
   1  continue
      return
      end
	subroutine chand (xphi,xmuv,xmus,xtau
     s			,xrray)
c input parameters: xphi,xmus,xmuv,xtau
c xphi: azimuthal difference between sun and observation (xphi=0,
c in backscattering) and expressed in degree (0.:360.)
c xmus: cosine of the sun zenith angle
c xmuv: cosine of the observation zenith angle
c xtau: molecular optical depth
c output parameter: xrray : molecular reflectance (0.:1.)
c constant : xdep: depolarization factor (0.0279)
	real xdep,pl(10)
	real fs0,fs1,fs2
	real as0(10),as1(2),as2(2)
        real xphi,xmus,fac,xmuv,xtau,xrray,pi,phios,xcosf1,xcosf2
        real xcosf3,xbeta2,xfd,xph1,xph2,xph3,xitm, xp1, xp2, xp3
        real cfonc1,cfonc2,cfonc3,xlntau,xitot1,xitot2,xitot3
        integer i
	data (as0(i),i=1,10) /.33243832,-6.777104e-02,.16285370
     s	,1.577425e-03,-.30924818,-1.240906e-02,-.10324388
     s	,3.241678e-02,.11493334,-3.503695e-02/
	data (as1(i),i=1,2) /.19666292, -5.439061e-02/
	data (as2(i),i=1,2) /.14545937,-2.910845e-02/
	pi=3.1415927
	fac=pi/180.
	phios=180.-xphi
	xcosf1=1.
	xcosf2=cos(phios*fac)
	xcosf3=cos(2*phios*fac)
	xbeta2=0.5
	xdep=0.0279
	xfd=xdep/(2-xdep)
	xfd=(1-xfd)/(1+2*xfd)
	xph1=1+(3*xmus*xmus-1)*(3*xmuv*xmuv-1)*xfd/8.
	xph2=-xmus*xmuv*sqrt(1-xmus*xmus)*sqrt(1-xmuv*xmuv)
	xph2=xph2*xfd*xbeta2*1.5
	xph3=(1-xmus*xmus)*(1-xmuv*xmuv)
	xph3=xph3*xfd*xbeta2*0.375
	xitm=(1-exp(-xtau*(1/xmus+1/xmuv)))*xmus/(4*(xmus+xmuv))
	xp1=xph1*xitm
	xp2=xph2*xitm
	xp3=xph3*xitm
	xitm=(1-exp(-xtau/xmus))*(1-exp(-xtau/xmuv))
	cfonc1=xph1*xitm
	cfonc2=xph2*xitm
	cfonc3=xph3*xitm
	xlntau=log(xtau)
	pl(1)=1.
	pl(2)=xlntau
	pl(3)=xmus+xmuv
	pl(4)=xlntau*pl(3)
	pl(5)=xmus*xmuv
	pl(6)=xlntau*pl(5)
	pl(7)=xmus*xmus+xmuv*xmuv
	pl(8)=xlntau*pl(7)
	pl(9)=xmus*xmus*xmuv*xmuv
	pl(10)=xlntau*pl(9)
	fs0=0.
	do i=1,10
	fs0=fs0+pl(i)*as0(i)
	enddo
	fs1=pl(1)*as1(1)+pl(2)*as1(2)
	fs2=pl(1)*as2(1)+pl(2)*as2(2)
	xitot1=xp1+cfonc1*fs0*xmus
	xitot2=xp2+cfonc2*fs1*xmus
	xitot3=xp3+cfonc3*fs2*xmus
	xrray=xitot1*xcosf1
	xrray=xrray+xitot2*xcosf2*2
	xrray=xrray+xitot3*xcosf3*2
	xrray=xrray/xmus
	return
	end
      subroutine csalbr(xtau,xalb)
      real xtau,xalb,fintexp3
      xalb=(3*xtau-fintexp3(xtau)*(4+2*xtau)+2*exp(-xtau))
      xalb=xalb/(4.+3*xtau)
      return
      end
      real function fintexp3(xtau)
      real xx,xtau,fintexp1
      xx=(exp(-xtau)*(1.-xtau)+xtau*xtau*fintexp1(xtau))/2.
      fintexp3=xx
      return
      end
      real function fintexp1(xtau)
c accuracy 2e-07... for 0<xtau<1
      real xx,a(0:5),xtau,xftau
      integer i
      data (a(i),i=0,5) /-.57721566,0.99999193,-0.24991055,
     c                  0.05519968,-0.00976004,0.00107857/
      xx=a(0)
      xftau=1.
      do i=1,5
      xftau=xftau*xtau
      xx=xx+a(i)*xftau
      enddo
      fintexp1=xx-log(xtau)
      return
      end
      subroutine discom (idatmp,iaer,xmus,xmuv,phi,
     a                   taer55,taer55p,palt,
     a                 phirad,nt,mu,np,rm,gb,rp,
     a                   ftray,xlm1,xlm2)
      integer mu,np
      real rm(-mu:mu),rp(np),gb(-mu:mu)
      real ftray,xlm1(-mu:mu,np),xlm2(-mu:mu,np)
      real xmus,xmuv,phi
      real taer55,taer55p,palt,phirad,ext,ome,gasym,phase,roatm
      real dtdir,dtdif,utdir,utdif,sphal,wldis,trayl,traypl,s
      real wlinf,wlsup,phasel,pdgs,cgaus,pha,betal,wl,tray,trayp,taer
      real taerp,piza,tamoy,tamoyp,pizmoy,rorayl
      real roaero,romix,ddirtt,ddiftt,udirtt,udiftt,sphalbt,ddirtr
      real ddiftr,udirtr,udiftr,sphalbr,ddirta,ddifta,udirta,udifta
      real sphalba,coeff
      integer idatmp,iaer,nt,l,k
      common /sixs_aer/ext(10),ome(10),gasym(10),phase(10)
      common /sixs_disc/ roatm(3,10),dtdir(3,10),dtdif(3,10),
     a utdir(3,10),utdif(3,10),sphal(3,10),wldis(10),trayl(10),
     a traypl(10)
      common /sixs_ffu/s(1501),wlinf,wlsup
      common /sixs_sos/phasel(10,83),cgaus(83),pdgs(83)
      common /sixs_trunc/pha(83),betal(0:80)

c     computation of all scattering parameters at wavelength
c     discrete values,so we
c     can interpolate at any wavelength
 
      do 50 l=1,10
      wl=wldis(l)
      if ((wlsup.lt.wldis(1)).and.(l.le.2)) goto 30
      if (wlinf.gt.wldis(10).and.(l.ge.9)) goto 30
      if ((l.lt.10).and.(wldis(l).lt.wlinf).and.
     a     (wldis(l+1).lt.wlinf))
     a     goto 50
      if ((l.gt.1).and.(wldis(l).gt.wlsup).and.
     a      (wldis(l-1).gt.wlsup))
     a     goto 50
 
c     computation of rayleigh optical depth at wl
 
 30   call odrayl(wl,
     a           tray)
 
c plane case discussed here above
 
      if (idatmp.eq.0.or.idatmp.eq.4) then
	  if (idatmp.eq.4) trayp=tray
	  if (idatmp.eq.0) trayp=0.
	  else
          trayp=tray*ftray
      endif
      trayl(l)=tray
      traypl(l)=trayp
 
c     computation of aerosol optical properties at wl
 
      taer=taer55*ext(l)/ext(4)
      taerp=taer55p*ext(l)/ext(4)
      piza=ome(l)
c
c     computation of atmospheric reflectances
c               rorayl is rayleigh ref
c               roaero is aerosol ref
c     call plegen to decompose aerosol phase function in Betal
      if (iaer.ne.0) then
      do k=1,83
      pha(k)=phasel(l,k)
      enddo
      call trunca(coeff)
      endif
c     write(6,*) 'truncation coefficient ',coeff
      tamoy=taer*(1.-piza*coeff)
      tamoyp=taerp*(1.-piza*coeff)
      pizmoy=piza*(1.-coeff)/(1.-piza*coeff)
c     write(6,*) 'tray,trayp,tamoy,tamoyp,pizmoy,piza,taer,taerp',
c    S            tray,trayp,tamoy,tamoyp,pizmoy,piza,taer,taerp
c     do i=0,80
c     write(6,'(A5,I2.2,1X,E13.7)') 'betal',i,betal(i)
c     enddo
c
      call atmref(iaer,tamoy,tray,pizmoy,tamoyp,trayp,palt,
     a               phi,xmus,xmuv,
     s               phirad,nt,mu,np,rm,gb,rp,
     a                   rorayl,roaero,romix,xlm1,xlm2)
c     write(6,*) 'wl,refrayl,refaero,refmix',wl,rorayl,roaero,romix
c     computation of scattering transmitances (direct and diffuse)
c     first time for rayleigh ,next total (rayleigh+aerosols)
      call scatra (tamoy,tamoyp,tray,trayp,pizmoy,
     a      palt,nt,mu,rm,gb,xmus,xmuv,
     a             ddirtt,ddiftt,udirtt,udiftt,sphalbt,
     a             ddirtr,ddiftr,udirtr,udiftr,sphalbr,
     a             ddirta,ddifta,udirta,udifta,sphalba)
      roatm(1,l)=rorayl
      roatm(2,l)=romix
      roatm(3,l)=roaero
      dtdir(1,l)=ddirtr
      dtdif(1,l)=ddiftr
      dtdir(2,l)=ddirtt
      dtdif(2,l)=ddiftt
      dtdir(3,l)=ddirta
      dtdif(3,l)=ddifta
      utdir(1,l)=udirtr
      utdif(1,l)=udiftr
      utdir(2,l)=udirtt
      utdif(2,l)=udiftt
      utdir(3,l)=udirta
      utdif(3,l)=udifta
      sphal(1,l)=sphalbr
      sphal(2,l)=sphalbt
      sphal(3,l)=sphalba
   50 continue
      return
      end
      subroutine discre(ta,ha,tr,hr,it,nt,yy,dd,ppp2,ppp1,
     s     zx)
      real ta,ha,tr,hr,yy,dd,ppp2,ppp1,zx,dt,ti,y1,y2,y3,x2
      real xd,delta,ecart
      double precision xx
      integer it,nt
      if (ha.ge.7.) then
          call print_error
     s    ('check aerosol measurements or plane altitude')
          return
          endif
      if (it.eq.0) then
         dt=1.e-17
         else
         dt=2.*(ta+tr-yy)/(nt-it+1.)
      endif
  99  dt=dt/2.
      ti=yy+dt
      y1=ppp2
      y3=ppp1
  706 y2=(y1+y3)*0.5
      xx=-y2/ha
      if (xx.lt.-18) then
         x2=tr*exp(-y2/hr)
         else
         x2=ta*dexp(xx)+tr*exp(-y2/hr)
         endif
      xd=abs(ti-x2)
      if(xd.lt.0.00001) go to 705
      if(ti-x2) 701,703,703
  701 y3=y2
      go to 706
  703 y1=y2
      go to 706
  705 zx=y2
      delta=1./(1.+ta*hr/tr/ha*exp((zx-ppp1)*(1./hr-1./ha)))
      ecart=0
      if(dd.ne.0) ecart=abs((dd-delta)/dd)
      if((ecart.gt.0.75).and.(it.ne.0)) go to 99
      return
      end
      subroutine   dust
      common /sixs_aerbas/ ph(10,83)
      real phr(10,83),ph
      integer i,j
c
c    model: dust-like
c
            DATA ((PHR(I,J),J=1,83),I=01,01) /
     *0.2021E+00,0.2079E+00,0.2462E+00,0.2310E+00,0.2069E+00,0.1883E+00,
     *0.1750E+00,0.1624E+00,0.1458E+00,0.1241E+00,0.1013E+00,0.8379E-01,
     *0.7097E-01,0.6207E-01,0.5595E-01,0.5174E-01,0.4879E-01,0.4675E-01,
     *0.4531E-01,0.4435E-01,0.4373E-01,0.4337E-01,0.4324E-01,0.4330E-01,
     *0.4353E-01,0.4392E-01,0.4449E-01,0.4522E-01,0.4612E-01,0.4721E-01,
     *0.4850E-01,0.5001E-01,0.5177E-01,0.5381E-01,0.5616E-01,0.5885E-01,
     *0.6191E-01,0.6540E-01,0.6936E-01,0.7383E-01,0.7889E-01,0.8168E-01,
     *0.8459E-01,0.9096E-01,0.9808E-01,0.1060E+00,0.1148E+00,0.1246E+00,
     *0.1355E+00,0.1474E+00,0.1605E+00,0.1750E+00,0.1910E+00,0.2088E+00,
     *0.2284E+00,0.2501E+00,0.2739E+00,0.3000E+00,0.3284E+00,0.3594E+00,
     *0.3935E+00,0.4308E+00,0.4718E+00,0.5172E+00,0.5670E+00,0.6222E+00,
     *0.6840E+00,0.7528E+00,0.8308E+00,0.9217E+00,0.1029E+01,0.1159E+01,
     *0.1327E+01,0.1553E+01,0.1878E+01,0.2386E+01,0.3253E+01,0.4937E+01,
     *0.8737E+01,0.1952E+02,0.6427E+02,0.4929E+03,0.5169E+05/
            DATA ((PHR(I,J),J=1,83),I=02,02) /
     *0.2467E+00,0.2483E+00,0.2871E+00,0.2722E+00,0.2454E+00,0.2231E+00,
     *0.2060E+00,0.1900E+00,0.1704E+00,0.1452E+00,0.1186E+00,0.9754E-01,
     *0.8182E-01,0.7067E-01,0.6284E-01,0.5734E-01,0.5345E-01,0.5070E-01,
     *0.4875E-01,0.4741E-01,0.4651E-01,0.4596E-01,0.4570E-01,0.4569E-01,
     *0.4589E-01,0.4631E-01,0.4693E-01,0.4776E-01,0.4879E-01,0.5005E-01,
     *0.5153E-01,0.5328E-01,0.5532E-01,0.5768E-01,0.6040E-01,0.6350E-01,
     *0.6704E-01,0.7104E-01,0.7559E-01,0.8071E-01,0.8648E-01,0.8967E-01,
     *0.9298E-01,0.1002E+00,0.1083E+00,0.1173E+00,0.1273E+00,0.1384E+00,
     *0.1507E+00,0.1641E+00,0.1790E+00,0.1954E+00,0.2134E+00,0.2335E+00,
     *0.2557E+00,0.2801E+00,0.3070E+00,0.3366E+00,0.3687E+00,0.4039E+00,
     *0.4427E+00,0.4850E+00,0.5316E+00,0.5834E+00,0.6402E+00,0.7032E+00,
     *0.7738E+00,0.8527E+00,0.9422E+00,0.1047E+01,0.1171E+01,0.1321E+01,
     *0.1516E+01,0.1780E+01,0.2160E+01,0.2753E+01,0.3768E+01,0.5728E+01,
     *0.1011E+02,0.2231E+02,0.7109E+02,0.5001E+03,0.3548E+05/
            DATA ((PHR(I,J),J=1,83),I=03,03) /
     *0.2599E+00,0.2602E+00,0.2986E+00,0.2838E+00,0.2563E+00,0.2330E+00,
     *0.2148E+00,0.1978E+00,0.1774E+00,0.1513E+00,0.1237E+00,0.1017E+00,
     *0.8513E-01,0.7333E-01,0.6499E-01,0.5912E-01,0.5494E-01,0.5198E-01,
     *0.4986E-01,0.4840E-01,0.4742E-01,0.4681E-01,0.4651E-01,0.4647E-01,
     *0.4667E-01,0.4708E-01,0.4772E-01,0.4858E-01,0.4965E-01,0.5094E-01,
     *0.5249E-01,0.5430E-01,0.5642E-01,0.5887E-01,0.6169E-01,0.6491E-01,
     *0.6858E-01,0.7273E-01,0.7744E-01,0.8274E-01,0.8872E-01,0.9201E-01,
     *0.9544E-01,0.1029E+00,0.1113E+00,0.1206E+00,0.1309E+00,0.1424E+00,
     *0.1550E+00,0.1689E+00,0.1842E+00,0.2011E+00,0.2198E+00,0.2404E+00,
     *0.2633E+00,0.2886E+00,0.3163E+00,0.3468E+00,0.3800E+00,0.4164E+00,
     *0.4565E+00,0.5002E+00,0.5485E+00,0.6020E+00,0.6608E+00,0.7261E+00,
     *0.7993E+00,0.8810E+00,0.9739E+00,0.1083E+01,0.1211E+01,0.1368E+01,
     *0.1571E+01,0.1846E+01,0.2242E+01,0.2860E+01,0.3918E+01,0.5956E+01,
     *0.1050E+02,0.2307E+02,0.7281E+02,0.4999E+03,0.3196E+05/
            DATA ((PHR(I,J),J=1,83),I=04,04) /
     *0.2765E+00,0.2752E+00,0.3129E+00,0.2981E+00,0.2697E+00,0.2452E+00,
     *0.2256E+00,0.2075E+00,0.1862E+00,0.1589E+00,0.1301E+00,0.1069E+00,
     *0.8939E-01,0.7677E-01,0.6780E-01,0.6145E-01,0.5690E-01,0.5366E-01,
     *0.5134E-01,0.4973E-01,0.4862E-01,0.4794E-01,0.4758E-01,0.4751E-01,
     *0.4769E-01,0.4811E-01,0.4877E-01,0.4965E-01,0.5076E-01,0.5212E-01,
     *0.5373E-01,0.5563E-01,0.5784E-01,0.6041E-01,0.6336E-01,0.6672E-01,
     *0.7055E-01,0.7488E-01,0.7979E-01,0.8532E-01,0.9155E-01,0.9497E-01,
     *0.9854E-01,0.1063E+00,0.1150E+00,0.1247E+00,0.1354E+00,0.1473E+00,
     *0.1604E+00,0.1748E+00,0.1907E+00,0.2083E+00,0.2276E+00,0.2491E+00,
     *0.2729E+00,0.2990E+00,0.3279E+00,0.3596E+00,0.3941E+00,0.4319E+00,
     *0.4735E+00,0.5191E+00,0.5693E+00,0.6251E+00,0.6864E+00,0.7545E+00,
     *0.8309E+00,0.9163E+00,0.1013E+01,0.1127E+01,0.1262E+01,0.1426E+01,
     *0.1640E+01,0.1928E+01,0.2345E+01,0.2995E+01,0.4106E+01,0.6242E+01,
     *0.1098E+02,0.2400E+02,0.7481E+02,0.4984E+03,0.2810E+05/
            DATA ((PHR(I,J),J=1,83),I=05,05) /
     *0.3140E+00,0.3090E+00,0.3440E+00,0.3291E+00,0.2988E+00,0.2716E+00,
     *0.2491E+00,0.2285E+00,0.2053E+00,0.1759E+00,0.1447E+00,0.1190E+00,
     *0.9926E-01,0.8484E-01,0.7446E-01,0.6700E-01,0.6162E-01,0.5774E-01,
     *0.5493E-01,0.5295E-01,0.5158E-01,0.5070E-01,0.5021E-01,0.5005E-01,
     *0.5019E-01,0.5060E-01,0.5129E-01,0.5224E-01,0.5344E-01,0.5492E-01,
     *0.5668E-01,0.5876E-01,0.6118E-01,0.6400E-01,0.6723E-01,0.7091E-01,
     *0.7509E-01,0.7981E-01,0.8516E-01,0.9117E-01,0.9793E-01,0.1016E+00,
     *0.1055E+00,0.1140E+00,0.1234E+00,0.1338E+00,0.1454E+00,0.1582E+00,
     *0.1724E+00,0.1879E+00,0.2051E+00,0.2241E+00,0.2449E+00,0.2681E+00,
     *0.2937E+00,0.3220E+00,0.3531E+00,0.3873E+00,0.4247E+00,0.4656E+00,
     *0.5108E+00,0.5603E+00,0.6149E+00,0.6756E+00,0.7425E+00,0.8168E+00,
     *0.9003E+00,0.9939E+00,0.1101E+01,0.1226E+01,0.1374E+01,0.1557E+01,
     *0.1793E+01,0.2114E+01,0.2577E+01,0.3299E+01,0.4529E+01,0.6879E+01,
     *0.1204E+02,0.2596E+02,0.7866E+02,0.4906E+03,0.2124E+05/
            DATA ((PHR(I,J),J=1,83),I=06,06) /
     *0.3397E+00,0.3323E+00,0.3646E+00,0.3493E+00,0.3179E+00,0.2889E+00,
     *0.2644E+00,0.2424E+00,0.2181E+00,0.1874E+00,0.1547E+00,0.1274E+00,
     *0.1062E+00,0.9063E-01,0.7928E-01,0.7107E-01,0.6509E-01,0.6076E-01,
     *0.5761E-01,0.5537E-01,0.5380E-01,0.5278E-01,0.5218E-01,0.5196E-01,
     *0.5206E-01,0.5246E-01,0.5317E-01,0.5415E-01,0.5542E-01,0.5697E-01,
     *0.5883E-01,0.6103E-01,0.6359E-01,0.6657E-01,0.6998E-01,0.7387E-01,
     *0.7829E-01,0.8327E-01,0.8891E-01,0.9524E-01,0.1024E+00,0.1063E+00,
     *0.1103E+00,0.1192E+00,0.1291E+00,0.1400E+00,0.1522E+00,0.1656E+00,
     *0.1805E+00,0.1968E+00,0.2148E+00,0.2346E+00,0.2565E+00,0.2807E+00,
     *0.3076E+00,0.3372E+00,0.3699E+00,0.4058E+00,0.4451E+00,0.4881E+00,
     *0.5357E+00,0.5878E+00,0.6454E+00,0.7094E+00,0.7800E+00,0.8586E+00,
     *0.9471E+00,0.1046E+01,0.1160E+01,0.1293E+01,0.1451E+01,0.1646E+01,
     *0.1899E+01,0.2242E+01,0.2738E+01,0.3509E+01,0.4820E+01,0.7310E+01,
     *0.1274E+02,0.2720E+02,0.8080E+02,0.4822E+03,0.1763E+05/
            DATA ((PHR(I,J),J=1,83),I=07,07) /
     *0.3665E+00,0.3585E+00,0.3853E+00,0.3705E+00,0.3386E+00,0.3093E+00,
     *0.2869E+00,0.2705E+00,0.2507E+00,0.2187E+00,0.1832E+00,0.1512E+00,
     *0.1258E+00,0.1065E+00,0.9217E-01,0.8162E-01,0.7386E-01,0.6812E-01,
     *0.6393E-01,0.6088E-01,0.5870E-01,0.5723E-01,0.5631E-01,0.5585E-01,
     *0.5579E-01,0.5612E-01,0.5681E-01,0.5783E-01,0.5918E-01,0.6088E-01,
     *0.6291E-01,0.6532E-01,0.6815E-01,0.7143E-01,0.7521E-01,0.7951E-01,
     *0.8439E-01,0.8988E-01,0.9607E-01,0.1030E+00,0.1108E+00,0.1151E+00,
     *0.1196E+00,0.1293E+00,0.1400E+00,0.1520E+00,0.1652E+00,0.1799E+00,
     *0.1961E+00,0.2140E+00,0.2338E+00,0.2557E+00,0.2799E+00,0.3069E+00,
     *0.3367E+00,0.3696E+00,0.4060E+00,0.4461E+00,0.4901E+00,0.5388E+00,
     *0.5927E+00,0.6520E+00,0.7180E+00,0.7913E+00,0.8725E+00,0.9634E+00,
     *0.1066E+01,0.1181E+01,0.1314E+01,0.1469E+01,0.1655E+01,0.1885E+01,
     *0.2183E+01,0.2586E+01,0.3166E+01,0.4061E+01,0.5568E+01,0.8386E+01,
     *0.1440E+02,0.2992E+02,0.8452E+02,0.4537E+03,0.1132E+05/
            DATA ((PHR(I,J),J=1,83),I=08,08) /
     *0.2248E+00,0.2041E+00,0.2013E+00,0.2015E+00,0.2038E+00,0.2142E+00,
     *0.2218E+00,0.2177E+00,0.2078E+00,0.1973E+00,0.1876E+00,0.1779E+00,
     *0.1666E+00,0.1530E+00,0.1377E+00,0.1221E+00,0.1078E+00,0.9531E-01,
     *0.8504E-01,0.7686E-01,0.7052E-01,0.6573E-01,0.6219E-01,0.5966E-01,
     *0.5794E-01,0.5689E-01,0.5645E-01,0.5656E-01,0.5718E-01,0.5825E-01,
     *0.5974E-01,0.6159E-01,0.6382E-01,0.6647E-01,0.6955E-01,0.7314E-01,
     *0.7723E-01,0.8187E-01,0.8711E-01,0.9302E-01,0.9976E-01,0.1035E+00,
     *0.1075E+00,0.1163E+00,0.1263E+00,0.1377E+00,0.1507E+00,0.1653E+00,
     *0.1819E+00,0.2008E+00,0.2222E+00,0.2467E+00,0.2745E+00,0.3060E+00,
     *0.3418E+00,0.3822E+00,0.4279E+00,0.4800E+00,0.5391E+00,0.6066E+00,
     *0.6838E+00,0.7715E+00,0.8718E+00,0.9864E+00,0.1117E+01,0.1268E+01,
     *0.1442E+01,0.1643E+01,0.1880E+01,0.2160E+01,0.2496E+01,0.2906E+01,
     *0.3423E+01,0.4095E+01,0.5014E+01,0.6356E+01,0.8465E+01,0.1211E+02,
     *0.1924E+02,0.3569E+02,0.8510E+02,0.3357E+03,0.3290E+04/
            DATA ((PHR(I,J),J=1,83),I=09,09) /
     *0.8649E-01,0.6705E-01,0.5195E-01,0.7001E-01,0.7008E-01,0.6002E-01,
     *0.5176E-01,0.4616E-01,0.4241E-01,0.3977E-01,0.3795E-01,0.3668E-01,
     *0.3583E-01,0.3535E-01,0.3514E-01,0.3524E-01,0.3565E-01,0.3638E-01,
     *0.3751E-01,0.3892E-01,0.4055E-01,0.4217E-01,0.4354E-01,0.4447E-01,
     *0.4473E-01,0.4432E-01,0.4334E-01,0.4196E-01,0.4043E-01,0.3895E-01,
     *0.3767E-01,0.3668E-01,0.3599E-01,0.3567E-01,0.3568E-01,0.3603E-01,
     *0.3675E-01,0.3782E-01,0.3929E-01,0.4119E-01,0.4354E-01,0.4489E-01,
     *0.4638E-01,0.4977E-01,0.5377E-01,0.5848E-01,0.6402E-01,0.7052E-01,
     *0.7819E-01,0.8720E-01,0.9780E-01,0.1103E+00,0.1250E+00,0.1423E+00,
     *0.1629E+00,0.1872E+00,0.2164E+00,0.2514E+00,0.2934E+00,0.3442E+00,
     *0.4055E+00,0.4799E+00,0.5709E+00,0.6824E+00,0.8200E+00,0.9912E+00,
     *0.1205E+01,0.1474E+01,0.1814E+01,0.2247E+01,0.2801E+01,0.3520E+01,
     *0.4460E+01,0.5710E+01,0.7406E+01,0.9765E+01,0.1318E+02,0.1847E+02,
     *0.2749E+02,0.4547E+02,0.9155E+02,0.2798E+03,0.1582E+04/
            DATA ((PHR(I,J),J=1,83),I=10,10) /
     *0.9344E-01,0.8261E-01,0.6680E-01,0.7550E-01,0.8962E-01,0.9095E-01,
     *0.8469E-01,0.7755E-01,0.7170E-01,0.6726E-01,0.6401E-01,0.6173E-01,
     *0.6034E-01,0.5974E-01,0.5979E-01,0.6028E-01,0.6096E-01,0.6155E-01,
     *0.6179E-01,0.6151E-01,0.6067E-01,0.5928E-01,0.5752E-01,0.5554E-01,
     *0.5354E-01,0.5165E-01,0.4997E-01,0.4858E-01,0.4752E-01,0.4683E-01,
     *0.4651E-01,0.4657E-01,0.4701E-01,0.4781E-01,0.4897E-01,0.5053E-01,
     *0.5250E-01,0.5493E-01,0.5787E-01,0.6137E-01,0.6550E-01,0.6782E-01,
     *0.7033E-01,0.7593E-01,0.8242E-01,0.8992E-01,0.9860E-01,0.1087E+00,
     *0.1203E+00,0.1339E+00,0.1497E+00,0.1682E+00,0.1896E+00,0.2147E+00,
     *0.2441E+00,0.2786E+00,0.3193E+00,0.3675E+00,0.4248E+00,0.4931E+00,
     *0.5747E+00,0.6726E+00,0.7902E+00,0.9324E+00,0.1105E+01,0.1316E+01,
     *0.1575E+01,0.1895E+01,0.2292E+01,0.2787E+01,0.3407E+01,0.4192E+01,
     *0.5195E+01,0.6498E+01,0.8221E+01,0.1057E+02,0.1389E+02,0.1886E+02,
     *0.2699E+02,0.4205E+02,0.7598E+02,0.1847E+03,0.5926E+03/
      do 1 i=1,10
      do 1 j=1,83
      ph(i,j)=phr(i,j)
    1 continue
      return
      end
      subroutine enviro (difr,difa,r,palt,xmuv,
     a                   fra,fae,fr)
      real difr, difa, r, palt
      real fae,fra,fr,fae0,fra0,xmuv,xlnv,a0,b0,a1,b1
      real zmin,zmax,xcfr1,xcfr2,xcfa1,xcfa2,xcfa3
      real alt(16),cfr1(16),cfr2(16),cfa1(16),cfa2(16),cfa3(16)
      integer i
      data (alt(i),i=1,16) /0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,
     s                      10.0,12.0,14.0,16.0,18.0,20.0,60.0/
      data (cfr1(i),i=1,16) /0.730,0.710,0.656,0.606,0.560,0.516,0.473,
     s 0.433,0.395,0.323,0.258,0.209,0.171,0.142,0.122,0.070/
      data (cfr2(i),i=1,16) /2.8,1.51,0.845,0.634,0.524,0.465,0.429,
     s 0.405,0.390,0.386,0.409,0.445,0.488,0.545,0.608,0.868/
      data (cfa1(i),i=1,16) /0.239,0.396,0.588,0.626,0.612,0.505,0.454,
     s 0.448,0.444,0.445,0.444,0.448,0.448,0.448,0.448,0.448/
      data (cfa2(i),i=1,16) /1.40,1.20,1.02,0.86,0.74,0.56,0.46,0.42,
     s 0.38,0.34,0.3,0.28,0.27,0.27,0.27,0.27/
      data (cfa3(i),i=1,16) /9.17,6.26,5.48,5.16,4.74,3.65,3.24,3.15,
     s 3.07,2.97,2.88,2.83,2.83,2.83,2.83,2.83/
c
c     calculation of the environmental function for
c     rayleigh and aerosols contribution.
c
c     this calculation have been done for nadir observation
c     and are corrected of the effect of the view zenith angle.
c
      a0=1.3347
      b0=0.57757
      a1=-1.479
      b1=-1.5275

      if (palt.ge.60.) then
      fae0=1-0.448*exp(-r*0.27)-0.552*exp(-r*2.83)
      fra0=1-0.930*exp(-r*0.080)-0.070*exp(-r*1.100)
      else
      i=0
 10   i=i+1
      if (palt.ge.alt(i)) goto 10
      if ((i.gt.1).and.(i.lt.16)) then
         zmin=alt(i-1)
         zmax=alt(i)
         xcfr1=cfr1(i-1)+(cfr1(i)-cfr1(i-1))*(zmax-palt)/(zmax-zmin)
         xcfr2=cfr2(i-1)+(cfr2(i)-cfr2(i-1))*(zmax-palt)/(zmax-zmin)
         xcfa1=cfa1(i-1)+(cfa1(i)-cfa1(i-1))*(zmax-palt)/(zmax-zmin)
         xcfa2=cfa2(i-1)+(cfa2(i)-cfa2(i-1))*(zmax-palt)/(zmax-zmin)
         xcfa3=cfa3(i-1)+(cfa3(i)-cfa3(i-1))*(zmax-palt)/(zmax-zmin)
         endif
      if (i.eq.1) then
         xcfr1=cfr1(1)
         xcfr2=cfr2(1)
         xcfa1=cfa1(1)
         xcfa2=cfa2(1)
         xcfa3=cfa3(1)
         endif
      fra0=1.-xcfr1*exp(-r*xcfr2)-(1.-xcfr1)*exp(-r*0.08)
      fae0=1.-xcfa1*exp(-r*xcfa2)-(1.-xcfa1)*exp(-r*xcfa3)
      endif
c correction of the effect of the view zenith angle
      xlnv=log(xmuv)
      fra=fra0*(xlnv*(1-fra0)+1)
      fae=fae0*((1+a0*xlnv+b0*xlnv*xlnv)+fae0*(a1*xlnv+b1*xlnv*xlnv)+
     sfae0*fae0*((-a1-a0)*xlnv+(-b1-b0)*xlnv*xlnv))
c
      if ((difa+difr).gt.1.E-03) then
         fr=(fae*difa+fra*difr)/(difa+difr)
         else
         fr=1.
         endif
      return
      end
      subroutine equivwl(iinf,isup,step,
     s               wlmoy)
      common /sixs_ffu/s(1501),wlinf,wlsup
      real step,wlmoy,s,wlinf,wlsup,seb,wlwave,sbor,wl,swl,coef
      integer iinf,isup,l
      seb=0.
      wlwave=0.
      do 50 l=iinf,isup
      sbor=s(l)
      if(l.eq.iinf.or.l.eq.isup) sbor=sbor*0.5
      wl=.25+(l-1)*step
C---      call solirr(wl,
C---     s            swl)
      swl = 1.0
C---
      coef=sbor*step*swl
      seb=seb+coef
      wlwave=wlwave+wl*coef
  50  continue
      wlmoy=wlwave/seb
      return
      end
      subroutine gauss(x1,x2,x,w,n)
      integer n
      real x1,x2,x(n),w(n)
      double precision xm,xl,z,p1,p2,p3,pp,z1
      integer m,i,j
      parameter (eps=3.d-14)
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.eps)go to 1
        if (abs(z).lt.eps) z=0.
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      end
      subroutine interp (iaer,idatmp,wl,taer55,taer55p,xmud,
     a                   romix,rorayl,roaero,phaa,phar,tsca,
     a                   tray,trayp,taer,taerp,dtott,utott,
     a                   astot,asray,asaer,
     a                   utotr,utota,dtotr,dtota)
 
      common /sixs_aer/ext(10),ome(10),gasym(10),phase(10)
      common /sixs_disc/ roatm(3,10),dtdir(3,10),dtdif(3,10),
     a utdir(3,10),utdif(3,10),sphal(3,10),wldis(10),trayl(10),
     a traypl(10)
      common /sixs_del/ delta,sigma
      Real wl,taer55,taer55p
      Real xmud,romix,rorayl,roaero,phaa,phar,tsca,tray
      Real trayp,taer,taerp,dtott,utott,astot,asray,asaer,utotr
      Real utota,dtotr,dtota,ext,ome,gasym,phase,roatm,dtdir
      Real dtdif,utdir,utdif,sphal,wldis,trayl,traypl,delta,sigma
      Real alphaa,betaa,alphar,betar,alphac,betac,coef,wlinf,d2
      Real drinf,drsup,dtinf,dtsup,dtotc,dainf,dasup,urinf,ursup
      Real utinf,utsup,utotc,uainf,uasup,arinf,arsup,atinf,atsup
      Real aainf,aasup
      Integer iaer,idatmp,linf,ll,lsup

 
c     that for the atmosphere :
c     the reflectances
c                     rayleigh                             = rorayl
c                     aerosols                             = roaero
c                     mixing                               = romix
c     the downward transmittances
c                     rayleigh                             = dtotr
c                     aerosols                             = dtota
c                     total                                = dtott
c     the upward transmittances
c                     rayleigh                             = utotr
c                     aerosols                             = utota
c                     total                                = utott
c     the spherical albedos
c                     rayleigh                             = asray
c                     aerosols                             = asaer
c                     total                                = astot
c     the optical thickness of total atmosphere
c                     rayleigh                             = tray
c                     aerosols                             = taer
c     the optical thickness of the atmosphere above the plane
c                     rayleigh                             = trayp
c                     aerosols                             = taerp
c     the tsca of the aerosols (god dammed it)
c                     total atmosphere                     = tsca
      
      linf=1
      do 81 ll=1,9
      if(wl.gt.wldis(ll).and.wl.le.wldis(ll+1)) linf=ll
   81 continue
      if(wl.gt.wldis(10)) linf=9
      lsup=linf+1
 
c     interpolation in function of wavelength for scattering
c     atmospheric functions from discrete values at wldis
 
      alphaa=0.
      betaa=0.
      alphar=0.
      betar=0.
      alphac=0.
      betac=0.
      phaa=0.
      roaero=0.
      dtota=1.
      utota=1.
      asaer=0.
      taer=0.
      taerp=0.
      coef=alog(wldis(lsup)/wldis(linf))
      wlinf=wldis(linf)
c
      if(iaer.eq.0) goto 1240
      alphaa=alog(phase(lsup)/phase(linf))/coef
      betaa=phase(linf)/(wlinf**(alphaa))
      phaa=betaa*(wl**alphaa)
 1240 d2=2.+delta
      phar=(2.*(1.-delta)/d2)*.75*(1.+xmud*xmud)+3.*delta/d2
      if (idatmp.eq.0) then
         betar=0.
         betaa=0.
         betac=0.
         goto 1234
      endif
      if(roatm(1,linf).lt..001) then
	rorayl=roatm(1,linf)+(roatm(1,lsup)-roatm(1,linf))
     s     *(wl-wldis(linf))/(wldis(lsup)-wldis(linf))
	else
        alphar=alog(roatm(1,lsup)/roatm(1,linf))/ coef
        betar=roatm(1,linf)/(wlinf**(alphar))
	rorayl=betar*(wl**alphar)
      endif
      if(roatm(2,linf).lt..001) then
        romix=roatm(2,linf)+(roatm(2,lsup)-roatm(2,linf))
     s     *(wl-wldis(linf))/(wldis(lsup)-wldis(linf))
        else
        alphac=alog(roatm(2,lsup)/roatm(2,linf))/coef
        betac=roatm(2,linf)/(wlinf**(alphac))
	romix=betac*(wl**alphac)
      endif
      if(iaer.eq.0) goto 1234
      if(roatm(3,linf).lt..001) then
	roaero=roatm(3,linf)+(roatm(3,lsup)-roatm(3,linf))
     s     *(wl-wldis(linf))/(wldis(lsup)-wldis(linf))
	else
        alphaa=alog(roatm(3,lsup)/roatm(3,linf))/coef
        betaa=roatm(3,linf)/(wlinf**(alphaa))
        roaero=betaa*(wl**alphaa)
      endif
 1234 continue
c
      alphar=alog(trayl(lsup)/trayl(linf))/coef
      betar=trayl(linf)/(wlinf**(alphar))
      tray=betar*(wl**alphar)
      if (idatmp.ne.0.) then
         alphar=alog(traypl(lsup)/traypl(linf))/coef
         betar=traypl(linf)/(wlinf**(alphar))
         trayp=betar*(wl**alphar)
         else
         trayp=0.
         endif
c
      if(iaer.eq.0) goto 1235
      alphaa=alog(ext(lsup)*ome(lsup)/(ext(linf)*ome(linf)))/coef
      betaa=ext(linf)*ome(linf)/(wlinf**(alphaa))
      tsca=taer55*betaa*(wl**alphaa)/ext(4)
      alphaa=alog(ext(lsup)/ext(linf))/coef
      betaa=ext(linf)/(wlinf**(alphaa))
      taerp=taer55p*betaa*(wl**alphaa)/ext(4)
      taer=taer55*betaa*(wl**alphaa)/ext(4)
c
 1235 drinf=dtdif(1,linf)+dtdir(1,linf)
      drsup=dtdif(1,lsup)+dtdir(1,lsup)
      alphar=alog(drsup/drinf)/coef
      betar=drinf/(wlinf**(alphar))
      dtotr=betar*(wl**alphar)
      dtinf=dtdif(2,linf)+dtdir(2,linf)
      dtsup=dtdif(2,lsup)+dtdir(2,lsup)
      alphac=alog((dtsup*drinf)/(dtinf*drsup))/coef
      betac=(dtinf/drinf)/(wlinf**(alphac))
      dtotc=betac*(wl**alphac)
      dainf=dtdif(3,linf)+dtdir(3,linf)
      dasup=dtdif(3,lsup)+dtdir(3,lsup)
      if(iaer.eq.0) goto 1236
      alphaa=alog(dasup/dainf)/coef
      betaa=dainf/(wlinf**(alphaa))
      dtota=betaa*(wl**alphaa)
 1236 dtott=dtotc*dtotr
      urinf=utdif(1,linf)+utdir(1,linf)
      ursup=utdif(1,lsup)+utdir(1,lsup)
      alphar=alog(ursup/urinf)/ coef
      betar=urinf/(wlinf**(alphar))
      utotr=betar*(wl**alphar)
      utinf=utdif(2,linf)+utdir(2,linf)
      utsup=utdif(2,lsup)+utdir(2,lsup)
      alphac=alog((utsup*urinf)/(utinf*ursup))/ coef
      betac=(utinf/urinf)/(wlinf**(alphac))
      utotc=betac*(wl**alphac)
      uainf=utdif(3,linf)+utdir(3,linf)
      uasup=utdif(3,lsup)+utdir(3,lsup)
      if(iaer.eq.0) goto 1237
      alphaa=alog(uasup/uainf)/ coef
      betaa=uainf/(wlinf**(alphaa))
      utota=betaa*(wl**alphaa)
 1237 utott=utotc*utotr
      arinf=sphal(1,linf)
      arsup=sphal(1,lsup)
      alphar=alog(arsup/arinf)/ coef
      betar=arinf/(wlinf**(alphar))
      asray=betar*(wl**alphar)
      atinf=sphal(2,linf)
      atsup=sphal(2,lsup)
      alphac=alog(atsup/atinf)/coef
      betac=atinf/(wlinf**(alphac))
      astot=betac*(wl**alphac)
      aainf=sphal(3,linf)
      aasup=sphal(3,lsup)
      if(iaer.eq.0) goto 1239
      alphaa=alog(aasup/aainf)/coef
      betaa=aainf/(wlinf**(alphaa))
      asaer=betaa*(wl**alphaa)
 1239 return
      end
      subroutine iso(tamoy,trmoy,pizmoy,tamoyp,trmoyp,palt,
     s               nt,mu,rm,gb,
     s                     xf)
c  dimension for gauss integration
      integer mu
      real rm(-mu:mu),gb(-mu:mu)
c  dimension for os computation
      real  xf(-1:1)
c array for sos computation
      real xpl(-25:25),psl(-1:80,-25:25),bp(0:25,-25:25),
     s xdel(0:30),ydel(0:30),h(0:30)
      real i1(0:30,-25:25),i2(0:30,-25:25),i3(-25:25),
     s   in(-25:25),inm1(-25:25),inm2(-25:25)
      real altc(0:30)
      real ii1,ii2
      Real tamoy,trmoy,pizmoy
      Real tamoyp,trmoyp,palt
      Real delta,sigma,pha,betal,accu,accu2,ta,piz
      Real tr,trp,tap,hr,ha,zx,yy,dd,ppp2,ppp1,ca
      Real cr,ratio,taup,th,xt1,xt2,aaaa,ron,beta0,beta2
      Real tavion0,tavion1,tavion2,tavion,zi1,xpk,ypk,x,y,xpj
      Real z,xi1,xi2,bpjk,bpjmk,f,a,b,c,d,xx,a1,d1,g1
      Real y1
      Double precision xxx
      integer snt,nt,iplane,ntp,j,it,itp,i,ig,k,index,iwr,m
      integer jj,l
      logical ier
      common/sixs_del/delta,sigma
      common/sixs_trunc/pha(83),betal(0:80)
      common/sixs_ier/iwr,ier
      snt=nt
      iplane=0
      accu=1.e-20
      accu2=1.e-3
      ta=tamoy
      piz=pizmoy
      tr=trmoy
      do 615 m=-1,1
 615  xf(m)=0.
c
c     molecular ratio within the layer
c     computations are performed assuming a scale of 8km for
c     molecules and 2km for aerosols
c
c the optical thickness above plane are recomputed to give o.t above pla
      trp=trmoy-trmoyp
      tap=tamoy-tamoyp
c     print *, 'tamoy,trmoy,pizmoy,tap,trp,palt,nt'
c     print *,tamoy,trmoy,pizmoy,tap,trp,palt,nt
      accu=1.e-20
c if plane observations recompute scale height for aerosol knowing:
c the aerosol optical depth as measure from the plane 	= tamoyp
c the rayleigh scale   height = 			= hr (8km)
c the rayleigh optical depth  at plane level 		= trmoyp
c the altitude of the plane 				= palt
c the rayleigh optical depth for total atmos		= trmoy
c the aerosol  optical depth for total atmos		= tamoy
c if not plane observations then ha is equal to 2.0km
c ntp local variable: if ntp=nt     no plane observation selected
c                        ntp=nt-1   plane observation selected
      hr=8.0
c     it's a mixing rayleigh+aerosol
      if(palt.le.900..and.palt.gt.0.0)then
      if (tap.gt.1.e-03) then
         ha=-palt/log(tap/ta)
         else
         ha=2.
         endif
      ntp=nt-1
      else
      ha=2.0
      ntp=nt
      endif
c
      ta=tamoy
      tr=trmoy
      piz=pizmoy
c
c compute mixing rayleigh, aerosol
c case 1: pure rayleigh
c case 2: pure aerosol
c case 3: mixing rayleigh-aerosol
c
      if((ta.le.accu2).and.(tr.gt.ta)) then
      do j=0,ntp
      h(j)=j*tr/ntp
      ydel(j)=1.0
      xdel(j)=0.0
      enddo
      endif
      if((tr.le.accu2).and.(ta.gt.tr)) then
      do j=0,ntp
      h(j)=j*ta/ntp
      ydel(j)=0.0
      xdel(j)=piz
      enddo
      endif
c
      if(tr.gt.accu2.and.ta.gt.accu2)then
      ydel(0)=1.0
      xdel(0)=0.0
      h(0)=0.
      altc(0)=300.
      zx=300.
      iplane=0
      do 14 it=0,ntp
      if (it.eq.0) then
         yy=0.
         dd=0.
         goto 111
      endif
      yy=h(it-1)
      dd=ydel(it-1)
 111  ppp2=300.0
      ppp1=0.0
      itp=it
      call discre(ta,ha,tr,hr,itp,ntp,yy,dd,ppp2,ppp1,
     s    zx)
      if(ier)return
      xxx=-zx/ha
      if (xxx.lt.-18) then
         ca=0.
         else
         ca=ta*dexp(xxx)
         endif
      xxx=-zx/hr
      cr=tr*dexp(xxx)
      h(it)=cr+ca
      altc(it)=zx
c     print *,it,cr,ca,h(it),zx
      cr=cr/hr
      ca=ca/ha
      ratio=cr/(cr+ca)
      xdel(it)=(1.e+00-ratio)*piz
      ydel(it)=ratio
  14  continue
      endif
c update plane layer if necessary
      if (ntp.eq.(nt-1)) then
c compute position of the plane layer
         taup=tap+trp
         iplane=-1
         do i=0,ntp
         if (taup.ge.h(i)) iplane=i
         enddo
c update the layer from the end to the position to update if necessary
         th=0.0005
         xt1=abs(h(iplane)-taup)
         xt2=abs(h(iplane+1)-taup)
         if ((xt1.gt.th).and.(xt2.gt.th)) then
         do i=nt,iplane+1,-1
            xdel(i)=xdel(i-1)
            ydel(i)=ydel(i-1)
            h(i)=h(i-1)
            altc(i)=altc(i-1)
         enddo
         else
         nt=ntp
         if (xt2.lt.xt1) iplane=iplane+1
         endif
         h(iplane)=taup
         if ( tr.gt.accu2.and.ta.gt.accu2) then
         ca=ta*exp(-palt/ha)
         cr=tr*exp(-palt/hr)
         cr=cr/hr
         ca=ca/ha
         ratio=cr/(cr+ca)
         xdel(iplane)=(1.e+00-ratio)*piz
         ydel(iplane)=ratio
         altc(iplane)=palt
         endif
         if ( tr.gt.accu2.and.ta.le.accu2) then
         ydel(iplane)=1.
         xdel(iplane)=0.
         altc(iplane)=palt
         endif
         if ( tr.le.accu2.and.ta.gt.accu2) then
         ydel(iplane)=0.
         xdel(iplane)=1.*piz
         altc(iplane)=palt
         endif
      endif
c
c     print *,ha,hr,palt,tamoy,trmoy,tap,trp
c     do i=0,nt
c     print *,i,h(i),xdel(i),ydel(i),altc(i)
c     enddo
c
      aaaa=delta/(2-delta)
      ron=(1-aaaa)/(1+2*aaaa)
c
c     rayleigh phase function
c
      beta0=1.
      beta2=0.5*ron
c
c    primary scattering
c
      ig=1
      tavion0=0.
      tavion1=0.
      tavion2=0.
      tavion=0.
      do 16 j=-mu,mu
      i3(j)=0.
   16 continue
c
c     kernel computations
c
      call kernel(0,mu,rm,xpl,psl,bp)
      do 100 j=-mu,mu
      do 101 k=0,nt
      i2(k,j)=0.0000
  101 continue
  100 continue
c
c     vertical integration, primary upward radiation
c
 
      do 108 k=1,mu
      i1(nt,k)=1.0
      zi1=i1(nt,k)
      yy=rm(k)
      do 108 i=nt-1,0,-1
      i1(i,k)=exp(-(ta+tr-h(i))/yy)
  108 continue
c
c     vertical integration, primary downward radiation
c
      do 109 k=-mu,-1
      do 109 i=0,nt
      i1(i,k)=0.00
  109 continue
c
c     inm2 is inialized with scattering computed at n-2
c     i3 is inialized with primary scattering
c
      do 20 k=-mu,mu
      if(k) 21,20,23
   21 index=nt
      go to 25
   23 index=0
   25 continue
      inm1(k)=i1(index,k)
      inm2(k)=i1(index,k)
      i3(k)=i1(index,k)
   20 continue
      tavion=i1(iplane,mu)
      tavion2=i1(iplane,mu)
c
c     loop on successive order
c
  503 ig=ig+1
c     write(6,*) 'ig ',ig
c
c     successive orders
c
c     multiple scattering source function at every level within the laye
c
c
      do455 k=1,mu
      xpk=xpl(k)
      ypk=xpl(-k)
      do 455 i=0,nt
      ii1=0.
      ii2=0.
      x=xdel(i)
      y=ydel(i)
      do477 j=1,mu
      xpj=xpl(j)
      z=gb(j)
      xi1=i1(i,j)
      xi2=i1(i,-j)
      bpjk=bp(j,k)*x+y*(beta0+beta2*xpj*xpk)
      bpjmk=bp(j,-k)*x+y*(beta0+beta2*xpj*ypk)
      ii2=ii2+z*(xi1*bpjk+xi2*bpjmk)
      ii1=ii1+z*(xi1*bpjmk+xi2*bpjk)
 477  continue
      i2(i,k)=ii2
      i2(i,-k)=ii1
 455  continue
c
c     vertical integration, upward radiation
c
      do 48 k=1,mu
      i1(nt,k)=0.0
      zi1=i1(nt,k)
      yy=rm(k)
      do 48 i=nt-1,0,-1
      jj=i+1
      f=h(jj)-h(i)
      a=(i2(jj,k)-i2(i,k))/f
      b=i2(i,k)-a*h(i)
      c=exp(-f/yy)
      d=1.e+00-c
      xx=h(i)-h(jj)*c
      zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5e+00
      i1(i,k)=zi1
   48 continue
c
c     vertical integration, downward radiation
c
      do 50 k=-mu,-1
      i1(0,k)=0.
      zi1=i1(0,k)
      yy=rm(k)
      do 50 i=1,nt
      jj=i-1
      f=h(i)-h(jj)
      c=exp(f/yy)
      d=1.e+00-c
      a=(i2(i,k)-i2(jj,k))/f
      b=i2(i,k)-a*h(i)
      xx=h(i)-h(jj)*c
      zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5e+00
      i1(i,k)=zi1
   50 continue
c
c     in is the nieme scattering order
c
      do 30 k=-mu,mu
      if(k) 31,30,33
   31 index=nt
      go to 34
   33 index=0
   34 continue
      in(k)=i1(index,k)
   30 continue
      tavion0=i1(iplane,mu)
c
c   convergence test (geometrical serie)
c
      if(ig.gt.2) then
      z=0.
      a1=tavion2
      d1=tavion1
      g1=tavion0
      if (a1.ge.accu.and.d1.ge.accu.and.tavion.ge.accu) then
         y=((g1/d1-d1/a1)/((1.-g1/d1)**2)*(g1/tavion))
         y=abs(y)
         z=amax1(y,z)
      endif
      do 99 l=-mu,mu
      if (l.eq.0) goto 99
      a1=inm2(l)
      d1=inm1(l)
      g1=in(l)
      if(a1.eq.0.) go to 99
      if(d1.eq.0.) go to 99
      if(i3(l).eq.0.) go to 99
      y=((g1/d1-d1/a1)/((1-g1/d1)**2)*(g1/i3(l)))
      y=abs(y)
      z=amax1(y,z)
  99  continue
      if(z.lt.0.0001) then
c
c     successful test (geometrical serie)
c
      do 606 l=-mu,mu
      if (l.eq.0) goto 606
      y1=1.
      d1=inm1(l)
      g1=in(l)
      if(d1.eq.0.0) go to 606
      y1=1-g1/d1
      g1=g1/y1
      i3(l)=i3(l)+g1
  606 continue
      d1=tavion1
      g1=tavion0
      y1=1.
      if (d1.ge.accu) then
      if (abs(g1-d1).ge.accu) then
         y1=1.-g1/d1
         g1=g1/y1
      endif
      tavion=tavion+g1
      endif
      go to 505
      endif
c
c     inm2 is the (n-2)ieme scattering order
c
      do 26 k=-mu,mu
      inm2(k)=inm1(k)
   26 continue
      tavion2=tavion1
      endif
c
c     inm1 is the (n-1)ieme scattering order
c
      do 27 k=-mu,mu
      inm1(k)=in(k)
   27 continue
      tavion1=tavion0
c
c     sum of the n-1 orders
c
      do 610 l=-mu,mu
      i3(l)=i3(l)+in(l)
  610 continue
      tavion=tavion+tavion0
c
c     stop if order n is less than 1% of the sum
c
      z=0.
      do 611 l=-mu,mu
      if(i3(l).ne.0)then
      y=abs(in(l)/i3(l))
      z=amax1(z,y)
      endif
  611 continue
      if(z.lt.0.00001) go to 505
c
c      stop if order n is greater than 20 in any case
c
      if(ig-20) 503,503,505
  505 continue
c
c
      xf(1)=xf(1)+i3(mu)
      xf(-1)=tavion
      do k=1,mu
      xf(0)=xf(0)+rm(k)*gb(k)*i3(-k)
      enddo
      nt=snt
      return
      end
      subroutine kernel(is,mu,rm,xpl,psl,bp)
      integer mu
      real rm(-mu:mu)
      real psl(-1:80,-25:25),xpl(-25:25),bp(0:25,-25:25)
      real pha,betal
      integer is,ip1,j,i,k,ip,ig,l,lp,lm,ij
      double precision xdb,a,b,c,xx,rac3,x,bt,sbp
      common /sixs_trunc/pha(83),betal(0:80)
      ip1=80
      rac3=dsqrt(3.D+00)
      if(is.ne.0)go to 700
      do 25 j=0,mu
      c=dble(rm(j))
      psl(0,-j)=1.
      psl(0,j)=1.
      psl(1,j)=c
      psl(1,-j)=-c
      xdb=(3.*c*c-1.)*0.5
      if (abs(xdb).lt.1.E-30) xdb=0.0
      psl(2,-j)=xdb
      psl(2,j)=xdb
   25 continue
      psl(1,0)=rm(0)
      goto 501
c
  700 if(is.ne.1)go to 701
      do 26 j=0,mu
      c=dble(rm(j))
      x=1.-c*c
      psl(0,j)=0.
      psl(0,-j)=0.
      psl(1,-j)=sqrt(x*0.5)
      psl(1,j)=sqrt(x*0.5)
      psl(2,j)=c*psl(1,j)*rac3
      psl(2,-j)=-psl(2,j)
   26 continue
      psl(2,0)=-psl(2,0)
      goto 501
c
  701 a=1
      do 27 i=1,is
      x=i
      a=a*sqrt((i+is)/x)*0.5
 27   continue
      b=a*sqrt(is/(is+1.))*sqrt((is-1.)/(is+2.))
      do 28 j=0,mu
      c=dble(rm(j))
      xx=1.-c*c
      psl(is-1,j)=0.
      xdb=a*xx**(is*0.5)
      if (abs(xdb).lt.1.E-30) xdb=0.0
      psl(is,-j)=xdb
      psl(is,j)=xdb
   28 continue
  501 k=2
      ip=ip1
      if(is.gt.2)k=is
      if(k.eq.ip)goto 502
      ig=-1
      if(is.eq.1)ig=1
      do 30 l=k,ip-1
      lp=l+1
      lm=l-1
      a=(2*l+1.)/sqrt((l+is+1.)*(l-is+1.))
      b=sqrt(float((l+is)*(l-is)))/(2.*l+1.)
      do 31 j=0,mu
      c=dble(rm(j))
      xdb=a*(c*psl(l,j)-b*psl(lm,j))
      if (abs(xdb).lt.1.E-30) xdb=0.
      psl(lp,j)=xdb
      if(j.eq.0) go to 31
      psl(lp,-j)=ig*psl(lp,j)
   31 continue
      ig=-ig
   30 continue
  502 continue
      do 1005 j=-mu,mu
      xpl(j)=psl(2,j)
 1005 continue
      ij=ip1
      do 32 j=0,mu
      do 32 k=-mu,mu
      sbp=0.
      if(is.gt.ij) goto 1
      do 33 l=is,ij
      bt=betal(l)
      sbp=sbp+dble(psl(l,j))*psl(l,k)*bt
  33  continue
 1    continue
      if (abs(sbp).lt.1.E-30) sbp=0.
      bp(j,k)=sbp
   32 continue
      return
      end
      subroutine   midsum
      common /sixs_atm/z(34),p(34),t(34),wh(34),wo(34)
      real z2(34),p2(34),t2(34),wh2(34),wo2(34)
      real z,p,t,wh,wo
      integer i
c
c     model: midlatitude summer mc clatchey
c
      data(z2(i),i=1, 34)/
     1    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,
     2    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,
     3   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,
     4   35.,   40.,   45.,   50.,   70.,  100.,99999./
      data (p2(i),i=1,34) /
     a1.013e+03,9.020e+02,8.020e+02,7.100e+02,6.280e+02,5.540e+02,
     a4.870e+02,4.260e+02,3.720e+02,3.240e+02,2.810e+02,2.430e+02,
     a2.090e+02,1.790e+02,1.530e+02,1.300e+02,1.110e+02,9.500e+01,
     a8.120e+01,6.950e+01,5.950e+01,5.100e+01,4.370e+01,3.760e+01,
     a3.220e+01,2.770e+01,1.320e+01,6.520e+00,3.330e+00,1.760e+00,
     a9.510e-01,6.710e-02,3.000e-04,0.000e+00/
      data (t2(i),i=1,34) /
     a2.940e+02,2.900e+02,2.850e+02,2.790e+02,2.730e+02,2.670e+02,
     a2.610e+02,2.550e+02,2.480e+02,2.420e+02,2.350e+02,2.290e+02,
     a2.220e+02,2.160e+02,2.160e+02,2.160e+02,2.160e+02,2.160e+02,
     a2.160e+02,2.170e+02,2.180e+02,2.190e+02,2.200e+02,2.220e+02,
     a2.230e+02,2.240e+02,2.340e+02,2.450e+02,2.580e+02,2.700e+02,
     a2.760e+02,2.180e+02,2.100e+02,2.100e+02/
      data (wh2(i),i=1,34) /
     a1.400e+01,9.300e+00,5.900e+00,3.300e+00,1.900e+00,1.000e+00,
     a6.100e-01,3.700e-01,2.100e-01,1.200e-01,6.400e-02,2.200e-02,
     a6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,
     a5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,
     a6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,
     a1.300e-06,1.400e-07,1.000e-09,0.000e+00/
      data (wo2(i),i=1,34) /
     a6.000e-05,6.000e-05,6.000e-05,6.200e-05,6.400e-05,6.600e-05,
     a6.900e-05,7.500e-05,7.900e-05,8.600e-05,9.000e-05,1.100e-04,
     a1.200e-04,1.500e-04,1.800e-04,1.900e-04,2.100e-04,2.400e-04,
     a2.800e-04,3.200e-04,3.400e-04,3.600e-04,3.600e-04,3.400e-04,
     a3.200e-04,3.000e-04,2.000e-04,9.200e-05,4.100e-05,1.300e-05,
     a4.300e-06,8.600e-08,4.300e-11,0.000e+00/
      do 1 i=1,34
      z(i)=z2(i)
      p(i)=p2(i)
      t(i)=t2(i)
      wh(i)=wh2(i)
      wo(i)=wo2(i)
    1 continue
      return
      end
      subroutine   midwin
      common /sixs_atm/z(34),p(34),t(34),wh(34),wo(34)
      real z3(34),p3(34),t3(34),wh3(34),wo3(34)
      real z,p,t,wh,wo
      integer i
c
c     model: midlatitude winter mc clatchey
c
      data(z3(i),i=1, 34)/
     1    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,
     2    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,
     3   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,
     4   35.,   40.,   45.,   50.,   70.,  100.,99999./
      data (p3(i),i=1,34) /
     a1.018e+03,8.973e+02,7.897e+02,6.938e+02,6.081e+02,5.313e+02,
     a4.627e+02,4.016e+02,3.473e+02,2.992e+02,2.568e+02,2.199e+02,
     a1.882e+02,1.610e+02,1.378e+02,1.178e+02,1.007e+02,8.610e+01,
     a7.350e+01,6.280e+01,5.370e+01,4.580e+01,3.910e+01,3.340e+01,
     a2.860e+01,2.430e+01,1.110e+01,5.180e+00,2.530e+00,1.290e+00,
     a6.820e-01,4.670e-02,3.000e-04,0.000e+00/
      data (t3(i),i=1,34) /
     a2.722e+02,2.687e+02,2.652e+02,2.617e+02,2.557e+02,2.497e+02,
     a2.437e+02,2.377e+02,2.317e+02,2.257e+02,2.197e+02,2.192e+02,
     a2.187e+02,2.182e+02,2.177e+02,2.172e+02,2.167e+02,2.162e+02,
     a2.157e+02,2.152e+02,2.152e+02,2.152e+02,2.152e+02,2.152e+02,
     a2.152e+02,2.152e+02,2.174e+02,2.278e+02,2.432e+02,2.585e+02,
     a2.657e+02,2.307e+02,2.102e+02,2.100e+02/
      data (wh3(i),i=1,34) /
     a3.500e+00,2.500e+00,1.800e+00,1.200e+00,6.600e-01,3.800e-01,
     a2.100e-01,8.500e-02,3.500e-02,1.600e-02,7.500e-03,6.900e-03,
     a6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,
     a5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,
     a6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,
     a6.300e-06,1.400e-07,1.000e-09,0.000e+00/
      data (wo3(i),i=1,34) /
     a6.000e-05,5.400e-05,4.900e-05,4.900e-05,4.900e-05,5.800e-05,
     a6.400e-05,7.700e-05,9.000e-05,1.200e-04,1.600e-04,2.100e-04,
     a2.600e-04,3.000e-04,3.200e-04,3.400e-04,3.600e-04,3.900e-04,
     a4.100e-04,4.300e-04,4.500e-04,4.300e-04,4.300e-04,3.900e-04,
     a3.600e-04,3.400e-04,1.900e-04,9.200e-05,4.100e-05,1.300e-05,
     a4.300e-06,8.600e-08,4.300e-11,0.000e+00/
      do 1 i=1,34
      z(i)=z3(i)
      p(i)=p3(i)
      t(i)=t3(i)
      wh(i)=wh3(i)
      wo(i)=wo3(i)
    1 continue
      return
      end
      subroutine mie(iaer,wldis,ex,sc,asy)

      double precision nr,p11(83),p1(10,4,83),ext(10,4),sca(10,4),np(4)
      double precision pi,r,rmind,rmaxd,r0,alpha,dr,xndpr2,Qext,Qsca
      double precision rlogpas
      real ex(4,10),sc(4,10),asy(4,10),wldis(10)
      real phasel,cgaus,pdgs,rmax,rmin,rn,ri,x1,x2,x3,rsunph,nrsunph
      real asy_n,asy_d,cij,ph
      integer nbmu,icp,i,j,l,k,iaer,irsunph
      double precision arg,ldexp
      
      common /sixs_sos/ phasel(10,83),cgaus(83),pdgs(83)
      common /mie_in/ rmax,rmin,icp,rn(10,4),ri(10,4),x1(4),x2(4),
     s x3(4),cij(4),irsunph,rsunph(50),nrsunph(50)
      common /sixs_aerbas/ ph(10,83)

      ldexp=-300.
      pi=4.D+00*datan(1.D+00)
      rlogpas=0.030
      nbmu=83
      do i=1,icp
        np(i)=0.D+00
        do l=1,10
          ex(i,l)=0.0
          sc(i,l)=0.0
          asy(i,l)=0.0
          ext(l,i)=0.D+00
          sca(l,i)=0.D+00
          do k=1,nbmu
            p1(l,i,k)=0.D+00                         
          enddo
        enddo
      enddo
      rmaxd=dble(rmax)
      rmind=dble(rmin)

c LOOPS ON THE NUMBER OF PARTICLE TYPE (4 max)
      do 600 i=1,icp
       r=rmind
       dr=r*(10**rlogpas-1.D+00)
 123   continue
C LOOPS ON THE RADIUS OF THE PARTICLE     

c call of the size distribution nr. For our computation, we need dn/dr for
c all functions except for sun-photometer inputs for which we need dV/dlog(r)
       goto(300,301,302,303)iaer-7
C --- Mixture of particles (Log-Normal distribution functions, up to 5)
 300   nr=DEXP(-5.D-01*(DLOG10(r/x1(i))/DLOG10(1.D+00*x2(i)))**2.D+00)
       nr=nr/dsqrt(2.D+00*pi)/DLOG10(1.D+00*x2(i))
       nr=nr/DLOG(10.D+00)/r
       goto 399

c --- Modified Gamma distribution function
 301   r0=1.00D+00   
       arg=-x2(i)*((r/r0)**x3(i))
       if (arg.gt.ldexp) then
          nr=((r/r0)**x1(i))*DEXP(arg)
          else
          nr=0.
          endif
       goto 399

C --- Junge power-law function
 302   r0=0.1000D+00
       nr= r0**(-x1(i))
       IF(r.GT.r0 ) nr= r**(-x1(i))
       goto 399
C
C --- from sun photometer
 303    nr=0.D+00
	do 299 j=2,irsunph
	if ((r-rsunph(j)).lt.0.000001)then
         nr=(r-rsunph(j-1))/(rsunph(j)-rsunph(j-1))
         nr=nrsunph(j-1)+nr*(nrsunph(j)-nrsunph(j-1))
	 goto 399
	endif
 299   continue
C
c The Mie's calculations have to be called several times (min=2, max=10 for
c each type of particle): at wavelengths bounding the range of the selected
c wavelengths,and at 0.550 microns to normalized the extinction coefficient 
c (if it's not in the selected range of wavelengths).
 399   continue
       xndpr2=nr*dr*pi*(r**2.D+00)
c relatif number of particle for each type of particle (has to be equal to 1)
       np(i)=np(i)+nr*dr
       do l=1,10

         if ((xndpr2*cij(i)).lt.(1.D-08/sqrt(wldis(l))))goto 599

	 alpha=2.D+00*pi*r/wldis(l)
         call EXSCPHASE(alpha,rn(l,i),ri(l,i),Qext,Qsca,p11)
         ext(l,i)=ext(l,i)+xndpr2*Qext
         sca(l,i)=sca(l,i)+xndpr2*Qsca
c phase function for each type of particle
         do k=1,nbmu
          p1(l,i,k)=p1(l,i,k)+p11(k)*xndpr2
         enddo
       enddo
  599  continue
       r=r+dr
       dr=r*(10**rlogpas-1.D+00)
       if(r.ge.rmaxd) goto 600
       goto 123
  600 continue

     
c NOW WE MIXTE THE DIFFERENT TYPES OF PARTICLE
c computation of the scattering and extinction coefficients. We first start
c at 0.550 micron (the extinction coefficient is normalized at 0.550 micron)
      do l=1,10
	do i=1,icp
          ext(l,i)=ext(l,i)/np(i)/1.D+03
          sca(l,i)=sca(l,i)/np(i)/1.D+03
          ex(1,l)=ex(1,l)+cij(i)*real(ext(l,i))
          sc(1,l)=sc(1,l)+cij(i)*real(sca(l,i))
        enddo
      enddo
c computation of the phase function and the asymetry coefficient
c of the mixture of particles
      do l=1,10
        asy_n=0.
        asy_d=0.
        do k=1,nbmu
          ph(l,k)=0.
          do i=1,icp
           ph(l,k)=ph(l,k)+real(cij(i)*p1(l,i,k)/np(i)/1.D+3)
          enddo
          ph(l,k)=ph(l,k)/sc(1,l)
	  asy_n=asy_n+cgaus(k)*ph(l,k)*pdgs(k)/10.
	  asy_d=asy_d+ph(l,k)*pdgs(k)/10.
        enddo
	asy(1,l)=asy_n/asy_d
      enddo

      return
      END                                                                       
C***************************************************************************
C Using the Mie's theory, this subroutine compute the scattering and 
C extinction efficiency factors (usually written Qsca and Qext) and it also 
C compute the scattering intensity efficiency
      subroutine EXSCPHASE(X,nr,ni,Qext,Qsca,p11)
      parameter (nser=10000)
      double precision Ren,Imn,X,Up,XnumRDnY,XnumIDnY
      double precision XdenDnY,coxj,Qsca,Qext,xJonH,XdenGNX
      double precision Xnum1An,Xnum2An,XdenAn,Xden1An,Xden2An,RAnb,IAnb
      double precision Xnum1Bn,Xnum2Bn,XdenBn,Xden1Bn,Xden2Bn,RBnb,IBnb
      double precision xmud,xpond,RS1,RS2,IS1,IS2,co_n,test
      double precision xj(0:nser),xy(-1:nser),Rn(0:nser)
      double precision IDnY(0:nser),RDnX(0:nser),RDnY(0:nser)
      double precision IGnX(0:nser),RGnX(0:nser)
      double precision RAn(0:nser),IAn(0:nser),RBn(0:nser),IBn(0:nser)
      double precision TAUn(0:nser),PIn(0:nser),p11(83)
      real nr,ni,cgaus,phasel,pdgs
      integer N,Np,mu,mub,mu1,mu2,k,nbmu,j

      common /sixs_sos/ phasel(10,83),cgaus(83),pdgs(83)

      nbmu=83

      Ren=nr/(nr*nr+ni*ni)
      Imn=ni/(nr*nr+ni*ni)

c ---Identification of the greater order of computation (=mu)
c    as defined by F.J. Corbato, J. Assoc. Computing Machinery, 1959,
c    6, 366-375
      N=int(0.5D+00*(-1.D+00+dsqrt(1.D+00+4.D+00*X*X)))+1
      if (N.eq.1)N=2

      mu2=1000000
      Np=N
      Up=2.D+00*X/(2.D+00*Np+1.D+00)
      mu1=int(Np+30.*(0.10+0.35*Up*(2-Up*Up)/2./(1-Up)))
      Np=int(X-0.5D+00+dsqrt(30.*0.35*X))
      if (Np.gt.N)then
       Up=2.D+00*X/(2.D+00*Np+1.D+00)
       mu2=int(Np+30.*(0.10+0.35*Up*(2-Up*Up)/2./(1-Up)))
      endif
      mu=min0(mu1,mu2)

c --- Identification of the transition line. Below this line the Bessel 
c     function j behaves as oscillating functions. Above the behavior 
c     becomes monotonic. We start at a order greater than this transition 
c     line (order max=mu) because a downward recursion is called for.
      Rn(mu)=0.D+00
      k=mu+1
 149  continue
      k=k-1
      xj(k)=0.D+00
      Rn(k-1)=X/(2.D+00*k+1.D+00-X*Rn(k))
      if (k.eq.2)then
	  mub=mu
	  xj(mub+1)=0.D+00
	  xj(mub)=1.D+00
	  goto 150
      endif
      if (Rn(k-1).gt.1.D+00)then
	  mub=k-1
	  xj(mub+1)=Rn(mub)
	  xj(mub)=1.D+00
	  goto 150
      endif
      goto 149
 150  continue

      do k=mub,1,-1
	xj(k-1)=(2.D+00*k+1.D+00)*xj(k)/X-xj(k+1)
      enddo
      coxj=(xj(0)-X*xj(1))*dcos(X)+X*xj(0)*sin(X)

c --- Computation Dn(alpha) and Dn(alpha*m) (cf MIE's theory) 
c     downward recursion    - real and imaginary parts
      RDnY(mu)=0.D+00
      IDnY(mu)=0.D+00
      RDnX(mu)=0.D+00
      do k=mu,1,-1
	 RDnX(k-1)=k/X-1.D+00/(RDnX(k)+k/X)
	 XnumRDnY=RDnY(k)+Ren*k/X
	 XnumIDnY=IDnY(k)+Imn*k/X
	 XdenDnY=XnumRDnY*XnumRDnY+XnumIDnY*XnumIDnY
	 RDnY(k-1)=k*Ren/X-XnumRDnY/XdenDnY
	 IDnY(k-1)=k*Imn/X+XnumIDnY/XdenDnY

      enddo

c --- Initialization of the upward recursions
      xy(-1)=dsin(x)/x
      xy(0)=-dcos(x)/x
      RGnX(0)=0.D+00
      IGnX(0)=-1.D+00
      Qsca=0.D+00
      Qext=0.D+00
      do k=1,mu
	 if (k.le.mub)then
	   xj(k)=xj(k)/coxj
	 else
	   xj(k)=Rn(k-1)*xj(k-1)
	 endif

c --- Computation of bessel's function y(alpha)
	 xy(k)=(2.D+00*k-1.D+00)*xy(k-1)/X-xy(k-2)
	 xJonH=xj(k)/(xj(k)*xj(k)+xy(k)*xy(k))

c --- Computation of Gn(alpha), Real and Imaginary part
         XdenGNX=(RGnX(k-1)-k/X)**2.D+00+IGnX(k-1)*IGnX(k-1)
	 RGnX(k)=(k/X-RGnX(k-1))/XdenGNX-k/X
	 IGnX(k)=IGnX(k-1)/XdenGNX

c --- Computation of An(alpha) and Bn(alpha), Real and Imaginary part
	 Xnum1An=RDnY(k)-nr*RDnX(k)
	 Xnum2An=IDnY(k)+ni*RDnX(k)
	 Xden1An=RDnY(k)-nr*RGnX(k)-ni*IGnX(k)
	 Xden2An=IDnY(k)+ni*RGnX(k)-nr*IGnX(k)
	 XdenAn=Xden1An*Xden1An+Xden2An*Xden2An
	 RAnb=(Xnum1An*Xden1An+Xnum2An*Xden2An)/XdenAn
	 IAnb=(-Xnum1An*Xden2An+Xnum2An*Xden1An)/XdenAn
	 RAn(k)=xJonH*(xj(k)*RAnb-xy(k)*IAnb)
	 IAn(k)=xJonH*(xy(k)*RAnb+xj(k)*IAnb)

	 Xnum1Bn=nr*RDnY(k)+ni*IDnY(k)-RDnX(k)
	 Xnum2Bn=nr*IDnY(k)-ni*RDnY(k)
	 Xden1Bn=nr*RDnY(k)+ni*IDnY(k)-RGnX(k)
	 Xden2Bn=nr*IDnY(k)-ni*RDnY(k)-IGnX(k)
	 XdenBn=Xden1Bn*Xden1Bn+Xden2Bn*Xden2Bn
	 RBnb=(Xnum1Bn*Xden1Bn+Xnum2Bn*Xden2Bn)/XdenBn
	 IBnb=(-Xnum1Bn*Xden2Bn+Xnum2Bn*Xden1Bn)/XdenBn
	 RBn(k)=xJonH*(xj(k)*RBnb-xy(k)*IBnb)
	 IBn(k)=xJonH*(xy(k)*RBnb+xj(k)*IBnb)

c ---Criterion on the recursion formulas as defined by D. Deirmendjian 
c    et al., J. Opt. Soc. Am., 1961, 51, 6, 620-633
 	 test=(RAn(k)**2.+IAn(k)**2.+RBn(k)**2.+IBn(k)**2.)/k
 	 if (test.lt.1.0D-14)then
           mu=k
           goto 400
         endif
c --- Computation of the scattering and extinction efficiency factor
         xpond=2.D+00/X/X*(2.D+00*k+1)
         Qsca=Qsca+xpond*(RAn(k)**2.+IAn(k)**2.+RBn(k)**2.+IBn(k)**2.)
         Qext=Qext+xpond*(RAn(k)+RBn(k))

      enddo
 400  continue

c --- Computation of the amplitude functions S1 and S2 (cf MIE's theory)
c     defined by PIn, TAUn, An and Bn with PIn and TAUn related to the 
c     Legendre polynomials.
      do j=1,nbmu
	 xmud=cgaus(j)
	 RS1=0.D+00
	 RS2=0.D+00
	 IS1=0.D+00
	 IS2=0.D+00
	 PIn(0)=0.D+00
	 PIn(1)=1.D+00
	 TAUn(1)=xmud
	 do k=1,mu
          co_n=(2.D+00*k+1.D+00)/k/(k+1.D+00)
	  RS1=RS1+co_n*(RAn(k)*PIn(k)+RBn(k)*TAUn(k))
	  RS2=RS2+co_n*(RAn(k)*TAUn(k)+RBn(k)*PIn(k))
	  IS1=IS1+co_n*(IAn(k)*PIn(k)+IBn(k)*TAUn(k))
	  IS2=IS2+co_n*(IAn(k)*TAUn(k)+IBn(k)*PIn(k))
          PIn(k+1)=((2.D+00*k+1)*xmud*PIn(k)-(k+1.D+00)*PIn(k-1))/k
          TAUn(k+1)=(k+1.D+00)*xmud*PIn(k+1)-(k+2.D+00)*PIn(k)
         enddo
C --- Computation of the scattering intensity efficiency
         p11(j)=2.D+00*(RS1*RS1+IS1*IS1+RS2*RS2+IS2*IS2)/X/X
      enddo
      return
      end

      block data aeroso_data
      common /sixs_sos/phasel(10,83),cgaus(83),pdgs(83)
      real phasel,cgaus,pdgs
      data cgaus/
     a-1.0000,-0.9996,-0.9976,-0.9942,-0.9893,-0.9828,-0.9749,-0.9655,
     a-0.9546,-0.9422,-0.9285,-0.9133,-0.8967,-0.8787,-0.8594,-0.8388,
     a-0.8170,-0.7938,-0.7695,-0.7440,-0.7174,-0.6896,-0.6609,-0.6311,
     a-0.6003,-0.5687,-0.5361,-0.5028,-0.4687,-0.4339,-0.3984,-0.3623,
     a-0.3257,-0.2885,-0.2510,-0.2130,-0.1747,-0.1362,-0.0974,-0.0585,
     a-0.0195, 0.0000, 0.0195, 0.0585, 0.0974, 0.1362, 0.1747, 0.2130,
     a 0.2510, 0.2885, 0.3257, 0.3623, 0.3984, 0.4339, 0.4687, 0.5028,
     a 0.5361, 0.5687, 0.6003, 0.6311, 0.6609, 0.6896, 0.7174, 0.7440,
     a 0.7695, 0.7938, 0.8170, 0.8388, 0.8594, 0.8787, 0.8967, 0.9133,
     a 0.9285, 0.9422, 0.9546, 0.9655, 0.9749, 0.9828, 0.9893, 0.9942,
     a 0.9976, 0.9996, 1.0000/
      data pdgs/
     a 0.0000, 0.0114, 0.0266, 0.0418, 0.0569, 0.0719, 0.0868, 0.1016,
     a 0.1162, 0.1307, 0.1449, 0.1590, 0.1727, 0.1863, 0.1995, 0.2124,
     a 0.2251, 0.2373, 0.2492, 0.2606, 0.2719, 0.2826, 0.2929, 0.3027,
     a 0.3121, 0.3210, 0.3294, 0.3373, 0.3447, 0.3516, 0.3579, 0.3637,
     a 0.3690, 0.3737, 0.3778, 0.3813, 0.3842, 0.3866, 0.3884, 0.3896,
     a 0.3902, 0.0000, 0.3902, 0.3896, 0.3884, 0.3866, 0.3842, 0.3813,
     a 0.3778, 0.3737, 0.3690, 0.3637, 0.3579, 0.3516, 0.3447, 0.3373,
     a 0.3294, 0.3210, 0.3121, 0.3027, 0.2929, 0.2826, 0.2719, 0.2606,
     a 0.2492, 0.2373, 0.2251, 0.2124, 0.1995, 0.1863, 0.1727, 0.1590,
     a 0.1449, 0.1307, 0.1162, 0.1016, 0.0868, 0.0719, 0.0569, 0.0418,
     a 0.0266, 0.0114, 0.0000/
      end
      subroutine   ocea
      common /sixs_aerbas/ ph(10,83)
      real phr(10,83),ph
      integer i,j
c
c    model: oceanic
c
            DATA ((PHR(I,J),J=1,83),I=01,01) /
     *0.7855E+00,0.6283E+00,0.5465E+00,0.4693E+00,0.4153E+00,0.3917E+00,
     *0.3657E+00,0.3378E+00,0.3161E+00,0.3025E+00,0.2972E+00,0.2990E+00,
     *0.3055E+00,0.3118E+00,0.3059E+00,0.2715E+00,0.2118E+00,0.1585E+00,
     *0.1230E+00,0.9913E-01,0.8327E-01,0.7292E-01,0.6585E-01,0.6171E-01,
     *0.5883E-01,0.5780E-01,0.5791E-01,0.5893E-01,0.6144E-01,0.6406E-01,
     *0.6717E-01,0.6966E-01,0.7130E-01,0.7291E-01,0.7434E-01,0.7626E-01,
     *0.7847E-01,0.8190E-01,0.8583E-01,0.9044E-01,0.9709E-01,0.1006E+00,
     *0.1045E+00,0.1128E+00,0.1239E+00,0.1360E+00,0.1497E+00,0.1667E+00,
     *0.1856E+00,0.2070E+00,0.2323E+00,0.2615E+00,0.2948E+00,0.3326E+00,
     *0.3772E+00,0.4263E+00,0.4840E+00,0.5492E+00,0.6242E+00,0.7103E+00,
     *0.8075E+00,0.9192E+00,0.1046E+01,0.1190E+01,0.1354E+01,0.1541E+01,
     *0.1756E+01,0.2002E+01,0.2277E+01,0.2603E+01,0.2976E+01,0.3416E+01,
     *0.3931E+01,0.4563E+01,0.5372E+01,0.6490E+01,0.8191E+01,0.1111E+02,
     *0.1692E+02,0.3097E+02,0.7524E+02,0.2992E+03,0.1697E+04/
            DATA ((PHR(I,J),J=1,83),I=02,02) /
     *0.7129E+00,0.5739E+00,0.5059E+00,0.4429E+00,0.4035E+00,0.3898E+00,
     *0.3678E+00,0.3416E+00,0.3195E+00,0.3042E+00,0.2975E+00,0.2961E+00,
     *0.2987E+00,0.2994E+00,0.2909E+00,0.2614E+00,0.2134E+00,0.1670E+00,
     *0.1336E+00,0.1100E+00,0.9363E-01,0.8252E-01,0.7480E-01,0.6967E-01,
     *0.6621E-01,0.6499E-01,0.6438E-01,0.6506E-01,0.6656E-01,0.6880E-01,
     *0.7108E-01,0.7332E-01,0.7497E-01,0.7681E-01,0.7860E-01,0.8093E-01,
     *0.8357E-01,0.8723E-01,0.9184E-01,0.9665E-01,0.1036E+00,0.1075E+00,
     *0.1112E+00,0.1200E+00,0.1316E+00,0.1436E+00,0.1580E+00,0.1748E+00,
     *0.1937E+00,0.2154E+00,0.2413E+00,0.2704E+00,0.3031E+00,0.3421E+00,
     *0.3856E+00,0.4356E+00,0.4928E+00,0.5586E+00,0.6333E+00,0.7196E+00,
     *0.8188E+00,0.9313E+00,0.1060E+01,0.1208E+01,0.1375E+01,0.1568E+01,
     *0.1791E+01,0.2047E+01,0.2340E+01,0.2679E+01,0.3075E+01,0.3547E+01,
     *0.4107E+01,0.4805E+01,0.5714E+01,0.6981E+01,0.8889E+01,0.1212E+02,
     *0.1839E+02,0.3283E+02,0.7515E+02,0.2626E+03,0.1134E+04/
            DATA ((PHR(I,J),J=1,83),I=03,03) /
     *0.6966E+00,0.5607E+00,0.4902E+00,0.4336E+00,0.3978E+00,0.3866E+00,
     *0.3674E+00,0.3412E+00,0.3187E+00,0.3039E+00,0.2960E+00,0.2945E+00,
     *0.2960E+00,0.2961E+00,0.2874E+00,0.2591E+00,0.2133E+00,0.1692E+00,
     *0.1362E+00,0.1129E+00,0.9630E-01,0.8484E-01,0.7707E-01,0.7190E-01,
     *0.6854E-01,0.6653E-01,0.6597E-01,0.6668E-01,0.6812E-01,0.7009E-01,
     *0.7216E-01,0.7425E-01,0.7580E-01,0.7758E-01,0.7959E-01,0.8174E-01,
     *0.8490E-01,0.8852E-01,0.9294E-01,0.9864E-01,0.1048E+00,0.1084E+00,
     *0.1128E+00,0.1220E+00,0.1325E+00,0.1453E+00,0.1596E+00,0.1762E+00,
     *0.1959E+00,0.2177E+00,0.2428E+00,0.2725E+00,0.3055E+00,0.3440E+00,
     *0.3882E+00,0.4382E+00,0.4953E+00,0.5613E+00,0.6365E+00,0.7225E+00,
     *0.8218E+00,0.9344E+00,0.1065E+01,0.1212E+01,0.1381E+01,0.1577E+01,
     *0.1801E+01,0.2059E+01,0.2360E+01,0.2701E+01,0.3107E+01,0.3586E+01,
     *0.4166E+01,0.4885E+01,0.5821E+01,0.7115E+01,0.9088E+01,0.1241E+02,
     *0.1877E+02,0.3323E+02,0.7480E+02,0.2523E+03,0.1018E+04/
            DATA ((PHR(I,J),J=1,83),I=04,04) /
     *0.6774E+00,0.5476E+00,0.4775E+00,0.4252E+00,0.3937E+00,0.3855E+00,
     *0.3684E+00,0.3432E+00,0.3209E+00,0.3059E+00,0.2974E+00,0.2950E+00,
     *0.2951E+00,0.2935E+00,0.2832E+00,0.2550E+00,0.2114E+00,0.1697E+00,
     *0.1380E+00,0.1153E+00,0.9882E-01,0.8737E-01,0.7952E-01,0.7423E-01,
     *0.7074E-01,0.6859E-01,0.6788E-01,0.6842E-01,0.6969E-01,0.7150E-01,
     *0.7349E-01,0.7557E-01,0.7720E-01,0.7911E-01,0.8125E-01,0.8356E-01,
     *0.8685E-01,0.9062E-01,0.9516E-01,0.1010E+00,0.1073E+00,0.1109E+00,
     *0.1154E+00,0.1247E+00,0.1352E+00,0.1482E+00,0.1626E+00,0.1793E+00,
     *0.1991E+00,0.2210E+00,0.2462E+00,0.2760E+00,0.3091E+00,0.3477E+00,
     *0.3920E+00,0.4422E+00,0.4994E+00,0.5656E+00,0.6410E+00,0.7275E+00,
     *0.8272E+00,0.9405E+00,0.1071E+01,0.1220E+01,0.1391E+01,0.1588E+01,
     *0.1815E+01,0.2077E+01,0.2382E+01,0.2731E+01,0.3145E+01,0.3636E+01,
     *0.4233E+01,0.4974E+01,0.5942E+01,0.7282E+01,0.9319E+01,0.1273E+02,
     *0.1919E+02,0.3364E+02,0.7414E+02,0.2397E+03,0.8914E+03/
            DATA ((PHR(I,J),J=1,83),I=05,05) /
     *0.6153E+00,0.5058E+00,0.4382E+00,0.3950E+00,0.3738E+00,0.3731E+00,
     *0.3585E+00,0.3354E+00,0.3139E+00,0.2983E+00,0.2892E+00,0.2849E+00,
     *0.2832E+00,0.2800E+00,0.2703E+00,0.2469E+00,0.2112E+00,0.1741E+00,
     *0.1442E+00,0.1219E+00,0.1054E+00,0.9356E-01,0.8531E-01,0.7966E-01,
     *0.7561E-01,0.7323E-01,0.7198E-01,0.7214E-01,0.7291E-01,0.7415E-01,
     *0.7601E-01,0.7747E-01,0.7901E-01,0.8091E-01,0.8293E-01,0.8564E-01,
     *0.8906E-01,0.9289E-01,0.9788E-01,0.1033E+00,0.1102E+00,0.1141E+00,
     *0.1181E+00,0.1275E+00,0.1385E+00,0.1511E+00,0.1660E+00,0.1823E+00,
     *0.2018E+00,0.2241E+00,0.2491E+00,0.2784E+00,0.3123E+00,0.3503E+00,
     *0.3942E+00,0.4451E+00,0.5020E+00,0.5684E+00,0.6448E+00,0.7319E+00,
     *0.8325E+00,0.9481E+00,0.1081E+01,0.1234E+01,0.1409E+01,0.1612E+01,
     *0.1846E+01,0.2118E+01,0.2440E+01,0.2809E+01,0.3249E+01,0.3773E+01,
     *0.4413E+01,0.5211E+01,0.6259E+01,0.7710E+01,0.9888E+01,0.1347E+02,
     *0.2009E+02,0.3435E+02,0.7217E+02,0.2130E+03,0.6728E+03/
            DATA ((PHR(I,J),J=1,83),I=06,06) /
     *0.5916E+00,0.4877E+00,0.4171E+00,0.3786E+00,0.3632E+00,0.3654E+00,
     *0.3546E+00,0.3335E+00,0.3124E+00,0.2967E+00,0.2869E+00,0.2822E+00,
     *0.2792E+00,0.2744E+00,0.2635E+00,0.2413E+00,0.2085E+00,0.1740E+00,
     *0.1459E+00,0.1244E+00,0.1084E+00,0.9682E-01,0.8822E-01,0.8243E-01,
     *0.7835E-01,0.7606E-01,0.7463E-01,0.7441E-01,0.7473E-01,0.7609E-01,
     *0.7739E-01,0.7905E-01,0.8078E-01,0.8256E-01,0.8474E-01,0.8745E-01,
     *0.9082E-01,0.9490E-01,0.9996E-01,0.1057E+00,0.1127E+00,0.1166E+00,
     *0.1207E+00,0.1301E+00,0.1412E+00,0.1539E+00,0.1686E+00,0.1858E+00,
     *0.2048E+00,0.2270E+00,0.2528E+00,0.2818E+00,0.3154E+00,0.3545E+00,
     *0.3980E+00,0.4487E+00,0.5067E+00,0.5728E+00,0.6491E+00,0.7374E+00,
     *0.8386E+00,0.9547E+00,0.1090E+01,0.1244E+01,0.1423E+01,0.1630E+01,
     *0.1870E+01,0.2149E+01,0.2477E+01,0.2862E+01,0.3316E+01,0.3862E+01,
     *0.4527E+01,0.5365E+01,0.6458E+01,0.7974E+01,0.1023E+02,0.1390E+02,
     *0.2058E+02,0.3459E+02,0.7042E+02,0.1961E+03,0.5608E+03/
            DATA ((PHR(I,J),J=1,83),I=07,07) /
     *0.5164E+00,0.4330E+00,0.3650E+00,0.3341E+00,0.3313E+00,0.3413E+00,
     *0.3356E+00,0.3182E+00,0.2998E+00,0.2844E+00,0.2744E+00,0.2677E+00,
     *0.2626E+00,0.2560E+00,0.2453E+00,0.2267E+00,0.2009E+00,0.1730E+00,
     *0.1485E+00,0.1291E+00,0.1141E+00,0.1028E+00,0.9425E-01,0.8828E-01,
     *0.8375E-01,0.8105E-01,0.7927E-01,0.7843E-01,0.7860E-01,0.7925E-01,
     *0.8010E-01,0.8165E-01,0.8331E-01,0.8499E-01,0.8754E-01,0.9034E-01,
     *0.9390E-01,0.9825E-01,0.1034E+00,0.1093E+00,0.1164E+00,0.1203E+00,
     *0.1246E+00,0.1342E+00,0.1452E+00,0.1582E+00,0.1728E+00,0.1896E+00,
     *0.2094E+00,0.2310E+00,0.2569E+00,0.2863E+00,0.3195E+00,0.3587E+00,
     *0.4030E+00,0.4534E+00,0.5122E+00,0.5794E+00,0.6565E+00,0.7463E+00,
     *0.8505E+00,0.9697E+00,0.1109E+01,0.1270E+01,0.1457E+01,0.1674E+01,
     *0.1929E+01,0.2226E+01,0.2578E+01,0.2997E+01,0.3495E+01,0.4096E+01,
     *0.4831E+01,0.5758E+01,0.6967E+01,0.8629E+01,0.1105E+02,0.1487E+02,
     *0.2152E+02,0.3465E+02,0.6548E+02,0.1595E+03,0.3700E+03/
            DATA ((PHR(I,J),J=1,83),I=08,08) /
     *0.3257E+00,0.2888E+00,0.2378E+00,0.2215E+00,0.2345E+00,0.2532E+00,
     *0.2578E+00,0.2504E+00,0.2390E+00,0.2282E+00,0.2194E+00,0.2123E+00,
     *0.2059E+00,0.1991E+00,0.1906E+00,0.1797E+00,0.1665E+00,0.1520E+00,
     *0.1379E+00,0.1254E+00,0.1147E+00,0.1061E+00,0.9917E-01,0.9373E-01,
     *0.8960E-01,0.8656E-01,0.8438E-01,0.8306E-01,0.8243E-01,0.8240E-01,
     *0.8294E-01,0.8394E-01,0.8543E-01,0.8740E-01,0.8990E-01,0.9302E-01,
     *0.9681E-01,0.1013E+00,0.1067E+00,0.1129E+00,0.1200E+00,0.1240E+00,
     *0.1283E+00,0.1379E+00,0.1490E+00,0.1618E+00,0.1764E+00,0.1932E+00,
     *0.2124E+00,0.2345E+00,0.2599E+00,0.2892E+00,0.3231E+00,0.3622E+00,
     *0.4072E+00,0.4593E+00,0.5195E+00,0.5895E+00,0.6711E+00,0.7664E+00,
     *0.8781E+00,0.1009E+01,0.1163E+01,0.1343E+01,0.1556E+01,0.1808E+01,
     *0.2107E+01,0.2464E+01,0.2891E+01,0.3405E+01,0.4025E+01,0.4779E+01,
     *0.5707E+01,0.6863E+01,0.8338E+01,0.1027E+02,0.1291E+02,0.1670E+02,
     *0.2248E+02,0.3211E+02,0.5001E+02,0.8772E+02,0.1334E+03/
            DATA ((PHR(I,J),J=1,83),I=09,09) /
     *0.2139E+00,0.1949E+00,0.1618E+00,0.1541E+00,0.1685E+00,0.1828E+00,
     *0.1856E+00,0.1800E+00,0.1718E+00,0.1642E+00,0.1581E+00,0.1534E+00,
     *0.1495E+00,0.1460E+00,0.1421E+00,0.1375E+00,0.1318E+00,0.1252E+00,
     *0.1178E+00,0.1105E+00,0.1036E+00,0.9754E-01,0.9237E-01,0.8811E-01,
     *0.8468E-01,0.8198E-01,0.7994E-01,0.7852E-01,0.7768E-01,0.7741E-01,
     *0.7767E-01,0.7843E-01,0.7969E-01,0.8144E-01,0.8373E-01,0.8662E-01,
     *0.9014E-01,0.9438E-01,0.9939E-01,0.1052E+00,0.1120E+00,0.1158E+00,
     *0.1198E+00,0.1289E+00,0.1394E+00,0.1514E+00,0.1653E+00,0.1813E+00,
     *0.1997E+00,0.2208E+00,0.2453E+00,0.2736E+00,0.3064E+00,0.3444E+00,
     *0.3886E+00,0.4400E+00,0.5000E+00,0.5703E+00,0.6528E+00,0.7502E+00,
     *0.8652E+00,0.1001E+01,0.1163E+01,0.1355E+01,0.1584E+01,0.1859E+01,
     *0.2188E+01,0.2586E+01,0.3067E+01,0.3649E+01,0.4358E+01,0.5222E+01,
     *0.6282E+01,0.7594E+01,0.9235E+01,0.1132E+02,0.1404E+02,0.1768E+02,
     *0.2278E+02,0.3033E+02,0.4233E+02,0.6237E+02,0.7953E+02/
            DATA ((PHR(I,J),J=1,83),I=10,10) /
     *0.2110E+00,0.2025E+00,0.1832E+00,0.1730E+00,0.1773E+00,0.1912E+00,
     *0.2055E+00,0.2138E+00,0.2152E+00,0.2113E+00,0.2040E+00,0.1946E+00,
     *0.1842E+00,0.1734E+00,0.1627E+00,0.1524E+00,0.1429E+00,0.1344E+00,
     *0.1268E+00,0.1203E+00,0.1149E+00,0.1104E+00,0.1068E+00,0.1040E+00,
     *0.1019E+00,0.1006E+00,0.9982E-01,0.9972E-01,0.1003E+00,0.1014E+00,
     *0.1031E+00,0.1054E+00,0.1084E+00,0.1119E+00,0.1162E+00,0.1212E+00,
     *0.1271E+00,0.1338E+00,0.1415E+00,0.1503E+00,0.1603E+00,0.1658E+00,
     *0.1717E+00,0.1847E+00,0.1995E+00,0.2163E+00,0.2354E+00,0.2571E+00,
     *0.2818E+00,0.3100E+00,0.3422E+00,0.3792E+00,0.4216E+00,0.4702E+00,
     *0.5261E+00,0.5903E+00,0.6644E+00,0.7500E+00,0.8493E+00,0.9645E+00,
     *0.1098E+01,0.1254E+01,0.1436E+01,0.1649E+01,0.1897E+01,0.2189E+01,
     *0.2531E+01,0.2934E+01,0.3408E+01,0.3968E+01,0.4630E+01,0.5415E+01,
     *0.6348E+01,0.7463E+01,0.8805E+01,0.1044E+02,0.1244E+02,0.1495E+02,
     *0.1816E+02,0.2237E+02,0.2799E+02,0.3517E+02,0.3934E+02/
      do 1 i=1,10
      do 1 j=1,83
      ph(i,j)=phr(i,j)
    1 continue
      return
      end
      subroutine oda550 (iaer,v,
     a                   taer55)
 
      double precision bnz,bnz1
      common /sixs_atm/ z(34),p(34),t(34),wh(34),wo(34)
      common /sixs_del/ delta,sigma
      real an5(34),an23(34)
      Real v,taer55,z,p,t,wh
      Real wo,delta,sigma,dz,bn5,bn51,bn23,bn231,az
      Real az1,bz,bz1,ev
      Integer iaer,k
c    aerosol optical depth at wl=550nm
c     vertical repartition of aerosol density for v=23km
c                     ( in nbr of part/cm3 )
 
      data an23 /2.828e+03,1.244e+03,5.371e+02,2.256e+02,1.192e+02
     a,8.987e+01,6.337e+01,5.890e+01,6.069e+01,5.818e+01,5.675e+01
     a,5.317e+01,5.585e+01,5.156e+01,5.048e+01,4.744e+01,4.511e+01
     a,4.458e+01,4.314e+01,3.634e+01,2.667e+01,1.933e+01,1.455e+01
     a,1.113e+01,8.826e+00,7.429e+00,2.238e+00,5.890e-01,1.550e-01
     a,4.082e-02,1.078e-02,5.550e-05,1.969e-08,0.000e+00/
 
c     vertical repartition of aerosol density for v=5km
c                     ( in nbr of part/cm3 )
 
      data  an5 /1.378e+04,5.030e+03,1.844e+03,6.731e+02,2.453e+02
     a,8.987e+01,6.337e+01,5.890e+01,6.069e+01,5.818e+01,5.675e+01
     a,5.317e+01,5.585e+01,5.156e+01,5.048e+01,4.744e+01,4.511e+01
     a,4.458e+01,4.314e+01,3.634e+01,2.667e+01,1.933e+01,1.455e+01
     a,1.113e+01,8.826e+00,7.429e+00,2.238e+00,5.890e-01,1.550e-01
     a,4.082e-02,1.078e-02,5.550e-05,1.969e-08,0.000e+00/
 
 
      taer55=0.
 
      if(abs(v).le.0.) return
      if(iaer.eq.0) return
 
      do 1 k=1,32
      dz=z(k+1)-z(k)
      bn5=an5(k)
      bn51=an5(k+1)
      bn23=an23(k)
      bn231=an23(k+1)
      az=(115./18.)*(bn5-bn23)
      az1=(115./18.)*(bn51-bn231)
      bz=(5.*bn5/18.)-(23.*bn23/18.)
      bz1=(5.*bn51/18.)-(23.*bn231/18.)
      bnz=az/v-bz
      bnz1=az1/v-bz1
      ev=dz*exp((dlog(bnz)+dlog(bnz1))*.5)
      taer55=taer55+ev*sigma*1.0e-03
    1 continue
      return
      end
      subroutine odrayl ( wl,
     a                   tray)
      double precision a1,a2,a3,a4,awl,an,a
      real wl,tray,z,p,t,wh,wo,delta,sigma,pi,ak,dppt,sr
      integer k
c     molecular optical depth
 
      common /sixs_atm/ z(34),p(34),t(34),wh(34),wo(34)
      common /sixs_del/ delta,sigma
      real ns
      data pi /3.1415926/
      ak=1/wl
      awl=wl
c     air refraction index edlen 1966 / metrologia,2,71-80  putting pw=0
      a1=130.-ak*ak
      a2=38.9-ak*ak
      a3=2406030./a1
      a4=15997./a2
      an=(8342.13+a3+a4)*1.0e-08
      an=an+1.d+00
      a=(24.*pi**3)*((an*an-1.)**2)*(6.+3.*delta)/(6.-7.*delta)
     s        /((an*an+2.)**2)
      tray=0.
      do k=1,33
      ns=2.54743e+19
      dppt=(288.15/1013.25)*(p(k)/t(k)+p(k+1)/t(k+1))/2.
      sr=(a*dppt/(awl**4)/ns*1.e+16)*1.e+05
      tray=tray+(z(k+1)-z(k))*sr
      enddo
      return
      end
      subroutine os (tamoy,trmoy,pizmoy,tamoyp,trmoyp,palt,
     s               phirad,nt,mu,np,rm,gb,rp,
     s                     xl)
c  dimension for gauss integration
      integer mu,np
      real rm(-mu:mu),gb(-mu:mu),rp(np)
c  dimension for os computation
      real  xl(-mu:mu,np)
c array for sos computation
      real xpl(-25:25),psl(-1:80,-25:25),bp(0:25,-25:25),
     s xdel(0:30),ydel(0:30),ch(0:30),h(0:30)
      real i1(0:30,-25:25),i2(0:30,-25:25),i3(-25:25),
     s   i4(-25:25),in(-25:25),inm1(-25:25),inm2(-25:25)
      real altc(0:30)
      Real tamoy,trmoy,pizmoy
      Real tamoyp,trmoyp,palt,phirad
      Real delta,sigma,pha,betal,hr,ta,tr,trp
      Real tap,piz,accu,accu2,ha,xmus,zx,yy,dd,ppp2,ppp1,ca,cr,ratio
      Real taup,th,xt1,xt2,pi,phi,aaaa,ron
      Real beta0,beta2,roavion0,roavion1,roavion2,roavion,spl,sa1
      Real sa2,c,zi1,f,d,xpk,y
      Real a1,d1,g1,y1,delta0s
      integer snt
      integer nt,iwr,iplane,mum1,ntp,j,it,itp,i,l,m,iborm
      integer is,isp,ig,k,jj,index
      logical ier
      common/sixs_del/delta,sigma
      common /sixs_trunc/pha(83),betal(0:80)
      common/sixs_ier/iwr,ier
      double precision xx,xdb,bpjk,bpjmk,z,xi1,xi2,x,xpj,ypk,a,b,ii1,ii2
c the optical thickness above plane are recomputed to give o.t above pla
c     write(6,*) 'tamoy,trmoy,tamoyp,trmoyp,palt,pizmoy'
c     write(6,*) tamoy,trmoy,tamoyp,trmoyp,palt,pizmoy
c     write(6,*) 'betal 0:80'
c     do i=0,80
c     write(6,*) i,betal(i)
c     enddo
c     write(6,*) 'phase function 83 terms'
c     do i=1,83
c     write(6,*) pha(i)
c     enddo
      snt=nt
      hr=8.0
      ta=tamoy
      tr=trmoy
      trp=trmoy-trmoyp
      tap=tamoy-tamoyp
      piz=pizmoy
c     print *, 'ta,tr,piz,tap,trp,palt,nt'
c     print *,ta,tr,piz,tap,trp,palt,nt
      iplane=0
      accu=1.e-20
      accu2=1.e-3
      mum1=mu-1
c if plane observations recompute scale height for aerosol knowing:
c the aerosol optical depth as measure from the plane 	= tamoyp
c the rayleigh scale   height = 			= hr (8km)
c the rayleigh optical depth  at plane level 		= trmoyp
c the altitude of the plane 				= palt
c the rayleigh optical depth for total atmos		= trmoy
c the aerosol  optical depth for total atmos		= tamoy
c if not plane observations then ha is equal to 2.0km
c ntp local variable: if ntp=nt     no plane observation selected
c                        ntp=nt-1   plane observation selected
c     it's a mixing rayleigh+aerosol
      if(palt.le.900..and.palt.gt.0.0) then
      if (tap.gt.1.e-03) then
         ha=-palt/log(tap/ta)
         else
         ha=2.
         endif
      ntp=nt-1
      else
      ha=2.0
      ntp=nt
      endif
c
      xmus=-rm(0)
c
c compute mixing rayleigh, aerosol
c case 1: pure rayleigh
c case 2: pure aerosol
c case 3: mixing rayleigh-aerosol
c
      if((ta.le.accu2).and.(tr.gt.ta)) then
      do j=0,ntp
      h(j)=j*tr/ntp
      ch(j)=exp(-h(j)/xmus)/2.
      ydel(j)=1.0
      xdel(j)=0.0
      if (j.eq.0) then
         altc(j)=300.
         else
         altc(j)=-log(h(j)/tr)*hr
         endif
      enddo
      endif
      if((tr.le.accu2).and.(ta.gt.tr)) then
      do j=0,ntp
      h(j)=j*ta/ntp
      ch(j)=exp(-h(j)/xmus)/2.
      ydel(j)=0.0
      xdel(j)=piz
      if (j.eq.0) then
         altc(j)=300.
         else
         altc(j)=-log(h(j)/ta)*ha
         endif
      enddo
      endif
c
      if(tr.gt.accu2.and.ta.gt.accu2)then
      ydel(0)=1.0
      xdel(0)=0.0
      h(0)=0.
      ch(0)=0.5
      altc(0)=300.
      zx=300.
      iplane=0
      do 14 it=0,ntp
      if (it.eq.0) then
         yy=0.
         dd=0.
         goto 111
      endif
      yy=h(it-1)
      dd=ydel(it-1)
 111  ppp2=300.0
      ppp1=0.0
      itp=it
      call discre(ta,ha,tr,hr,itp,ntp,yy,dd,ppp2,ppp1,
     s    zx)
      if(ier)return
      xx=-zx/ha
      if (xx.le.-20) then
         ca=0.
         else
         ca=ta*dexp(xx)
         endif
      xx=-zx/hr
      cr=tr*dexp(xx)
      h(it)=cr+ca
      altc(it)=zx
      ch(it)=exp(-h(it)/xmus)/2.
      cr=cr/hr
      ca=ca/ha
      ratio=cr/(cr+ca)
      xdel(it)=(1.e+00-ratio)*piz
      ydel(it)=ratio
c     print *,'discre ',it,cr,ca,xdel(it),ydel(it),zx
  14  continue
      endif
c update plane layer if necessary
      if (ntp.eq.(nt-1)) then
c compute position of the plane layer
         taup=tap+trp
         iplane=-1
         do i=0,ntp
         if (taup.ge.h(i)) iplane=i
         enddo
c update the layer from the end to the position to update if necessary
         th=0.005
         xt1=abs(h(iplane)-taup)
         xt2=abs(h(iplane+1)-taup)
         if ((xt1.gt.th).and.(xt2.gt.th)) then
         do i=nt,iplane+1,-1
            xdel(i)=xdel(i-1)
            ydel(i)=ydel(i-1)
            h(i)=h(i-1)
            altc(i)=altc(i-1)
            ch(i)=ch(i-1)
         enddo
         else
         nt=ntp
         if (xt2.lt.xt1) iplane=iplane+1
         endif
         h(iplane)=taup
         if ( tr.gt.accu2.and.ta.gt.accu2) then
         ca=ta*exp(-palt/ha)
         cr=tr*exp(-palt/hr)
         h(iplane)=ca+cr
         cr=cr/hr
         ca=ca/ha
         ratio=cr/(cr+ca)
         xdel(iplane)=(1.e+00-ratio)*piz
         ydel(iplane)=ratio
         altc(iplane)=palt
         ch(iplane)=exp(-h(iplane)/xmus)/2.
         endif
         if ( tr.gt.accu2.and.ta.le.accu2) then
         ydel(iplane)=1.
         xdel(iplane)=0.
         altc(iplane)=palt
         endif
         if ( tr.le.accu2.and.ta.gt.accu2) then
         ydel(iplane)=0.
         xdel(iplane)=1.*piz
         altc(iplane)=palt
         endif
      endif
c
c
c     print *,ha,hr,palt,ta,tr,tap,trp
c     do i=0,nt
c     print *,i,h(i),ch(i),xdel(i),ydel(i),altc(i)
c     enddo
c
      pi=acos(-1.)
      phi=phirad
      do 615 l=1,np
      do 615 m=-mu,mu
 615  xl(m,l)=0.
c
c     ************ incident angle mus *******
c
c
      aaaa=delta/(2-delta)
      ron=(1-aaaa)/(1+2*aaaa)
c     write(6,*) 'ron ',ron
c
c     rayleigh phase function
c
      beta0=1.
      beta2=0.5*ron
c
c     fourier decomposition
c
      do 17 j=-mu,mu
      i4(j)=0.
   17 continue
      iborm=80
      if( abs (xmus-1.000000) .lt.1.e-06)iborm=0
      do 24 is=0,iborm
c
c    primary scattering
c
      ig=1
      roavion0=0.
      roavion1=0.
      roavion2=0.
      roavion=0.
      do 16 j=-mu,mu
      i3(j)=0.
   16 continue
c
c     kernel computations
c
      isp=is
      call kernel(isp,mu,rm,xpl,psl,bp)
      if(is.gt.0)beta0=0.
      do 100 j=-mu,mu
      if(is-2)200,200,201
 200  spl=xpl(0)
      sa1=beta0+beta2*xpl(j)*spl
      sa2=bp(0,j)
      goto 202
 201  sa2=bp(0,j)
      sa1=0.
c
c     primary scattering source function at every level within the layer
c
 202  do 101 k=0,nt
      c=ch(k)
      a=ydel(k)
      b=xdel(k)
      i2(k,j)=c*(sa2*b+sa1*a)
  101 continue
  100 continue
c
c     vertical integration, primary upward radiation
c
 
      do 108 k=1,mu
      i1(nt,k)=0.
      zi1=i1(nt,k)
      yy=rm(k)
      do 108 i=nt-1,0,-1
      jj=i+1
      f=h(jj)-h(i)
      a=(i2(jj,k)-i2(i,k))/f
      b=i2(i,k)-a*h(i)
      c=exp(-f/yy)
      d=1.0e+00-c
      xx=h(i)-h(jj)*c
      zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5e+00
      i1(i,k)=zi1
  108 continue
c
c     vertical integration, primary downward radiation
c
      do 109 k=-mu,-1
      i1(0,k)=0.
      zi1=i1(0,k)
      yy=rm(k)
      do 109 i=1,nt
      jj=i-1
      f=h(i)-h(jj)
      c=exp(f/yy)
      d=1.0e+00-c
      a=(i2(i,k)-i2(jj,k))/f
      b=i2(i,k)-a*h(i)
      xx=h(i)-h(jj)*c
      zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5e+00
      i1(i,k)=zi1
  109 continue
c
c     inm2 is inialized with scattering computed at n-2
c     i3 is inialized with primary scattering
c
      do 20 k=-mu,mu
      if(k) 21,20,23
   21 index=nt
      go to 25
   23 index=0
   25 continue
      inm1(k)=i1(index,k)
      inm2(k)=i1(index,k)
      i3(k)=i1(index,k)
   20 continue
      roavion2=i1(iplane,mu)
      roavion=i1(iplane,mu)
c
c     loop on successive order
c
  503 ig=ig+1
c     write(6,*) 'ig ',ig
c
c     successive orders
c
c     multiple scattering source function at every level within the laye
c
c     if is < ou = 2 kernels are a mixing of aerosols and molecules kern
c     if is >2 aerosols kernels only
c
      if(is-2)210,210,211
  210 do455 k=1,mu
      xpk=xpl(k)
      ypk=xpl(-k)
      do 455 i=0,nt
      ii1=0.
      ii2=0.
      x=xdel(i)
      y=ydel(i)
      do477 j=1,mu
      xpj=xpl(j)
      z=gb(j)
      xi1=i1(i,j)
      xi2=i1(i,-j)
      bpjk=bp(j,k)*x+y*(beta0+beta2*xpj*xpk)
      bpjmk=bp(j,-k)*x+y*(beta0+beta2*xpj*ypk)
      xdb=z*(xi1*bpjk+xi2*bpjmk)
      ii2=ii2+xdb
      xdb=z*(xi1*bpjmk+xi2*bpjk)
      ii1=ii1+xdb
 477  continue
      if (ii2.lt.1.E-30) ii2=0.
      if (ii1.lt.1.E-30) ii1=0.
      i2(i,k)=ii2
      i2(i,-k)=ii1
 455  continue
      goto 213
 211  do45 k=1,mu
      do 45 i=0,nt
      ii1=0.
      ii2=0.
      x=xdel(i)
      do47 j=1,mu
      z=gb(j)
      xi1=i1(i,j)
      xi2=i1(i,-j)
      bpjk=bp(j,k)*x
      bpjmk=bp(j,-k)*x
      xdb=z*(xi1*bpjk+xi2*bpjmk)
      ii2=ii2+xdb
      xdb=z*(xi1*bpjmk+xi2*bpjk)
      ii1=ii1+xdb
   47 continue
      if (ii2.lt.1.E-30) ii2=0.
      if (ii1.lt.1.E-30) ii1=0.
      i2(i,k)=ii2
      i2(i,-k)=ii1
   45 continue
c
c     vertical integration, upward radiation
c
 213  do 48 k=1,mu
      i1(nt,k)=0.
      zi1=i1(nt,k)
      yy=rm(k)
      do 48 i=nt-1,0,-1
      jj=i+1
      f=h(jj)-h(i)
      a=(i2(jj,k)-i2(i,k))/f
      b=i2(i,k)-a*h(i)
      c=exp(-f/yy)
      d=1.e+00-c
      xx=h(i)-h(jj)*c
      zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5e+00
      if (abs(zi1).le.1.E-20) zi1=0.
      i1(i,k)=zi1
   48 continue
c
c     vertical integration, downward radiation
c
      do 50 k=-mu,-1
      i1(0,k)=0.
      zi1=i1(0,k)
      yy=rm(k)
      do 50 i=1,nt
      jj=i-1
      f=h(i)-h(jj)
      c=exp(f/yy)
      d=1.e+00-c
      a=(i2(i,k)-i2(jj,k))/f
      b=i2(i,k)-a*h(i)
      xx=h(i)-h(jj)*c
      zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5e+00
      if (abs(zi1).le.1.E-20) zi1=0.
      i1(i,k)=zi1
   50 continue
c
c     in is the nieme scattering order
c
      do 30 k=-mu,mu
      if(k) 31,30,33
   31 index=nt
      go to 34
   33 index=0
   34 continue
      in(k)=i1(index,k)
   30 continue
      roavion0=i1(iplane,mu)
c
c   convergence test (geometrical serie)
c
      if(ig.gt.2) then
      z=0.
      a1=roavion2
      d1=roavion1
      g1=roavion0
      if(a1.ge.accu.and.d1.ge.accu.and.roavion.ge.accu) then
      y=((g1/d1-d1/a1)/((1-g1/d1)**2)*(g1/roavion))
      y=abs(y)
      z=dmax1(dble(y),z)
      endif
      do 99 l=-mu,mu
      if (l.eq.0) goto 99
      a1=inm2(l)
      d1=inm1(l)
      g1=in(l)
      if(a1.le.accu) go to 99
      if(d1.le.accu) go to 99
      if(i3(l).le.accu) go to 99
      y=((g1/d1-d1/a1)/((1-g1/d1)**2)*(g1/i3(l)))
      y=abs(y)
      z=dmax1(dble(y),z)
  99  continue
      if(z.lt.0.0001) then
c
c     successful test (geometrical serie)
c
      do 606 l=-mu,mu
      y1=1.
      d1=inm1(l)
      g1=in(l)
      if(d1.le.accu) go to 606
      y1=1-g1/d1
      if(abs(g1-d1).le.accu) then
      go to 606
      endif
      g1=g1/y1
      i3(l)=i3(l)+g1
  606 continue
      d1=roavion1
      g1=roavion0
      y1=1.
      if(d1.ge.accu) then
      if(abs(g1-d1).ge.accu) then
      y1=1-g1/d1
      g1=g1/y1
      endif
      roavion=roavion+g1
      endif
      go to 505
      endif
c
c     inm2 is the (n-2)ieme scattering order
c
      do 26 k=-mu,mu
      inm2(k)=inm1(k)
   26 continue
      roavion2=roavion1
      endif
c
c     inm1 is the (n-1)ieme scattering order
c
      do 27 k=-mu,mu
      inm1(k)=in(k)
   27 continue
      roavion1=roavion0
c
c     sum of the n-1 orders
c
      do 610 l=-mu,mu
      i3(l)=i3(l)+in(l)
  610 continue
      roavion=roavion+roavion0
c
c     stop if order n is less than 1% of the sum
c
      z=0.
      do 611 l=-mu,mu
      if (abs(i3(l)).ge.accu) then
      y=abs(in(l)/i3(l))
      z=dmax1(z,dble(y))
      endif
  611 continue
      if(z.lt.0.00001) go to 505
c
c      stop if order n is greater than 20 in any case
c
      if(ig-20) 503,503,505
  505 continue
c
c     sum of the fourier component s
c
      delta0s=1
      if(is.ne.0) delta0s=2
      do 612 l=-mu,mu
      i4(l)=i4(l)+delta0s*i3(l)
  612 continue
c
c     stop of the fourier decomposition
c
      do 614 l=1,np
      phi=rp(l)
      do 614 m=-mum1,mum1
      if(m.gt.0) then
      xl(m,l)=xl(m,l)+delta0s*i3(m)*cos(is*(phi+pi))
      else
      xl(m,l)=xl(m,l)+delta0s*i3(m)*cos(is*phi)
      endif
 614  continue
      if(is.eq.0) then
      do k=1,mum1
      xl(0,1)=xl(0,1)+rm(k)*gb(k)*i3(-k)
      enddo
      endif
      xl(mu,1)=xl(mu,1)+delta0s*i3(mu)*cos(is*(phirad+pi))
      xl(-mu,1)=xl(-mu,1)+delta0s*roavion*cos(is*(phirad+pi))
      z=0.
      do 613 l=-mu,mu
       if (abs(i4(l)).lt.accu) goto 613
      x=abs(i3(l)/i4(l))
      z=dmax1(z,x)
  613 continue
      if(z.gt.0.001) go to 24
      goto 243
   24 continue
  243 continue
      nt=snt
c     write(6,*) 'reflectance ', xl(mu,1)/xmus
      return
      end
      subroutine possol (month,jday,tu,xlon,xlat,
     a                   asol,phi0)
 
      real    tu,xlon,xlat,asol,phi0
      integer month,jday,ia,nojour

c     solar position (zenithal angle asol,azimuthal angle phi0
c                     in degrees)
c     jday is the number of the day in the month
 
      ia = 0
      call day_number(jday,month,ia,nojour)
 
      call  pos_fft (nojour, tu, xlon, xlat, asol, phi0)
 
      if(asol.gt.90) call print_error(
     s 'The sun is not raised')
      return
      end

      subroutine day_number(jday,month,ia,j)
      integer jday, month, ia, j

      if (month.le.2) then
                      j=31*(month-1)+jday
		      return
		      endif
      if (month.gt.8) then
                      j=31*(month-1)-((month-2)/2)-2+jday
		      else
                      j=31*(month-1)-((month-1)/2)-2+jday
		      endif
      if(ia.ne.0 .and. mod(ia,4).eq.0) j=j+1
      return
      end

      subroutine pos_fft (j,tu,xlon,xlat,asol,phi0)
      real    tu, xlat, asol,phi0, tsm, xlon,xla, xj, tet,
     a	      a1, a2, a3, a4, a5, et, tsv, ah, b1, b2, b3, b4,
     a	      b5, b6, b7, delta, amuzero, elev, az, caz, azim, pi2
      integer j
      parameter (pi=3.14159265,fac=pi/180.)
c     solar position (zenithal angle asol,azimuthal angle phi0
c                     in degrees)
c     j is the day number in the year
c
c    mean solar time (heure decimale)
 
      tsm=tu+xlon/15.
      xla=xlat*fac
      xj=float(j)
      tet=2.*pi*xj/365.
 
c    time equation (in mn.dec)
      a1=.000075
      a2=.001868
      a3=.032077
      a4=.014615
      a5=.040849
      et=a1+a2*cos(tet)-a3*sin(tet)-a4*cos(2.*tet)-a5*sin(2.*tet)
      et=et*12.*60./pi
 
c     true solar time
 
      tsv=tsm+et/60.
      tsv=(tsv-12.)
 
c     hour angle
 
      ah=tsv*15.*fac
 
c     solar declination   (in radian)
 
      b1=.006918
      b2=.399912
      b3=.070257
      b4=.006758
      b5=.000907
      b6=.002697
      b7=.001480
      delta=b1-b2*cos(tet)+b3*sin(tet)-b4*cos(2.*tet)+b5*sin(2.*tet)-
     &b6*cos(3.*tet)+b7*sin(3.*tet)
 
c     elevation,azimuth
 
      amuzero=sin(xla)*sin(delta)+cos(xla)*cos(delta)*cos(ah)
      elev=asin(amuzero)
      az=cos(delta)*sin(ah)/cos(elev)
      if ( (abs(az)-1.000).gt.0.00000) az = sign(1.,az)
      caz=(-cos(xla)*sin(delta)+sin(xla)*cos(delta)*cos(ah))/cos(elev)
      azim=asin(az)
      if(caz.le.0.) azim=pi-azim
      if(caz.gt.0.and.az.le.0) azim=2*pi+azim
      azim=azim+pi
      pi2=2*pi
      if(azim.gt.pi2) azim=azim-pi2
      elev=elev*180./pi
 
c     conversion in degrees
 
      asol=90.-elev
      phi0=azim/fac
      return
      end
       subroutine presplane(uw,uo3,xpp,ftray)
       real z,p,t,wh,wo,zpl,ppl,tpl,whpl,wopl,xa,xb,xalt
       real xtemp,xwo,xwh,g,air,ro3,rt,rp,roair,ds
       integer i,isup,iinf,k
       common /sixs_atm/z(34),p(34),t(34),wh(34),wo(34)
       common /sixs_planesim/zpl(34),ppl(34),tpl(34),whpl(34),wopl(34)
       real rmo3(34),rmwh(34)
       real ps,xpp,uo3,uw,ftray

C--      print*,'just get into PRESPLANE.f ...'
C--      do k=1,33
C--      write(6,301) z(k),p(k),t(k),wh(k),wo(k),k
C-- 301  format(1x,'z=',e11.4,1x,'p=',e11.4,1x,'t=',e11.4,1x,
C--     @ 'wh=',e11.4,1x,'wo=',e11.4, 'k=',i3)
C--      end do

C--      print*,'before interpolating with PRESPLANE.f ...'
C--      do k=1,34
C--      write(6,302) zpl(k),ppl(k),tpl(k),whpl(k),wopl(k),k
C-- 302  format(1x,'zPL=',e11.4,1x,'p=',e11.4,1x,'t=',e11.4,1x,
C--     @ 'wh=',e11.4,1x,'wo=',e11.4, 'k=',i3)
C--      end do

c      log linear interpolation
           xpp=xpp+z(1)
	   if (xpp.ge.100.) xpp=1000.
           i=0
 10        i=i+1
           if (z(i).le.xpp) goto 10
           isup=i
           iinf=i-1
           xa=(z(isup)-z(iinf))/alog(p(isup)/p(iinf))
           xb=z(isup)-xa*alog(p(isup))
           ps=exp((xpp-xb)/xa)
c interpolating temperature wator vapor and ozone profile versus altitud
	   xalt=xpp
	   xtemp=(t(isup)-t(iinf))/(z(isup)-z(iinf))
	   xtemp=xtemp*(xalt-z(iinf))+t(iinf)
	   xwo=(wo(isup)-wo(iinf))/(z(isup)-z(iinf))
	   xwo=xwo*(xalt-z(iinf))+wo(iinf)
	   xwh=(wh(isup)-wh(iinf))/(z(isup)-z(iinf))
	   xwh=xwh*(xalt-z(iinf))+wh(iinf)
c uptading atmospheric profile
c  last level: plane     , complete to 34
c  with interpolated layers
      do i=1,iinf
      zpl(i)=z(i)
      ppl(i)=p(i)
      tpl(i)=t(i)
      whpl(i)=wh(i)
      wopl(i)=wo(i)
      enddo
      zpl(iinf+1)=xalt
      ppl(iinf+1)=ps
      tpl(iinf+1)=xtemp
      whpl(iinf+1)=xwh
      wopl(iinf+1)=xwo
      do i=iinf+2,34
      zpl(i)=zpl(iinf+1)
      ppl(i)=ppl(iinf+1)
      tpl(i)=tpl(iinf+1)
      whpl(i)=whpl(iinf+1)
      wopl(i)=wopl(iinf+1)
      enddo
c compute modified h2o and o3 integrated content
c compute conversion factor for rayleigh optical thickness computation
c ftray=rp/rt
      uw=0.
      uo3=0.
      g=98.1
      air=0.028964/0.0224
      ro3=0.048/0.0224
      rt=0.
      rp=0.
      do k=1,33
      roair=air*273.16*ppl(k)/(1013.25*tpl(k))
      rmwh(k)=wh(k)/(roair*1000.)
      rmo3(k)=wo(k)/(roair*1000.)
      rt=rt+(p(k+1)/t(k+1)+p(k)/t(k))*(z(k+1)-z(k))
      rp=rp+(ppl(k+1)/tpl(k+1)+ppl(k)/tpl(k))*(zpl(k+1)-zpl(k))
      enddo
      ftray=rp/rt
      do k=2,33
      ds=(ppl(k-1)-ppl(k))/ppl(1)
      uw=uw+((rmwh(k)+rmwh(k-1))/2.)*ds
      uo3=uo3+((rmo3(k)+rmo3(k-1))/2.)*ds
      enddo
      uw=uw*ppl(1)*100./g
      uo3=uo3*ppl(1)*100./g
      uo3=1000.*uo3/ro3

C--      print*,'after interpolating with PRESPLANE.f ...'
C--      do k=1,34
C--      write(6,302) zpl(k),ppl(k),tpl(k),whpl(k),wopl(k),k
C--      end do

      return
      end
       subroutine pressure(uw,uo3,xps)
       common /sixs_atm/z(34),p(34),t(34),wh(34),wo(34)
       real z,p,t,wh,wo,xa,xb,xalt,xtemp,xwo,xwh,g
       real air,ro3,roair,ds
       integer i,isup,iinf,l,k
       real rmo3(34),rmwh(34)
       real ps,xps,uo3,uw
c      log linear interpolation
C---           xps=-xps
	   if (xps.ge.100.) xps=99.99
           i=0
 10        i=i+1
           if (z(i).le.xps) goto 10
           isup=i
           iinf=i-1
           xa=(z(isup)-z(iinf))/alog(p(isup)/p(iinf))
           xb=z(isup)-xa*alog(p(isup))
           ps=exp((xps-xb)/xa)
c interpolating temperature wator vapor and ozone profile versus altitud
	   xalt=xps
	   xtemp=(t(isup)-t(iinf))/(z(isup)-z(iinf))
	   xtemp=xtemp*(xalt-z(iinf))+t(iinf)
	   xwo=(wo(isup)-wo(iinf))/(z(isup)-z(iinf))
	   xwo=xwo*(xalt-z(iinf))+wo(iinf)
	   xwh=(wh(isup)-wh(iinf))/(z(isup)-z(iinf))
	   xwh=xwh*(xalt-z(iinf))+wh(iinf)
c uptading atmospheric profile
c  1rst level: target     , complete to 34
c  with interpolated layers
      z(1)=xalt
      p(1)=ps
      t(1)=xtemp
      wh(1)=xwh
      wo(1)=xwo
      do i=2,33-iinf+1
      z(i)=z(i+iinf-1)
      p(i)=p(i+iinf-1)
      t(i)=t(i+iinf-1)
      wh(i)=wh(i+iinf-1)
      wo(i)=wo(i+iinf-1)
      enddo
      l=33-iinf+1
      do i=l+1,34
      z(i)=(z(34)-z(l))*(i-l)/(34-l)+z(l)
      p(i)=(p(34)-p(l))*(i-l)/(34-l)+p(l)
      t(i)=(t(34)-t(l))*(i-l)/(34-l)+t(l)
      wh(i)=(wh(34)-wh(l))*(i-l)/(34-l)+wh(l)
      wo(i)=(wo(34)-wo(l))*(i-l)/(34-l)+wo(l)
      enddo
c compute modified h2o and o3 integrated content
      uw=0.
      uo3=0.
      g=98.1
      air=0.028964/0.0224
      ro3=0.048/0.0224
      do k=1,33
      roair=air*273.16*p(k)/(1013.25*t(k))
      rmwh(k)=wh(k)/(roair*1000.)
      rmo3(k)=wo(k)/(roair*1000.)
      enddo
      do k=2,33
      ds=(p(k-1)-p(k))/p(1)
      uw=uw+((rmwh(k)+rmwh(k-1))/2.)*ds
      uo3=uo3+((rmo3(k)+rmo3(k-1))/2.)*ds
      enddo
      uw=uw*p(1)*100./g
      uo3=uo3*p(1)*100./g
      uo3=1000.*uo3/ro3
      return
      end
      subroutine print_error(tex)
      character *(*) tex
      logical ier
      integer iwr
      common/sixs_ier/iwr,ier
      ier = .TRUE.
      write(iwr,'(a)')tex
      return
      end
      subroutine scatra (taer,taerp,tray,trayp,piza,
     a      palt,nt,mu,rm,gb,xmus,xmuv,
     a                   ddirtt,ddiftt,udirtt,udiftt,sphalbt,
     a                   ddirtr,ddiftr,udirtr,udiftr,sphalbr,
     a                   ddirta,ddifta,udirta,udifta,sphalba)
 
      integer mu
      real rm(-mu:mu),gb(-mu:mu)
c     computations of the direct and diffuse transmittances
c     for downward and upward paths , and spherical albedo
      real xtrans(-1:1)
      real taer,taerp,tray,trayp,piza,palt,xmus,xmuv
      real udiftt,sphalbt,ddirtr,ddiftr,udirtr,udiftr,sphalbr
      real ddirtt,ddiftt,udirtt,ddirta,ddifta,udirta,udifta
      real sphalba,tamol,tamolp
      integer nt,it
c
      ddirtt=1.
      ddiftt=0.
      udirtt=1.
      udiftt=0.
      sphalbt=0.
      ddirtr=1.
      ddiftr=0.
      udirtr=1.
      udiftr=0.
      sphalbr=0.
      ddirta=1.
      ddifta=0.
      udirta=1.
      udifta=0.
      sphalba=0.
 
      do 1 it=1,3
c it=1 rayleigh only, it=2 aerosol only, it=3 rayleigh+aerosol
      if (it.eq.2.and.taer.le.0.) goto 1
c     compute upward,downward diffuse transmittance for rayleigh,aerosol
      if (it.eq.1) then
         if (palt.gt.900) then
         udiftt=(2./3.+xmuv)+(2./3.-xmuv)*exp(-tray/xmuv)
         udiftt=udiftt/((4./3.)+tray)-exp(-tray/xmuv)
         ddiftt=(2./3.+xmus)+(2./3.-xmus)*exp(-tray/xmus)
         ddiftt=ddiftt/((4./3.)+tray)-exp(-tray/xmus)
         ddirtt=exp(-tray/xmus)
         udirtt=exp(-tray/xmuv)
         call csalbr(tray,sphalbt)
         endif
         if (palt.lt.900) then
           tamol=0.
           tamolp=0.
           rm(-mu)=-xmuv
           rm(mu)=xmuv
           rm(0)=xmus
         call iso(tamol,tray,piza,tamolp,trayp,palt,
     a         nt,mu,rm,gb,xtrans)
         udiftt=xtrans(-1)-exp(-trayp/xmuv)
         udirtt=exp(-trayp/xmuv)
         rm(-mu)=-xmus
         rm(mu)=xmus
         rm(0)=xmus
         ddiftt=(2./3.+xmus)+(2./3.-xmus)*exp(-tray/xmus)
         ddiftt=ddiftt/((4./3.)+tray)-exp(-tray/xmus)
         ddirtt=exp(-tray/xmus)
         udirtt=exp(-tray/xmuv)
         call csalbr(tray,sphalbt)
         endif
         if (palt.le.0.) then
            udiftt=0.
            udirtt=1.
         endif
      endif
      if (it.eq.2) then
      tamol=0.
      tamolp=0.
      rm(-mu)=-xmuv
      rm(mu)=xmuv
      rm(0)=xmus
      call iso(taer,tamol,piza,taerp,tamolp,palt,
     a         nt,mu,rm,gb,xtrans)
      udiftt=xtrans(-1)-exp(-taerp/xmuv)
      udirtt=exp(-taerp/xmuv)
      rm(-mu)=-xmus
      rm(mu)=xmus
      rm(0)=xmus
      call iso(taer,tamol,piza,taerp,tamolp,999.,
     a         nt,mu,rm,gb,xtrans)
      ddirtt=exp(-taer/xmus)
      ddiftt=xtrans(1)-exp(-taer/xmus)
      sphalbt=xtrans(0)*2.
      if (palt.le.0.) then
         udiftt=0.
         udirtt=1.
       endif
      endif
      if (it.eq.3) then
      rm(-mu)=-xmuv
      rm(mu)=xmuv
      rm(0)=xmus
      call iso(taer,tray,piza,taerp,trayp,palt,
     a         nt,mu,rm,gb,xtrans)
      udirtt=exp(-(taerp+trayp)/xmuv)
      udiftt=xtrans(-1)-exp(-(taerp+trayp)/xmuv)
      rm(-mu)=-xmus
      rm(mu)=xmus
      rm(0)=xmus
      call iso(taer,tray,piza,taerp,trayp,999.,
     a         nt,mu,rm,gb,xtrans)
      ddiftt=xtrans(1)-exp(-(taer+tray)/xmus)
      ddirtt=exp(-(taer+tray)/xmus)
      sphalbt=xtrans(0)*2.
      if (palt.le.0.) then
         udiftt=0.
         udirtt=1.
       endif
      endif
c     write(6,*) ddirtt,ddiftt,it,tray,taer,trayp,taerp
 
      if (it.eq.2) goto 2
      if (it.eq.3) goto 1
      ddirtr=ddirtt
      ddiftr=ddiftt
      udirtr=udirtt
      udiftr=udiftt
      sphalbr=sphalbt
      goto 1
    2 ddirta=ddirtt
      ddifta=ddiftt
      udirta=udirtt
      udifta=udiftt
      sphalba=sphalbt
    1 continue
      return
      end
      subroutine   soot
      real ph,phr
      integer i,j
      common /sixs_aerbas/ ph(10,83)
      dimension phr(10,83)
c
c    model: soot
c
            DATA ((PHR(I,J),J=1,83),I=01,01) /
     *0.4897E+00,0.4896E+00,0.4890E+00,0.4881E+00,0.4867E+00,0.4849E+00,
     *0.4827E+00,0.4802E+00,0.4773E+00,0.4743E+00,0.4709E+00,0.4675E+00,
     *0.4638E+00,0.4601E+00,0.4563E+00,0.4526E+00,0.4489E+00,0.4453E+00,
     *0.4419E+00,0.4388E+00,0.4359E+00,0.4334E+00,0.4312E+00,0.4296E+00,
     *0.4285E+00,0.4281E+00,0.4283E+00,0.4293E+00,0.4312E+00,0.4341E+00,
     *0.4380E+00,0.4430E+00,0.4494E+00,0.4571E+00,0.4663E+00,0.4771E+00,
     *0.4896E+00,0.5041E+00,0.5206E+00,0.5392E+00,0.5603E+00,0.5717E+00,
     *0.5838E+00,0.6101E+00,0.6392E+00,0.6714E+00,0.7069E+00,0.7459E+00,
     *0.7886E+00,0.8352E+00,0.8860E+00,0.9411E+00,0.1001E+01,0.1065E+01,
     *0.1135E+01,0.1210E+01,0.1290E+01,0.1376E+01,0.1468E+01,0.1566E+01,
     *0.1670E+01,0.1781E+01,0.1897E+01,0.2019E+01,0.2148E+01,0.2282E+01,
     *0.2421E+01,0.2565E+01,0.2713E+01,0.2865E+01,0.3019E+01,0.3173E+01,
     *0.3327E+01,0.3479E+01,0.3625E+01,0.3765E+01,0.3894E+01,0.4011E+01,
     *0.4111E+01,0.4192E+01,0.4250E+01,0.4284E+01,0.4292E+01/
            DATA ((PHR(I,J),J=1,83),I=02,02) /
     *0.5620E+00,0.5618E+00,0.5611E+00,0.5599E+00,0.5582E+00,0.5560E+00,
     *0.5533E+00,0.5502E+00,0.5467E+00,0.5428E+00,0.5387E+00,0.5342E+00,
     *0.5295E+00,0.5246E+00,0.5197E+00,0.5146E+00,0.5096E+00,0.5046E+00,
     *0.4998E+00,0.4951E+00,0.4907E+00,0.4866E+00,0.4829E+00,0.4797E+00,
     *0.4771E+00,0.4751E+00,0.4738E+00,0.4734E+00,0.4738E+00,0.4753E+00,
     *0.4779E+00,0.4817E+00,0.4868E+00,0.4934E+00,0.5016E+00,0.5114E+00,
     *0.5231E+00,0.5367E+00,0.5524E+00,0.5704E+00,0.5908E+00,0.6019E+00,
     *0.6137E+00,0.6393E+00,0.6678E+00,0.6993E+00,0.7340E+00,0.7720E+00,
     *0.8136E+00,0.8589E+00,0.9081E+00,0.9613E+00,0.1019E+01,0.1080E+01,
     *0.1147E+01,0.1218E+01,0.1293E+01,0.1373E+01,0.1459E+01,0.1549E+01,
     *0.1643E+01,0.1743E+01,0.1847E+01,0.1956E+01,0.2069E+01,0.2185E+01,
     *0.2305E+01,0.2428E+01,0.2553E+01,0.2679E+01,0.2806E+01,0.2931E+01,
     *0.3055E+01,0.3174E+01,0.3289E+01,0.3396E+01,0.3495E+01,0.3582E+01,
     *0.3656E+01,0.3716E+01,0.3758E+01,0.3782E+01,0.3788E+01/
            DATA ((PHR(I,J),J=1,83),I=03,03) /
     *0.5834E+00,0.5832E+00,0.5825E+00,0.5813E+00,0.5795E+00,0.5771E+00,
     *0.5743E+00,0.5710E+00,0.5673E+00,0.5632E+00,0.5587E+00,0.5540E+00,
     *0.5490E+00,0.5438E+00,0.5384E+00,0.5330E+00,0.5275E+00,0.5221E+00,
     *0.5168E+00,0.5117E+00,0.5068E+00,0.5023E+00,0.4981E+00,0.4944E+00,
     *0.4913E+00,0.4889E+00,0.4871E+00,0.4862E+00,0.4862E+00,0.4872E+00,
     *0.4894E+00,0.4928E+00,0.4975E+00,0.5037E+00,0.5115E+00,0.5210E+00,
     *0.5324E+00,0.5457E+00,0.5611E+00,0.5788E+00,0.5988E+00,0.6098E+00,
     *0.6215E+00,0.6468E+00,0.6749E+00,0.7061E+00,0.7405E+00,0.7781E+00,
     *0.8193E+00,0.8641E+00,0.9127E+00,0.9652E+00,0.1022E+01,0.1083E+01,
     *0.1148E+01,0.1217E+01,0.1291E+01,0.1370E+01,0.1453E+01,0.1541E+01,
     *0.1633E+01,0.1730E+01,0.1831E+01,0.1936E+01,0.2045E+01,0.2157E+01,
     *0.2272E+01,0.2390E+01,0.2509E+01,0.2629E+01,0.2749E+01,0.2867E+01,
     *0.2984E+01,0.3096E+01,0.3203E+01,0.3304E+01,0.3395E+01,0.3476E+01,
     *0.3545E+01,0.3599E+01,0.3638E+01,0.3660E+01,0.3666E+01/
            DATA ((PHR(I,J),J=1,83),I=04,04) /
     *0.6060E+00,0.6059E+00,0.6051E+00,0.6038E+00,0.6019E+00,0.5994E+00,
     *0.5964E+00,0.5929E+00,0.5889E+00,0.5846E+00,0.5798E+00,0.5747E+00,
     *0.5693E+00,0.5637E+00,0.5580E+00,0.5521E+00,0.5462E+00,0.5403E+00,
     *0.5345E+00,0.5289E+00,0.5235E+00,0.5185E+00,0.5138E+00,0.5096E+00,
     *0.5059E+00,0.5029E+00,0.5007E+00,0.4993E+00,0.4988E+00,0.4993E+00,
     *0.5010E+00,0.5040E+00,0.5083E+00,0.5142E+00,0.5216E+00,0.5307E+00,
     *0.5418E+00,0.5548E+00,0.5699E+00,0.5873E+00,0.6071E+00,0.6180E+00,
     *0.6295E+00,0.6546E+00,0.6825E+00,0.7134E+00,0.7474E+00,0.7848E+00,
     *0.8255E+00,0.8699E+00,0.9179E+00,0.9698E+00,0.1026E+01,0.1085E+01,
     *0.1150E+01,0.1218E+01,0.1290E+01,0.1367E+01,0.1448E+01,0.1534E+01,
     *0.1623E+01,0.1717E+01,0.1815E+01,0.1916E+01,0.2020E+01,0.2128E+01,
     *0.2237E+01,0.2349E+01,0.2462E+01,0.2576E+01,0.2688E+01,0.2800E+01,
     *0.2909E+01,0.3013E+01,0.3113E+01,0.3206E+01,0.3290E+01,0.3364E+01,
     *0.3427E+01,0.3477E+01,0.3512E+01,0.3532E+01,0.3537E+01/
            DATA ((PHR(I,J),J=1,83),I=05,05) /
     *0.6604E+00,0.6602E+00,0.6593E+00,0.6578E+00,0.6556E+00,0.6528E+00,
     *0.6494E+00,0.6454E+00,0.6409E+00,0.6358E+00,0.6304E+00,0.6245E+00,
     *0.6182E+00,0.6117E+00,0.6050E+00,0.5981E+00,0.5911E+00,0.5841E+00,
     *0.5771E+00,0.5703E+00,0.5636E+00,0.5573E+00,0.5513E+00,0.5458E+00,
     *0.5409E+00,0.5366E+00,0.5331E+00,0.5305E+00,0.5288E+00,0.5281E+00,
     *0.5287E+00,0.5305E+00,0.5338E+00,0.5385E+00,0.5450E+00,0.5532E+00,
     *0.5633E+00,0.5754E+00,0.5897E+00,0.6062E+00,0.6252E+00,0.6356E+00,
     *0.6467E+00,0.6710E+00,0.6980E+00,0.7280E+00,0.7610E+00,0.7972E+00,
     *0.8367E+00,0.8797E+00,0.9261E+00,0.9762E+00,0.1030E+01,0.1087E+01,
     *0.1149E+01,0.1214E+01,0.1283E+01,0.1355E+01,0.1432E+01,0.1512E+01,
     *0.1595E+01,0.1682E+01,0.1772E+01,0.1865E+01,0.1961E+01,0.2058E+01,
     *0.2157E+01,0.2257E+01,0.2358E+01,0.2458E+01,0.2557E+01,0.2654E+01,
     *0.2748E+01,0.2838E+01,0.2923E+01,0.3001E+01,0.3072E+01,0.3134E+01,
     *0.3187E+01,0.3228E+01,0.3257E+01,0.3273E+01,0.3277E+01/
            DATA ((PHR(I,J),J=1,83),I=06,06) /
     *0.6993E+00,0.6991E+00,0.6982E+00,0.6965E+00,0.6942E+00,0.6911E+00,
     *0.6874E+00,0.6830E+00,0.6781E+00,0.6726E+00,0.6666E+00,0.6601E+00,
     *0.6533E+00,0.6461E+00,0.6387E+00,0.6310E+00,0.6232E+00,0.6154E+00,
     *0.6076E+00,0.5998E+00,0.5923E+00,0.5851E+00,0.5782E+00,0.5717E+00,
     *0.5659E+00,0.5607E+00,0.5562E+00,0.5526E+00,0.5500E+00,0.5485E+00,
     *0.5482E+00,0.5491E+00,0.5515E+00,0.5555E+00,0.5611E+00,0.5686E+00,
     *0.5779E+00,0.5893E+00,0.6028E+00,0.6187E+00,0.6369E+00,0.6470E+00,
     *0.6577E+00,0.6812E+00,0.7074E+00,0.7366E+00,0.7687E+00,0.8040E+00,
     *0.8425E+00,0.8843E+00,0.9295E+00,0.9781E+00,0.1030E+01,0.1086E+01,
     *0.1145E+01,0.1208E+01,0.1274E+01,0.1344E+01,0.1417E+01,0.1494E+01,
     *0.1573E+01,0.1656E+01,0.1741E+01,0.1828E+01,0.1918E+01,0.2009E+01,
     *0.2101E+01,0.2194E+01,0.2287E+01,0.2380E+01,0.2470E+01,0.2559E+01,
     *0.2645E+01,0.2726E+01,0.2803E+01,0.2873E+01,0.2937E+01,0.2992E+01,
     *0.3038E+01,0.3075E+01,0.3100E+01,0.3115E+01,0.3118E+01/
            DATA ((PHR(I,J),J=1,83),I=07,07) /
     *0.7916E+00,0.7914E+00,0.7903E+00,0.7883E+00,0.7855E+00,0.7818E+00,
     *0.7773E+00,0.7721E+00,0.7662E+00,0.7595E+00,0.7522E+00,0.7444E+00,
     *0.7360E+00,0.7272E+00,0.7180E+00,0.7085E+00,0.6988E+00,0.6889E+00,
     *0.6790E+00,0.6692E+00,0.6595E+00,0.6500E+00,0.6408E+00,0.6321E+00,
     *0.6239E+00,0.6164E+00,0.6097E+00,0.6038E+00,0.5989E+00,0.5952E+00,
     *0.5926E+00,0.5915E+00,0.5918E+00,0.5936E+00,0.5972E+00,0.6027E+00,
     *0.6101E+00,0.6195E+00,0.6311E+00,0.6451E+00,0.6614E+00,0.6705E+00,
     *0.6803E+00,0.7017E+00,0.7259E+00,0.7529E+00,0.7828E+00,0.8156E+00,
     *0.8514E+00,0.8903E+00,0.9323E+00,0.9774E+00,0.1026E+01,0.1077E+01,
     *0.1131E+01,0.1189E+01,0.1249E+01,0.1312E+01,0.1378E+01,0.1447E+01,
     *0.1518E+01,0.1590E+01,0.1665E+01,0.1741E+01,0.1819E+01,0.1897E+01,
     *0.1976E+01,0.2054E+01,0.2132E+01,0.2209E+01,0.2284E+01,0.2356E+01,
     *0.2426E+01,0.2491E+01,0.2552E+01,0.2607E+01,0.2657E+01,0.2700E+01,
     *0.2736E+01,0.2764E+01,0.2783E+01,0.2795E+01,0.2797E+01/
            DATA ((PHR(I,J),J=1,83),I=08,08) /
     *0.1041E+01,0.1040E+01,0.1038E+01,0.1036E+01,0.1031E+01,0.1026E+01,
     *0.1019E+01,0.1011E+01,0.1002E+01,0.9924E+00,0.9814E+00,0.9694E+00,
     *0.9566E+00,0.9431E+00,0.9288E+00,0.9140E+00,0.8988E+00,0.8832E+00,
     *0.8673E+00,0.8513E+00,0.8353E+00,0.8194E+00,0.8038E+00,0.7885E+00,
     *0.7737E+00,0.7596E+00,0.7462E+00,0.7338E+00,0.7223E+00,0.7121E+00,
     *0.7031E+00,0.6955E+00,0.6895E+00,0.6852E+00,0.6827E+00,0.6820E+00,
     *0.6833E+00,0.6868E+00,0.6924E+00,0.7003E+00,0.7105E+00,0.7165E+00,
     *0.7232E+00,0.7383E+00,0.7559E+00,0.7760E+00,0.7987E+00,0.8240E+00,
     *0.8518E+00,0.8821E+00,0.9149E+00,0.9501E+00,0.9877E+00,0.1028E+01,
     *0.1069E+01,0.1113E+01,0.1159E+01,0.1207E+01,0.1256E+01,0.1306E+01,
     *0.1358E+01,0.1410E+01,0.1463E+01,0.1517E+01,0.1570E+01,0.1623E+01,
     *0.1676E+01,0.1727E+01,0.1778E+01,0.1827E+01,0.1873E+01,0.1918E+01,
     *0.1960E+01,0.1999E+01,0.2035E+01,0.2067E+01,0.2096E+01,0.2120E+01,
     *0.2140E+01,0.2156E+01,0.2167E+01,0.2173E+01,0.2174E+01/
            DATA ((PHR(I,J),J=1,83),I=09,09) /
     *0.1182E+01,0.1181E+01,0.1179E+01,0.1176E+01,0.1171E+01,0.1164E+01,
     *0.1156E+01,0.1147E+01,0.1136E+01,0.1124E+01,0.1110E+01,0.1096E+01,
     *0.1080E+01,0.1064E+01,0.1046E+01,0.1028E+01,0.1009E+01,0.9903E+00,
     *0.9708E+00,0.9510E+00,0.9312E+00,0.9114E+00,0.8919E+00,0.8726E+00,
     *0.8539E+00,0.8357E+00,0.8184E+00,0.8019E+00,0.7866E+00,0.7724E+00,
     *0.7595E+00,0.7481E+00,0.7383E+00,0.7302E+00,0.7239E+00,0.7195E+00,
     *0.7171E+00,0.7168E+00,0.7188E+00,0.7229E+00,0.7294E+00,0.7335E+00,
     *0.7382E+00,0.7494E+00,0.7630E+00,0.7790E+00,0.7974E+00,0.8182E+00,
     *0.8414E+00,0.8668E+00,0.8944E+00,0.9242E+00,0.9561E+00,0.9898E+00,
     *0.1025E+01,0.1063E+01,0.1101E+01,0.1141E+01,0.1183E+01,0.1225E+01,
     *0.1268E+01,0.1311E+01,0.1355E+01,0.1399E+01,0.1442E+01,0.1485E+01,
     *0.1528E+01,0.1569E+01,0.1609E+01,0.1648E+01,0.1685E+01,0.1720E+01,
     *0.1753E+01,0.1783E+01,0.1811E+01,0.1836E+01,0.1858E+01,0.1876E+01,
     *0.1891E+01,0.1903E+01,0.1911E+01,0.1916E+01,0.1917E+01/
            DATA ((PHR(I,J),J=1,83),I=10,10) /
     *0.1325E+01,0.1324E+01,0.1322E+01,0.1318E+01,0.1312E+01,0.1304E+01,
     *0.1294E+01,0.1283E+01,0.1270E+01,0.1256E+01,0.1240E+01,0.1222E+01,
     *0.1204E+01,0.1184E+01,0.1163E+01,0.1142E+01,0.1119E+01,0.1096E+01,
     *0.1073E+01,0.1049E+01,0.1025E+01,0.1001E+01,0.9776E+00,0.9541E+00,
     *0.9312E+00,0.9088E+00,0.8872E+00,0.8666E+00,0.8471E+00,0.8287E+00,
     *0.8118E+00,0.7963E+00,0.7825E+00,0.7704E+00,0.7602E+00,0.7519E+00,
     *0.7457E+00,0.7415E+00,0.7396E+00,0.7399E+00,0.7424E+00,0.7446E+00,
     *0.7473E+00,0.7545E+00,0.7640E+00,0.7758E+00,0.7899E+00,0.8063E+00,
     *0.8248E+00,0.8455E+00,0.8681E+00,0.8928E+00,0.9192E+00,0.9473E+00,
     *0.9771E+00,0.1008E+01,0.1041E+01,0.1074E+01,0.1109E+01,0.1144E+01,
     *0.1179E+01,0.1215E+01,0.1252E+01,0.1288E+01,0.1324E+01,0.1359E+01,
     *0.1393E+01,0.1427E+01,0.1460E+01,0.1491E+01,0.1521E+01,0.1549E+01,
     *0.1575E+01,0.1599E+01,0.1622E+01,0.1641E+01,0.1658E+01,0.1673E+01,
     *0.1685E+01,0.1694E+01,0.1701E+01,0.1704E+01,0.1705E+01/
      do 1 i=1,10
      do 1 j=1,83
      ph(i,j)=phr(i,j)
    1 continue
      return
      end
      subroutine specinterp(wl,taer55,taer55p,
     s     tamoy,tamoyp,pizmoy,pizmoyp)
      real wl,taer55,taer55p,tamoy,tamoyp,pizmoy,pizmoyp,roatm
      real dtdir,dtdif,utdir,utdif,sphal,wldis,trayl,traypl
      real ext,ome,gasym,phase,pha,betal,phasel,cgaus,pdgs,coef
      real wlinf,alphaa,betaa,tsca,coeff
      integer linf,ll,lsup,k
      common /sixs_disc/ roatm(3,10),dtdir(3,10),dtdif(3,10),
     s utdir(3,10),utdif(3,10),sphal(3,10),wldis(10),trayl(10),
     s traypl(10)
      common /sixs_aer/ext(10),ome(10),gasym(10),phase(10)
      common /sixs_trunc/pha(83),betal(0:80)
      common /sixs_sos/phasel(10,83),cgaus(83),pdgs(83)
      linf=1
      do 80 ll=1,9
      if(wl.ge.wldis(ll).and.wl.le.wldis(ll+1)) linf=ll
   80 continue
      if(wl.gt.wldis(10)) linf=9
      lsup=linf+1
      coef=alog(wldis(lsup)/wldis(linf))
      wlinf=wldis(linf)
      alphaa=alog(ext(lsup)*ome(lsup)/(ext(linf)*ome(linf)))/coef
      betaa=ext(linf)*ome(linf)/(wlinf**(alphaa))
      tsca=taer55*betaa*(wl**alphaa)/ext(4)
      alphaa=alog(ext(lsup)/(ext(linf)))/coef
      betaa=ext(linf)/(wlinf**(alphaa))
      tamoy=taer55*betaa*(wl**alphaa)/ext(4)
      tamoyp=taer55p*betaa*(wl**alphaa)/ext(4)
      pizmoy=tsca/tamoy
      pizmoyp=pizmoy
      do 81 k=1,83
      alphaa=alog(phasel(lsup,k)/phasel(linf,k))/coef
      betaa=phasel(linf,k)/(wlinf**(alphaa))
 81   pha(k)=betaa*(wl**alphaa)
      call trunca(coeff)
      tamoy=tamoy*(1.-pizmoy*coeff)
      tamoyp=tamoyp*(1.-pizmoyp*coeff)
      pizmoy=pizmoy*(1.-coeff)/(1.-pizmoy*coeff)
      return
      end
      subroutine splie2(x2a,ya,m,n,y2a)

      parameter (nn=100)
      integer m,n,j,k
      real x2a(n),ya(m,n),y2a(m,n),ytmp(nn),y2tmp(nn)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline(x2a,ytmp,n,1.e30,1.e30,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      return
      end
      subroutine splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      parameter (nn=100)
      integer m,n,j,k
      real x1,x2,y
      real x1a(m),x2a(n),ya(m,n),y2a(m,n),ytmp(nn),y2tmp(nn)
      real yytmp(nn)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12    continue
      call spline(x1a,yytmp,m,1.e30,1.e30,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      end

      subroutine  stm
      common/sixs_aerbas/ph(10,83)
      real phr(10,83)
      real ph
      integer i,j
c
c  model: Stratospheric aerosol as follow king's model 
CJournal of Climate and Applied Meteorology, Vol23, No7, pp=1121-1137, 1984
      data ((phr(i,j),j=1,83),i= 1, 1)/
     &  .4482,  .4378,  .3984,  .3460,  .3030,  .2864,  .3011,  .3393,
     &  .3852,  .4224,  .4395,  .4332,  .4068,  .3674,  .3232,  .2806,
     &  .2436,  .2137,  .1909,  .1740,  .1615,  .1523,  .1453,  .1398,
     &  .1356,  .1324,  .1300,  .1284,  .1277,  .1278,  .1286,  .1303,
     &  .1328,  .1362,  .1404,  .1455,  .1515,  .1585,  .1666,  .1759,
     &  .1864,  .1922,  .1984,  .2119,  .2272,  .2444,  .2638,  .2856,
     &  .3103,  .3381,  .3696,  .4052,  .4454,  .4911,  .5429,  .6018,
     &  .6687,  .7447,  .8309,  .9284, 1.0383, 1.1614, 1.2985, 1.4500,
     & 1.6169, 1.8014, 2.0088, 2.2506, 2.5487, 2.9404, 3.4830, 4.2562,
     & 5.3583, 6.8944, 8.9537,11.5772,14.7221,18.2338,21.8390,25.1693,
     &27.8195,29.4297,29.8220/
      data ((phr(i,j),j=1,83),i= 2, 2)/
     &  .3066,  .3025,  .2862,  .2621,  .2369,  .2173,  .2078,  .2095,
     &  .2201,  .2355,  .2504,  .2607,  .2637,  .2589,  .2472,  .2305,
     &  .2114,  .1919,  .1736,  .1577,  .1445,  .1340,  .1261,  .1203,
     &  .1162,  .1134,  .1117,  .1109,  .1108,  .1113,  .1124,  .1141,
     &  .1165,  .1194,  .1230,  .1273,  .1324,  .1384,  .1452,  .1531,
     &  .1620,  .1669,  .1722,  .1838,  .1969,  .2117,  .2285,  .2475,
     &  .2691,  .2936,  .3213,  .3528,  .3886,  .4293,  .4754,  .5278,
     &  .5872,  .6543,  .7304,  .8164,  .9142, 1.0260, 1.1554, 1.3080,
     & 1.4922, 1.7208, 2.0120, 2.3907, 2.8891, 3.5464, 4.4062, 5.5124,
     & 6.9014, 8.5929,10.5796,12.8175,15.2199,17.6577,19.9678,21.9699,
     &23.4901,24.3864,24.6019/
      data ((phr(i,j),j=1,83),i= 3, 3)/
     &  .2797,  .2765,  .2636,  .2440,  .2227,  .2045,  .1934,  .1907,
     &  .1956,  .2056,  .2171,  .2266,  .2316,  .2310,  .2247,  .2138,
     &  .1998,  .1843,  .1689,  .1546,  .1422,  .1319,  .1237,  .1176,
     &  .1131,  .1101,  .1082,  .1073,  .1072,  .1077,  .1089,  .1107,
     &  .1130,  .1160,  .1196,  .1239,  .1289,  .1347,  .1413,  .1490,
     &  .1577,  .1625,  .1676,  .1789,  .1916,  .2061,  .2225,  .2410,
     &  .2621,  .2859,  .3130,  .3437,  .3785,  .4179,  .4626,  .5133,
     &  .5706,  .6356,  .7094,  .7936,  .8904, 1.0031, 1.1367, 1.2984,
     & 1.4985, 1.7518, 2.0779, 2.5018, 3.0542, 3.7695, 4.6834, 5.8280,
     & 7.2258, 8.8822,10.7776,12.8620,15.0519,17.2327,19.2671,21.0083,
     &22.3181,23.0858,23.2698/
      data ((phr(i,j),j=1,83),i= 4, 4)/
     &  .2523,  .2499,  .2401,  .2249,  .2075,  .1914,  .1795,  .1736,
     &  .1735,  .1782,  .1854,  .1928,  .1984,  .2005,  .1988,  .1932,
     &  .1846,  .1739,  .1623,  .1506,  .1398,  .1303,  .1223,  .1159,
     &  .1110,  .1076,  .1054,  .1042,  .1039,  .1044,  .1055,  .1073,
     &  .1097,  .1127,  .1163,  .1205,  .1255,  .1312,  .1378,  .1453,
     &  .1539,  .1586,  .1636,  .1746,  .1871,  .2013,  .2173,  .2354,
     &  .2559,  .2792,  .3055,  .3352,  .3689,  .4070,  .4502,  .4990,
     &  .5545,  .6178,  .6905,  .7747,  .8738,  .9921, 1.1363, 1.3153,
     & 1.5410, 1.8295, 2.2003, 2.6770, 3.2861, 4.0549, 5.0090, 6.1680,
     & 7.5404, 9.1188,10.8752,12.7575,14.6902,16.5769,18.3076,19.7692,
     &20.8579,21.4919,21.6435/
      data ((phr(i,j),j=1,83),i= 5, 5)/
     &  .2099,  .2085,  .2029,  .1937,  .1824,  .1705,  .1597,  .1512,
     &  .1457,  .1433,  .1435,  .1455,  .1484,  .1511,  .1529,  .1533,
     &  .1519,  .1489,  .1445,  .1391,  .1331,  .1270,  .1212,  .1158,
     &  .1112,  .1075,  .1048,  .1029,  .1020,  .1019,  .1027,  .1041,
     &  .1063,  .1092,  .1128,  .1170,  .1220,  .1278,  .1344,  .1419,
     &  .1505,  .1551,  .1601,  .1710,  .1833,  .1971,  .2127,  .2303,
     &  .2501,  .2724,  .2976,  .3260,  .3583,  .3950,  .4371,  .4857,
     &  .5424,  .6092,  .6892,  .7862,  .9053, 1.0531, 1.2379, 1.4701,
     & 1.7619, 2.1272, 2.5813, 3.1398, 3.8174, 4.6261, 5.5735, 6.6598,
     & 7.8763, 9.2034,10.6092,12.0501,13.4719,14.8129,16.0082,16.9948,
     &17.7172,18.1334,18.2325/
      data ((phr(i,j),j=1,83),i= 6, 6)/
     &  .1911,  .1901,  .1861,  .1793,  .1706,  .1610,  .1516,  .1432,
     &  .1365,  .1318,  .1292,  .1284,  .1289,  .1301,  .1316,  .1328,
     &  .1333,  .1330,  .1317,  .1295,  .1266,  .1232,  .1196,  .1160,
     &  .1126,  .1096,  .1072,  .1054,  .1043,  .1040,  .1044,  .1056,
     &  .1075,  .1102,  .1136,  .1177,  .1227,  .1285,  .1351,  .1427,
     &  .1513,  .1560,  .1610,  .1719,  .1842,  .1981,  .2136,  .2311,
     &  .2509,  .2732,  .2986,  .3275,  .3607,  .3992,  .4441,  .4973,
     &  .5608,  .6374,  .7309,  .8458,  .9877, 1.1636, 1.3815, 1.6506,
     & 1.9812, 2.3839, 2.8694, 3.4473, 4.1253, 4.9077, 5.7944, 6.7794,
     & 7.8497, 8.9848,10.1567,11.3301,12.4643,13.5152,14.4381,15.1909,
     &15.7373,16.0504,16.1247/
      data ((phr(i,j),j=1,83),i= 7, 7)/
     &  .1657,  .1652,  .1631,  .1595,  .1546,  .1488,  .1424,  .1358,
     &  .1294,  .1235,  .1183,  .1141,  .1107,  .1084,  .1070,  .1063,
     &  .1062,  .1066,  .1072,  .1080,  .1088,  .1096,  .1103,  .1108,
     &  .1113,  .1117,  .1121,  .1126,  .1133,  .1142,  .1155,  .1172,
     &  .1193,  .1221,  .1255,  .1296,  .1345,  .1402,  .1469,  .1547,
     &  .1636,  .1686,  .1739,  .1856,  .1991,  .2147,  .2326,  .2534,
     &  .2775,  .3058,  .3392,  .3787,  .4256,  .4818,  .5491,  .6299,
     &  .7270,  .8435,  .9830, 1.1494, 1.3469, 1.5800, 1.8530, 2.1701,
     & 2.5350, 2.9507, 3.4187, 3.9394, 4.5111, 5.1299, 5.7894, 6.4806,
     & 7.1921, 7.9098, 8.6176, 9.2978, 9.9320,10.5016,10.9891,11.3786,
     &11.6571,11.8152,11.8525/
      data ((phr(i,j),j=1,83),i= 8, 8)/
     &  .1867,  .1866,  .1860,  .1850,  .1836,  .1819,  .1797,  .1773,
     &  .1746,  .1717,  .1687,  .1655,  .1624,  .1593,  .1563,  .1535,
     &  .1509,  .1487,  .1469,  .1455,  .1447,  .1444,  .1449,  .1460,
     &  .1480,  .1509,  .1547,  .1596,  .1656,  .1729,  .1814,  .1915,
     &  .2031,  .2164,  .2315,  .2488,  .2683,  .2902,  .3149,  .3426,
     &  .3736,  .3904,  .4081,  .4466,  .4894,  .5369,  .5895,  .6476,
     &  .7117,  .7821,  .8593,  .9436, 1.0355, 1.1354, 1.2434, 1.3598,
     & 1.4848, 1.6183, 1.7604, 1.9108, 2.0693, 2.2352, 2.4081, 2.5870,
     & 2.7711, 2.9591, 3.1498, 3.3417, 3.5332, 3.7226, 3.9080, 4.0876,
     & 4.2594, 4.4215, 4.5720, 4.7090, 4.8308, 4.9359, 5.0228, 5.0905,
     & 5.1379, 5.1645, 5.1708/
      data ((phr(i,j),j=1,83),i= 9, 9)/
     &  .4829,  .4828,  .4824,  .4816,  .4804,  .4790,  .4772,  .4751,
     &  .4728,  .4701,  .4673,  .4643,  .4611,  .4578,  .4544,  .4511,
     &  .4477,  .4444,  .4413,  .4384,  .4358,  .4335,  .4317,  .4304,
     &  .4298,  .4299,  .4308,  .4327,  .4356,  .4397,  .4452,  .4520,
     &  .4605,  .4708,  .4829,  .4971,  .5135,  .5323,  .5536,  .5776,
     &  .6045,  .6190,  .6344,  .6674,  .7038,  .7435,  .7869,  .8338,
     &  .8845,  .9390,  .9973, 1.0594, 1.1253, 1.1949, 1.2682, 1.3449,
     & 1.4249, 1.5080, 1.5939, 1.6823, 1.7728, 1.8650, 1.9584, 2.0527,
     & 2.1472, 2.2414, 2.3347, 2.4266, 2.5162, 2.6031, 2.6866, 2.7660,
     & 2.8408, 2.9103, 2.9739, 3.0312, 3.0815, 3.1246, 3.1599, 3.1873,
     & 3.2064, 3.2170, 3.2195/
      data ((phr(i,j),j=1,83),i=10,10)/
     & 1.0488, 1.0485, 1.0470, 1.0443, 1.0405, 1.0355, 1.0295, 1.0223,
     & 1.0141, 1.0049,  .9948,  .9838,  .9719,  .9594,  .9461,  .9323,
     &  .9180,  .9032,  .8882,  .8730,  .8577,  .8425,  .8273,  .8125,
     &  .7981,  .7841,  .7709,  .7584,  .7469,  .7364,  .7271,  .7191,
     &  .7126,  .7077,  .7045,  .7031,  .7036,  .7062,  .7109,  .7179,
     &  .7271,  .7326,  .7387,  .7527,  .7692,  .7881,  .8096,  .8335,
     &  .8599,  .8886,  .9198,  .9532,  .9888, 1.0265, 1.0661, 1.1075,
     & 1.1505, 1.1949, 1.2406, 1.2872, 1.3346, 1.3825, 1.4307, 1.4789,
     & 1.5267, 1.5741, 1.6205, 1.6659, 1.7098, 1.7521, 1.7924, 1.8305,
     & 1.8661, 1.8989, 1.9289, 1.9557, 1.9792, 1.9992, 2.0156, 2.0282,
     & 2.0370, 2.0419, 2.0431/
      do 1 i=1,10
      do 1 j=1,83
      ph(i,j)=phr(i,j)
   1  continue
      return
      end
      subroutine   subsum
      integer i
      real z4(34),p4(34),t4(34),wh4(34),wo4(34)
      real z,p,t,wh,wo
      common /sixs_atm/z(34),p(34),t(34),wh(34),wo(34)
c
c     model: subarctique summer mc clatchey
c
      data(z4(i),i=1, 34)/
     1    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,
     2    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,
     3   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,
     4   35.,   40.,   45.,   50.,   70.,  100.,99999./
      data (p4(i),i=1,34) /
     a1.010e+03,8.960e+02,7.929e+02,7.000e+02,6.160e+02,5.410e+02,
     a4.730e+02,4.130e+02,3.590e+02,3.107e+02,2.677e+02,2.300e+02,
     a1.977e+02,1.700e+02,1.460e+02,1.250e+02,1.080e+02,9.280e+01,
     a7.980e+01,6.860e+01,5.890e+01,5.070e+01,4.360e+01,3.750e+01,
     a3.227e+01,2.780e+01,1.340e+01,6.610e+00,3.400e+00,1.810e+00,
     a9.870e-01,7.070e-02,3.000e-04,0.000e+00/
      data (t4(i),i=1,34) /
     a2.870e+02,2.820e+02,2.760e+02,2.710e+02,2.660e+02,2.600e+02,
     a2.530e+02,2.460e+02,2.390e+02,2.320e+02,2.250e+02,2.250e+02,
     a2.250e+02,2.250e+02,2.250e+02,2.250e+02,2.250e+02,2.250e+02,
     a2.250e+02,2.250e+02,2.250e+02,2.250e+02,2.250e+02,2.250e+02,
     a2.260e+02,2.280e+02,2.350e+02,2.470e+02,2.620e+02,2.740e+02,
     a2.770e+02,2.160e+02,2.100e+02,2.100e+02/
      data (wh4(i),i=1,34) /
     a9.100e+00,6.000e+00,4.200e+00,2.700e+00,1.700e+00,1.000e+00,
     a5.400e-01,2.900e-01,1.300e-01,4.200e-02,1.500e-02,9.400e-03,
     a6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,
     a5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,
     a6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,
     a6.300e-06,1.400e-07,1.000e-09,0.000e+00/
      data (wo4(i),i=1,34) /
     a4.900e-05,5.400e-05,5.600e-05,5.800e-05,6.000e-05,6.400e-05,
     a7.100e-05,7.500e-05,7.900e-05,1.100e-04,1.300e-04,1.800e-04,
     a2.100e-04,2.600e-04,2.800e-04,3.200e-04,3.400e-04,3.900e-04,
     a4.100e-04,4.100e-04,3.900e-04,3.600e-04,3.200e-04,3.000e-04,
     a2.800e-04,2.600e-04,1.400e-04,9.200e-05,4.100e-05,1.300e-05,
     a4.300e-06,8.600e-08,4.300e-11,0.000e+00/
      do 1 i=1,34
      z(i)=z4(i)
      p(i)=p4(i)
      t(i)=t4(i)
      wh(i)=wh4(i)
      wo(i)=wo4(i)
    1 continue
      return
      end
      subroutine   subwin

      real z5(34),p5(34),t5(34),wh5(34),wo5(34)
      real z,p,t,wh,wo
      common /sixs_atm/z(34),p(34),t(34),wh(34),wo(34)
      integer i
c
c     model: subarctique winter mc clatchey
c
      data(z5(i),i=1, 34)/
     1    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,
     2    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,
     3   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,
     4   35.,   40.,   45.,   50.,   70.,  100.,99999./
      data (p5(i),i=1,34) /
     a1.013e+03,8.878e+02,7.775e+02,6.798e+02,5.932e+02,5.158e+02,
     a4.467e+02,3.853e+02,3.308e+02,2.829e+02,2.418e+02,2.067e+02,
     a1.766e+02,1.510e+02,1.291e+02,1.103e+02,9.431e+01,8.058e+01,
     a6.882e+01,5.875e+01,5.014e+01,4.277e+01,3.647e+01,3.109e+01,
     a2.649e+01,2.256e+01,1.020e+01,4.701e+00,2.243e+00,1.113e+00,
     a5.719e-01,4.016e-02,3.000e-04,0.000e+00/
      data (t5(i),i=1,34) /
     a2.571e+02,2.591e+02,2.559e+02,2.527e+02,2.477e+02,2.409e+02,
     a2.341e+02,2.273e+02,2.206e+02,2.172e+02,2.172e+02,2.172e+02,
     a2.172e+02,2.172e+02,2.172e+02,2.172e+02,2.166e+02,2.160e+02,
     a2.154e+02,2.148e+02,2.141e+02,2.136e+02,2.130e+02,2.124e+02,
     a2.118e+02,2.112e+02,2.160e+02,2.222e+02,2.347e+02,2.470e+02,
     a2.593e+02,2.457e+02,2.100e+02,2.100e+02/
      data (wh5(i),i=1,34) /
     a1.200e+00,1.200e+00,9.400e-01,6.800e-01,4.100e-01,2.000e-01,
     a9.800e-02,5.400e-02,1.100e-02,8.400e-03,5.500e-03,3.800e-03,
     a2.600e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,
     a5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,
     a6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,
     a6.300e-06,1.400e-07,1.000e-09,0.000e+00/
      data (wo5(i),i=1,34) /
     a4.100e-05,4.100e-05,4.100e-05,4.300e-05,4.500e-05,4.700e-05,
     a4.900e-05,7.100e-05,9.000e-05,1.600e-04,2.400e-04,3.200e-04,
     a4.300e-04,4.700e-04,4.900e-04,5.600e-04,6.200e-04,6.200e-04,
     a6.200e-04,6.000e-04,5.600e-04,5.100e-04,4.700e-04,4.300e-04,
     a3.600e-04,3.200e-04,1.500e-04,9.200e-05,4.100e-05,1.300e-05,
     a4.300e-06,8.600e-08,4.300e-11,0.000e+00/
      do 1 i=1,34
      z(i)=z5(i)
      p(i)=p5(i)
      t(i)=t5(i)
      wh(i)=wh5(i)
      wo(i)=wo5(i)
    1 continue
      return
      end
      subroutine   tropic
      integer i
      real z1(34),p1(34),t1(34),wh1(34),wo1(34)
      real z,p,t,wh,wo
      common /sixs_atm/z(34),p(34),t(34),wh(34),wo(34)
c
c     model: tropical mc clatchey
c
      data(z1(i),i=1, 34)/
     1    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,
     2    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,
     3   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,
     4   35.,   40.,   45.,   50.,   70.,  100.,99999./
      data (p1(i),i=1,34)/
     a1.013e+03,9.040e+02,8.050e+02,7.150e+02,6.330e+02,5.590e+02,
     a4.920e+02,4.320e+02,3.780e+02,3.290e+02,2.860e+02,2.470e+02,
     a2.130e+02,1.820e+02,1.560e+02,1.320e+02,1.110e+02,9.370e+01,
     a7.890e+01,6.660e+01,5.650e+01,4.800e+01,4.090e+01,3.500e+01,
     a3.000e+01,2.570e+01,1.220e+01,6.000e+00,3.050e+00,1.590e+00,
     a8.540e-01,5.790e-02,3.000e-04,0.000e+00/
      data (t1(i),i=1,34)/
     a3.000e+02,2.940e+02,2.880e+02,2.840e+02,2.770e+02,2.700e+02,
     a2.640e+02,2.570e+02,2.500e+02,2.440e+02,2.370e+02,2.300e+02,
     a2.240e+02,2.170e+02,2.100e+02,2.040e+02,1.970e+02,1.950e+02,
     a1.990e+02,2.030e+02,2.070e+02,2.110e+02,2.150e+02,2.170e+02,
     a2.190e+02,2.210e+02,2.320e+02,2.430e+02,2.540e+02,2.650e+02,
     a2.700e+02,2.190e+02,2.100e+02,2.100e+02/
      data (wh1(i),i=1,34)/
     a1.900e+01,1.300e+01,9.300e+00,4.700e+00,2.200e+00,1.500e+00,
     a8.500e-01,4.700e-01,2.500e-01,1.200e-01,5.000e-02,1.700e-02,
     a6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,
     a5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,
     a6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,
     a6.300e-06,1.400e-07,1.000e-09,0.000e+00/
      data (wo1(i),i=1,34)/
     a5.600e-05,5.600e-05,5.400e-05,5.100e-05,4.700e-05,4.500e-05,
     a4.300e-05,4.100e-05,3.900e-05,3.900e-05,3.900e-05,4.100e-05,
     a4.300e-05,4.500e-05,4.500e-05,4.700e-05,4.700e-05,6.900e-05,
     a9.000e-05,1.400e-04,1.900e-04,2.400e-04,2.800e-04,3.200e-04,
     a3.400e-04,3.400e-04,2.400e-04,9.200e-05,4.100e-05,1.300e-05,
     a4.300e-06,8.600e-08,4.300e-11,0.000e+00/
      do 1 i=1,34
      z(i)=z1(i)
      p(i)=p1(i)
      t(i)=t1(i)
      wh(i)=wh1(i)
      wo(i)=wo1(i)
    1 continue
      return
      end
      subroutine trunca(coeff)
      real aa,x1,x2,a,x,rm,z1
      real cosang(80),weight(80),ptemp(83),pl(-1:81)
      real rmu(83),ga(83)
      integer nbmu,nang,k,j,kk,i
      real pha,betal,coeff
      common /sixs_trunc/pha(1:83),betal(0:80)
      nbmu=83
      nang=80
      do k=1,nbmu
      ptemp(k)=pha(k)
      enddo
      call gauss(-1.,1.,cosang,weight,nang)
      do 1 j=1,40
      rmu(j+1)=cosang(j)
      ga(j+1)=weight(j)
   1  continue
      rmu(1)=-1.0
      ga(1)=0.
      rmu(42)=0.
      ga(42)=0.
      do 2 j=41,80
      rmu(j+2)=cosang(j)
      ga(j+2)=weight(j)
   2  continue
      rmu(83)=1.0
      ga(83)=0.
      do 3 j=1,nbmu
      if((rmu(j).gt.0.8)) then
      go to 20
      else
      k=j-1
      endif
   3  continue
  20  continue
      do 4 j=1,nbmu
      if((rmu(j).gt.0.94)) then
      go to 21
      else
      kk=j-1
      endif
   4  continue
  21  continue
      aa=(alog10(pha(kk))-alog10(pha(k)))/
     a       (acos(rmu(kk))-acos(rmu(k)))
      x1=alog10(pha(kk))
      x2=acos(rmu(kk))
      do 5 j=kk+1,nbmu
      if(abs(rmu(j)-1.).le.1d-08) a=x1-aa*x2
      a=x1+aa*(acos(rmu(j))-x2)
      ptemp(j)=10**a
    5 continue
      do i=1,83
      pha(i)=ptemp(i)
      enddo
c
      do 10 k=0,80
      betal(k)=0.
   10 continue
      do 11 j=1,83
      x=pha(j)*ga(j)
      rm=rmu(j)
      pl(-1)=0.
      pl(0)=1.
      do 12 k=0,80
      pl(k+1)=((2*k+1.)*rm*pl(k)-k*pl(k-1))/(k+1.)
      betal(k)=betal(k)+x*pl(k)
  12  continue
  11  continue
      do 13 k=0,80
      betal(k)=(2*k+1.)*0.5*betal(k)
  13  continue
      z1=betal(0)
      coeff=1.-z1
      do k=0,80
      betal(k)=betal(k)/z1
      enddo
      return
      end
      subroutine   us62

      integer i
      real z6(34),p6(34),t6(34),wh6(34),wo6(34)
      real z,p,t,wh,wo
      common /sixs_atm/z(34),p(34),t(34),wh(34),wo(34)
c
c     model: us standard 62 mc clatchey
c
      data(z6(i),i=1, 34)/
     1    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,
     2    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,
     3   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,
     4   35.,   40.,   45.,   50.,   70.,  100.,99999./
      data (p6(i),i=1,34) /
     a1.013e+03,8.986e+02,7.950e+02,7.012e+02,6.166e+02,5.405e+02,
     a4.722e+02,4.111e+02,3.565e+02,3.080e+02,2.650e+02,2.270e+02,
     a1.940e+02,1.658e+02,1.417e+02,1.211e+02,1.035e+02,8.850e+01,
     a7.565e+01,6.467e+01,5.529e+01,4.729e+01,4.047e+01,3.467e+01,
     a2.972e+01,2.549e+01,1.197e+01,5.746e+00,2.871e+00,1.491e+00,
     a7.978e-01,5.520e-02,3.008e-04,0.000e+00/
      data (t6(i),i=1,34) /
     a2.881e+02,2.816e+02,2.751e+02,2.687e+02,2.622e+02,2.557e+02,
     a2.492e+02,2.427e+02,2.362e+02,2.297e+02,2.232e+02,2.168e+02,
     a2.166e+02,2.166e+02,2.166e+02,2.166e+02,2.166e+02,2.166e+02,
     a2.166e+02,2.166e+02,2.166e+02,2.176e+02,2.186e+02,2.196e+02,
     a2.206e+02,2.216e+02,2.265e+02,2.365e+02,2.534e+02,2.642e+02,
     a2.706e+02,2.197e+02,2.100e+02,2.100e+02/
      data (wh6(i),i=1,34) /
     a5.900e+00,4.200e+00,2.900e+00,1.800e+00,1.100e+00,6.400e-01,
     a3.800e-01,2.100e-01,1.200e-01,4.600e-02,1.800e-02,8.200e-03,
     a3.700e-03,1.800e-03,8.400e-04,7.200e-04,6.100e-04,5.200e-04,
     a4.400e-04,4.400e-04,4.400e-04,4.800e-04,5.200e-04,5.700e-04,
     a6.100e-04,6.600e-04,3.800e-04,1.600e-04,6.700e-05,3.200e-05,
     a1.200e-05,1.500e-07,1.000e-09,0.000e+00/
      data (wo6(i),i=1,34) /
     a5.400e-05,5.400e-05,5.400e-05,5.000e-05,4.600e-05,4.600e-05,
     a4.500e-05,4.900e-05,5.200e-05,7.100e-05,9.000e-05,1.300e-04,
     a1.600e-04,1.700e-04,1.900e-04,2.100e-04,2.400e-04,2.800e-04,
     a3.200e-04,3.500e-04,3.800e-04,3.800e-04,3.900e-04,3.800e-04,
     a3.600e-04,3.400e-04,2.000e-04,1.100e-04,4.900e-05,1.700e-05,
     a4.000e-06,8.600e-08,4.300e-11,0.000e+00/
      do 1 i=1,34
      z(i)=z6(i)
      p(i)=p6(i)
      t(i)=t6(i)
      wh(i)=wh6(i)
      wo(i)=wo6(i)
    1 continue
      return
      end
      subroutine varsol (jday,month,
     a                   dsol)
 
      real dsol,pi,om
      integer jday,month,j

c     calculation of the variability of the solar constant during the
c     year.
c     jday is the number of the day in the month
c     dsol is a multiplicative factor to apply to the mean value of
c     solar constant
 
      if (month.le.2) goto 1
      if (month.gt.8) goto 2
      j=31*(month-1)-((month-1)/2)-2+jday
      goto 3
    1 j=31*(month-1)+jday
      goto 3
    2 j=31*(month-1)-((month-2)/2)-2+jday
 
    3 pi=2.*acos (0.)
      om=(.9856*float(j-4))*pi/180.
      dsol=1./((1.-.01673*cos(om))**2)
      return
      end
       subroutine   wate
       integer i,j
       real phr(10,83)
       real ph
       common /sixs_aerbas/ ph(10,83)
c
c    model: water-soluble
c
            DATA ((PHR(I,J),J=1,83),I=01,01) /
     *0.4115E+00,0.4045E+00,0.3805E+00,0.3495E+00,0.3192E+00,0.2943E+00,
     *0.2768E+00,0.2659E+00,0.2592E+00,0.2538E+00,0.2479E+00,0.2411E+00,
     *0.2336E+00,0.2255E+00,0.2175E+00,0.2098E+00,0.2026E+00,0.1961E+00,
     *0.1903E+00,0.1854E+00,0.1812E+00,0.1778E+00,0.1752E+00,0.1734E+00,
     *0.1723E+00,0.1719E+00,0.1724E+00,0.1736E+00,0.1756E+00,0.1784E+00,
     *0.1820E+00,0.1866E+00,0.1920E+00,0.1985E+00,0.2061E+00,0.2149E+00,
     *0.2249E+00,0.2363E+00,0.2492E+00,0.2638E+00,0.2803E+00,0.2893E+00,
     *0.2988E+00,0.3195E+00,0.3428E+00,0.3688E+00,0.3979E+00,0.4306E+00,
     *0.4671E+00,0.5079E+00,0.5537E+00,0.6048E+00,0.6622E+00,0.7264E+00,
     *0.7985E+00,0.8794E+00,0.9701E+00,0.1072E+01,0.1186E+01,0.1315E+01,
     *0.1460E+01,0.1622E+01,0.1805E+01,0.2011E+01,0.2242E+01,0.2503E+01,
     *0.2796E+01,0.3125E+01,0.3496E+01,0.3913E+01,0.4383E+01,0.4912E+01,
     *0.5510E+01,0.6185E+01,0.6951E+01,0.7825E+01,0.8828E+01,0.9991E+01,
     *0.1136E+02,0.1297E+02,0.1491E+02,0.1711E+02,0.1834E+02/
            DATA ((PHR(I,J),J=1,83),I=02,02) /
     *0.3918E+00,0.3859E+00,0.3654E+00,0.3384E+00,0.3117E+00,0.2895E+00,
     *0.2736E+00,0.2635E+00,0.2571E+00,0.2522E+00,0.2470E+00,0.2411E+00,
     *0.2345E+00,0.2275E+00,0.2204E+00,0.2135E+00,0.2071E+00,0.2012E+00,
     *0.1959E+00,0.1914E+00,0.1875E+00,0.1844E+00,0.1820E+00,0.1804E+00,
     *0.1794E+00,0.1792E+00,0.1797E+00,0.1810E+00,0.1831E+00,0.1860E+00,
     *0.1898E+00,0.1945E+00,0.2001E+00,0.2068E+00,0.2146E+00,0.2236E+00,
     *0.2339E+00,0.2456E+00,0.2589E+00,0.2739E+00,0.2909E+00,0.3001E+00,
     *0.3099E+00,0.3312E+00,0.3552E+00,0.3820E+00,0.4119E+00,0.4455E+00,
     *0.4830E+00,0.5249E+00,0.5718E+00,0.6243E+00,0.6829E+00,0.7486E+00,
     *0.8221E+00,0.9045E+00,0.9968E+00,0.1100E+01,0.1216E+01,0.1346E+01,
     *0.1492E+01,0.1655E+01,0.1839E+01,0.2045E+01,0.2275E+01,0.2534E+01,
     *0.2824E+01,0.3149E+01,0.3513E+01,0.3920E+01,0.4375E+01,0.4884E+01,
     *0.5454E+01,0.6092E+01,0.6807E+01,0.7611E+01,0.8516E+01,0.9543E+01,
     *0.1071E+02,0.1205E+02,0.1357E+02,0.1518E+02,0.1599E+02/
            DATA ((PHR(I,J),J=1,83),I=03,03) /
     *0.3872E+00,0.3816E+00,0.3620E+00,0.3360E+00,0.3102E+00,0.2887E+00,
     *0.2732E+00,0.2633E+00,0.2571E+00,0.2522E+00,0.2471E+00,0.2414E+00,
     *0.2350E+00,0.2283E+00,0.2214E+00,0.2148E+00,0.2085E+00,0.2028E+00,
     *0.1976E+00,0.1932E+00,0.1894E+00,0.1864E+00,0.1840E+00,0.1824E+00,
     *0.1815E+00,0.1813E+00,0.1819E+00,0.1832E+00,0.1853E+00,0.1883E+00,
     *0.1920E+00,0.1968E+00,0.2024E+00,0.2092E+00,0.2170E+00,0.2261E+00,
     *0.2364E+00,0.2483E+00,0.2617E+00,0.2768E+00,0.2939E+00,0.3032E+00,
     *0.3131E+00,0.3346E+00,0.3587E+00,0.3857E+00,0.4159E+00,0.4497E+00,
     *0.4875E+00,0.5297E+00,0.5769E+00,0.6297E+00,0.6887E+00,0.7547E+00,
     *0.8286E+00,0.9114E+00,0.1004E+01,0.1108E+01,0.1224E+01,0.1354E+01,
     *0.1500E+01,0.1664E+01,0.1847E+01,0.2053E+01,0.2284E+01,0.2542E+01,
     *0.2831E+01,0.3154E+01,0.3515E+01,0.3919E+01,0.4370E+01,0.4874E+01,
     *0.5436E+01,0.6064E+01,0.6765E+01,0.7549E+01,0.8430E+01,0.9422E+01,
     *0.1054E+02,0.1182E+02,0.1324E+02,0.1472E+02,0.1544E+02/
            DATA ((PHR(I,J),J=1,83),I=04,04) /
     *0.3737E+00,0.3687E+00,0.3509E+00,0.3269E+00,0.3030E+00,0.2830E+00,
     *0.2686E+00,0.2593E+00,0.2535E+00,0.2490E+00,0.2444E+00,0.2393E+00,
     *0.2335E+00,0.2273E+00,0.2210E+00,0.2148E+00,0.2089E+00,0.2036E+00,
     *0.1987E+00,0.1945E+00,0.1910E+00,0.1881E+00,0.1859E+00,0.1844E+00,
     *0.1836E+00,0.1835E+00,0.1842E+00,0.1855E+00,0.1877E+00,0.1907E+00,
     *0.1945E+00,0.1993E+00,0.2051E+00,0.2118E+00,0.2198E+00,0.2289E+00,
     *0.2394E+00,0.2513E+00,0.2649E+00,0.2802E+00,0.2974E+00,0.3068E+00,
     *0.3168E+00,0.3385E+00,0.3628E+00,0.3901E+00,0.4206E+00,0.4547E+00,
     *0.4928E+00,0.5353E+00,0.5829E+00,0.6361E+00,0.6955E+00,0.7620E+00,
     *0.8363E+00,0.9195E+00,0.1013E+01,0.1117E+01,0.1233E+01,0.1364E+01,
     *0.1510E+01,0.1674E+01,0.1858E+01,0.2063E+01,0.2293E+01,0.2550E+01,
     *0.2838E+01,0.3160E+01,0.3518E+01,0.3919E+01,0.4365E+01,0.4863E+01,
     *0.5416E+01,0.6033E+01,0.6719E+01,0.7483E+01,0.8337E+01,0.9292E+01,
     *0.1036E+02,0.1156E+02,0.1289E+02,0.1423E+02,0.1486E+02/
            DATA ((PHR(I,J),J=1,83),I=05,05) /
     *0.3651E+00,0.3607E+00,0.3449E+00,0.3233E+00,0.3016E+00,0.2832E+00,
     *0.2697E+00,0.2609E+00,0.2552E+00,0.2509E+00,0.2465E+00,0.2418E+00,
     *0.2364E+00,0.2307E+00,0.2249E+00,0.2191E+00,0.2137E+00,0.2086E+00,
     *0.2041E+00,0.2001E+00,0.1968E+00,0.1940E+00,0.1919E+00,0.1905E+00,
     *0.1898E+00,0.1897E+00,0.1904E+00,0.1919E+00,0.1941E+00,0.1971E+00,
     *0.2011E+00,0.2059E+00,0.2118E+00,0.2187E+00,0.2267E+00,0.2361E+00,
     *0.2467E+00,0.2589E+00,0.2727E+00,0.2883E+00,0.3059E+00,0.3155E+00,
     *0.3257E+00,0.3478E+00,0.3726E+00,0.4004E+00,0.4315E+00,0.4662E+00,
     *0.5050E+00,0.5483E+00,0.5967E+00,0.6507E+00,0.7110E+00,0.7783E+00,
     *0.8536E+00,0.9376E+00,0.1032E+01,0.1137E+01,0.1254E+01,0.1385E+01,
     *0.1531E+01,0.1695E+01,0.1878E+01,0.2083E+01,0.2311E+01,0.2566E+01,
     *0.2850E+01,0.3166E+01,0.3518E+01,0.3910E+01,0.4344E+01,0.4825E+01,
     *0.5358E+01,0.5947E+01,0.6597E+01,0.7314E+01,0.8106E+01,0.8978E+01,
     *0.9939E+01,0.1099E+02,0.1211E+02,0.1319E+02,0.1367E+02/
            DATA ((PHR(I,J),J=1,83),I=06,06) /
     *0.3540E+00,0.3501E+00,0.3360E+00,0.3166E+00,0.2969E+00,0.2801E+00,
     *0.2677E+00,0.2594E+00,0.2541E+00,0.2500E+00,0.2461E+00,0.2417E+00,
     *0.2369E+00,0.2317E+00,0.2263E+00,0.2211E+00,0.2160E+00,0.2113E+00,
     *0.2070E+00,0.2033E+00,0.2001E+00,0.1976E+00,0.1956E+00,0.1943E+00,
     *0.1937E+00,0.1937E+00,0.1945E+00,0.1960E+00,0.1982E+00,0.2013E+00,
     *0.2053E+00,0.2102E+00,0.2162E+00,0.2232E+00,0.2313E+00,0.2408E+00,
     *0.2516E+00,0.2639E+00,0.2779E+00,0.2937E+00,0.3115E+00,0.3213E+00,
     *0.3315E+00,0.3540E+00,0.3791E+00,0.4073E+00,0.4387E+00,0.4739E+00,
     *0.5131E+00,0.5569E+00,0.6057E+00,0.6603E+00,0.7211E+00,0.7890E+00,
     *0.8647E+00,0.9493E+00,0.1044E+01,0.1149E+01,0.1267E+01,0.1398E+01,
     *0.1545E+01,0.1708E+01,0.1891E+01,0.2095E+01,0.2322E+01,0.2575E+01,
     *0.2856E+01,0.3169E+01,0.3517E+01,0.3902E+01,0.4328E+01,0.4799E+01,
     *0.5318E+01,0.5890E+01,0.6519E+01,0.7208E+01,0.7963E+01,0.8788E+01,
     *0.9685E+01,0.1065E+02,0.1166E+02,0.1261E+02,0.1301E+02/
            DATA ((PHR(I,J),J=1,83),I=07,07) /
     *0.3121E+00,0.3097E+00,0.3008E+00,0.2882E+00,0.2753E+00,0.2643E+00,
     *0.2562E+00,0.2509E+00,0.2473E+00,0.2445E+00,0.2417E+00,0.2384E+00,
     *0.2348E+00,0.2307E+00,0.2265E+00,0.2223E+00,0.2182E+00,0.2144E+00,
     *0.2109E+00,0.2078E+00,0.2052E+00,0.2030E+00,0.2014E+00,0.2004E+00,
     *0.2000E+00,0.2002E+00,0.2011E+00,0.2027E+00,0.2051E+00,0.2082E+00,
     *0.2123E+00,0.2173E+00,0.2232E+00,0.2303E+00,0.2386E+00,0.2482E+00,
     *0.2591E+00,0.2717E+00,0.2859E+00,0.3019E+00,0.3201E+00,0.3300E+00,
     *0.3404E+00,0.3633E+00,0.3889E+00,0.4176E+00,0.4496E+00,0.4854E+00,
     *0.5253E+00,0.5699E+00,0.6196E+00,0.6749E+00,0.7367E+00,0.8055E+00,
     *0.8822E+00,0.9677E+00,0.1063E+01,0.1169E+01,0.1288E+01,0.1419E+01,
     *0.1566E+01,0.1730E+01,0.1912E+01,0.2115E+01,0.2341E+01,0.2591E+01,
     *0.2869E+01,0.3177E+01,0.3518E+01,0.3895E+01,0.4309E+01,0.4765E+01,
     *0.5265E+01,0.5811E+01,0.6405E+01,0.7049E+01,0.7744E+01,0.8489E+01,
     *0.9280E+01,0.1010E+02,0.1093E+02,0.1165E+02,0.1192E+02/
            DATA ((PHR(I,J),J=1,83),I=08,08) /
     *0.3070E+00,0.3061E+00,0.3027E+00,0.2975E+00,0.2918E+00,0.2865E+00,
     *0.2821E+00,0.2787E+00,0.2760E+00,0.2735E+00,0.2711E+00,0.2684E+00,
     *0.2656E+00,0.2626E+00,0.2594E+00,0.2562E+00,0.2530E+00,0.2500E+00,
     *0.2471E+00,0.2446E+00,0.2423E+00,0.2404E+00,0.2390E+00,0.2380E+00,
     *0.2375E+00,0.2377E+00,0.2385E+00,0.2400E+00,0.2422E+00,0.2453E+00,
     *0.2493E+00,0.2543E+00,0.2604E+00,0.2677E+00,0.2762E+00,0.2861E+00,
     *0.2976E+00,0.3108E+00,0.3258E+00,0.3428E+00,0.3620E+00,0.3725E+00,
     *0.3836E+00,0.4079E+00,0.4351E+00,0.4655E+00,0.4993E+00,0.5371E+00,
     *0.5791E+00,0.6258E+00,0.6776E+00,0.7351E+00,0.7988E+00,0.8694E+00,
     *0.9476E+00,0.1034E+01,0.1130E+01,0.1236E+01,0.1353E+01,0.1482E+01,
     *0.1625E+01,0.1783E+01,0.1957E+01,0.2148E+01,0.2359E+01,0.2590E+01,
     *0.2844E+01,0.3121E+01,0.3424E+01,0.3754E+01,0.4112E+01,0.4498E+01,
     *0.4913E+01,0.5356E+01,0.5826E+01,0.6320E+01,0.6833E+01,0.7358E+01,
     *0.7884E+01,0.8390E+01,0.8846E+01,0.9187E+01,0.9295E+01/
            DATA ((PHR(I,J),J=1,83),I=09,09) /
     *0.3321E+00,0.3315E+00,0.3294E+00,0.3266E+00,0.3238E+00,0.3214E+00,
     *0.3192E+00,0.3169E+00,0.3142E+00,0.3111E+00,0.3075E+00,0.3036E+00,
     *0.2994E+00,0.2950E+00,0.2905E+00,0.2860E+00,0.2817E+00,0.2775E+00,
     *0.2735E+00,0.2698E+00,0.2665E+00,0.2635E+00,0.2609E+00,0.2587E+00,
     *0.2571E+00,0.2561E+00,0.2556E+00,0.2558E+00,0.2568E+00,0.2586E+00,
     *0.2613E+00,0.2650E+00,0.2697E+00,0.2756E+00,0.2827E+00,0.2913E+00,
     *0.3013E+00,0.3131E+00,0.3267E+00,0.3422E+00,0.3600E+00,0.3698E+00,
     *0.3802E+00,0.4030E+00,0.4287E+00,0.4575E+00,0.4899E+00,0.5261E+00,
     *0.5665E+00,0.6115E+00,0.6617E+00,0.7175E+00,0.7795E+00,0.8484E+00,
     *0.9248E+00,0.1010E+01,0.1103E+01,0.1208E+01,0.1323E+01,0.1451E+01,
     *0.1592E+01,0.1749E+01,0.1922E+01,0.2113E+01,0.2324E+01,0.2557E+01,
     *0.2813E+01,0.3095E+01,0.3403E+01,0.3740E+01,0.4106E+01,0.4502E+01,
     *0.4928E+01,0.5383E+01,0.5863E+01,0.6364E+01,0.6878E+01,0.7395E+01,
     *0.7898E+01,0.8366E+01,0.8764E+01,0.9041E+01,0.9119E+01/
            DATA ((PHR(I,J),J=1,83),I=10,10) /
     *0.4248E+00,0.4242E+00,0.4221E+00,0.4189E+00,0.4153E+00,0.4116E+00,
     *0.4081E+00,0.4045E+00,0.4006E+00,0.3964E+00,0.3918E+00,0.3869E+00,
     *0.3818E+00,0.3764E+00,0.3709E+00,0.3654E+00,0.3600E+00,0.3547E+00,
     *0.3495E+00,0.3446E+00,0.3401E+00,0.3359E+00,0.3321E+00,0.3288E+00,
     *0.3260E+00,0.3239E+00,0.3224E+00,0.3218E+00,0.3219E+00,0.3230E+00,
     *0.3251E+00,0.3282E+00,0.3326E+00,0.3383E+00,0.3455E+00,0.3542E+00,
     *0.3646E+00,0.3768E+00,0.3911E+00,0.4075E+00,0.4263E+00,0.4366E+00,
     *0.4476E+00,0.4717E+00,0.4989E+00,0.5293E+00,0.5633E+00,0.6011E+00,
     *0.6431E+00,0.6896E+00,0.7410E+00,0.7977E+00,0.8603E+00,0.9291E+00,
     *0.1005E+01,0.1088E+01,0.1179E+01,0.1278E+01,0.1387E+01,0.1506E+01,
     *0.1636E+01,0.1778E+01,0.1933E+01,0.2100E+01,0.2283E+01,0.2480E+01,
     *0.2693E+01,0.2923E+01,0.3169E+01,0.3433E+01,0.3713E+01,0.4009E+01,
     *0.4319E+01,0.4642E+01,0.4973E+01,0.5308E+01,0.5640E+01,0.5962E+01,
     *0.6262E+01,0.6528E+01,0.6740E+01,0.6876E+01,0.6911E+01/
c
       do 1 i=1,10
       do 1 j=1,83
       ph(i,j)=phr(i,j)
  1    continue
       return
       end
