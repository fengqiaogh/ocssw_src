	subroutine extem101_64(tt,slat)
c * Extrapolate temperature profile up to 0.005 mb
c .... version of 24.04.02

c   Input:
c      tt   = top-down 101-level temperature profile with -1. at missing upper levels
c	slat = latitude (deg,+N,-S)
c   Output:
c      tt = same profile with missing levels filled in

        implicit real*8 (A-H,O-Z)

	parameter (nx=9,ny=35,nz=3)
	dimension tt(*),tx(nx),lx(nx),cc(nz,0:nx,ny),plat(nz),cfl(0:nx)
	common/extcoeff_64/coef(0:nx,ny,nz)
	data init/1/,lx/36,40,45,51,56,64,70,76,86/
	data x1/75./,x2/45./,x3/15./

        save

	if(init.ne.0) then
	   do j=1,ny
	      do i=0,nx
	         y1=coef(i,j,1)
	         y2=coef(i,j,2)
	         y3=coef(i,j,3)
	         call cofit3_64(x1,x2,x3,y1,y2,y3,c1,c2,c3)
	         cc(1,i,j)=c1
	         cc(2,i,j)=c2
	         cc(3,i,j)=c3
	      enddo
	   enddo
	   init=0
	endif

	alat=abs(slat)
	xx=1.
	do k=1,nz
	   plat(k)=xx
	   xx=xx*alat
	enddo

c * load the predictor array
	do i=1,nx
	   l=lx(i)
	   tx(i)=tt(l)
	enddo

	do j=1,ny
	   if(tt(j).gt.0.) return
	   do i=0,nx
	      cy=0.
	      do k=1,nz
	         cy=cy+cc(k,i,j)*plat(k)
	      enddo
	      cfl(i)=cy
	   enddo
	   sum=cfl(0)
	   do i=1,nx
	      sum=sum+cfl(i)*tx(i)
	   enddo
	   tt(j)=sum
	enddo

	return
	end

	block data extrap_coeffs_64
	parameter (nx=9,ny=35,nz=3)
        real*8 coef(0:nx,ny,nz)
	common/extcoeff_64/coef

c * Temperature-Extrapolation Coefficients ...
c	derived from NESDIS PROF1200                                

c +++ Zone 1

c *   0.0050 mb:
      data (coef(i, 1,1),i=0,nx)/
     +  768.137268,   -0.294463,   -1.938120,    1.445014,    1.010484,
     +   -1.501080,   -0.333852,   -0.460916,   -0.289109,   -0.185883/
c *   0.0161 mb:
      data (coef(i, 2,1),i=0,nx)/
     +  729.598755,   -0.048072,   -2.192048,    1.556883,    1.041545,
     +   -1.527906,   -0.179924,   -0.397472,   -0.527085,   -0.022882/
c *   0.0384 mb:
      data (coef(i, 3,1),i=0,nx)/
     +  580.007996,    0.256773,   -1.953103,    1.262224,    0.855626,
     +   -1.226672,    0.021504,   -0.272799,   -0.672050,    0.167226/
c *   0.0769 mb:
      data (coef(i, 4,1),i=0,nx)/
     +  378.418610,    0.481383,   -1.547368,    0.949718,    0.557970,
     +   -0.797499,    0.236666,   -0.090441,   -0.724859,    0.324352/
c *   0.1370 mb:
      data (coef(i, 5,1),i=0,nx)/
     +  227.726013,    0.399807,   -0.767517,    0.339597,    0.334890,
     +   -0.304244,    0.242581,    0.200534,   -0.710240,    0.359780/
c *   0.2244 mb:
      data (coef(i, 6,1),i=0,nx)/
     +  131.540237,    0.116343,    0.181285,   -0.404714,    0.194162,
     +    0.168128,    0.112673,    0.505991,   -0.649249,    0.319506/
c *   0.3454 mb:
      data (coef(i, 7,1),i=0,nx)/
     +   36.948456,   -0.219132,    1.068293,   -1.006580,    0.041164,
     +    0.596757,    0.016176,    0.724429,   -0.524707,    0.280223/
c *   0.5064 mb:
      data (coef(i, 8,1),i=0,nx)/
     +  -47.064129,   -0.502069,    1.832718,   -1.521932,   -0.102700,
     +    0.973025,   -0.065917,    0.916755,   -0.417642,    0.247873/
c *   0.7140 mb:
      data (coef(i, 9,1),i=0,nx)/
     + -125.058113,   -0.371252,    1.929798,   -1.497271,   -0.445293,
     +    1.206479,   -0.047493,    1.051199,   -0.411198,    0.284401/
c *   0.9753 mb:
      data (coef(i,10,1),i=0,nx)/
     + -195.826279,   -0.252554,    2.017885,   -1.474895,   -0.756146,
     +    1.418304,   -0.030777,    1.173188,   -0.405352,    0.317545/
c *   1.2972 mb:
      data (coef(i,11,1),i=0,nx)/
     + -211.840836,   -0.087453,    1.892584,   -1.498378,   -0.888045,
     +    1.557249,   -0.013888,    1.153322,   -0.274422,    0.228628/
c *   1.6872 mb:
      data (coef(i,12,1),i=0,nx)/
     + -222.487411,    0.077679,    1.743851,   -1.517234,   -0.992018,
     +    1.676924,    0.003067,    1.119941,   -0.144600,    0.140581/
c *   2.1526 mb:
      data (coef(i,13,1),i=0,nx)/
     + -220.345123,    0.284957,    1.478203,   -1.409153,   -1.098716,
     +    1.732799,    0.005315,    1.062679,   -0.051710,    0.083439/
c *   2.7009 mb:
      data (coef(i,14,1),i=0,nx)/
     + -191.945847,    0.575167,    0.992410,   -1.054349,   -1.230327,
     +    1.674808,   -0.024700,    0.961055,   -0.019512,    0.071966/
c *   3.3398 mb:
      data (coef(i,15,1),i=0,nx)/
     + -163.447739,    0.789576,    0.662758,   -0.848763,   -1.165136,
     +    1.486376,   -0.005209,    0.856173,    0.036504,    0.000516/
c *   4.0770 mb:
      data (coef(i,16,1),i=0,nx)/
     + -135.140701,    0.938526,    0.466544,   -0.770585,   -0.932123,
     +    1.188257,    0.055934,    0.750415,    0.112469,   -0.122608/
c *   4.9204 mb:
      data (coef(i,17,1),i=0,nx)/
     + -110.534798,    1.079643,    0.269122,   -0.685459,   -0.724919,
     +    0.927309,    0.105440,    0.666502,    0.180002,   -0.240317/
c *   5.8776 mb:
      data (coef(i,18,1),i=0,nx)/
     +  -86.496407,    1.197924,    0.180591,   -0.680802,   -0.560434,
     +    0.725634,    0.114811,    0.524036,    0.254776,   -0.292023/
c *   6.9567 mb:
      data (coef(i,19,1),i=0,nx)/
     +  -63.631668,    1.308646,    0.105890,   -0.683525,   -0.407437,
     +    0.538658,    0.120171,    0.383015,    0.326700,   -0.335437/
c *   8.1655 mb:
      data (coef(i,20,1),i=0,nx)/
     +  -53.883522,    1.446655,   -0.084576,   -0.584522,   -0.321630,
     +    0.418490,    0.138251,    0.323507,    0.266197,   -0.297677/
c *   9.5119 mb:
      data (coef(i,21,1),i=0,nx)/
     +  -45.056870,    1.579389,   -0.270613,   -0.486306,   -0.242174,
     +    0.306220,    0.155973,    0.269677,    0.203612,   -0.258671/
c *  11.0038 mb:
      data (coef(i,22,1),i=0,nx)/
     +  -37.774498,    1.605633,   -0.336989,   -0.423904,   -0.179642,
     +    0.233162,    0.156346,    0.226620,    0.162179,   -0.225690/
c *  12.6492 mb:
      data (coef(i,23,1),i=0,nx)/
     +  -31.381693,    1.580467,   -0.344828,   -0.379910,   -0.126498,
     +    0.180357,    0.148424,    0.189608,    0.131713,   -0.196274/
c *  14.4559 mb:
      data (coef(i,24,1),i=0,nx)/
     +  -25.256893,    1.556356,   -0.352337,   -0.337762,   -0.075582,
     +    0.129767,    0.140834,    0.154147,    0.102526,   -0.168092/
c *  16.4318 mb:
      data (coef(i,25,1),i=0,nx)/
     +  -19.332832,    1.532327,   -0.358523,   -0.296751,   -0.023514,
     +    0.076359,    0.135129,    0.118184,    0.075598,   -0.141048/
c *  18.5847 mb:
      data (coef(i,26,1),i=0,nx)/
     +  -13.621389,    1.508887,   -0.364070,   -0.257120,    0.027774,
     +    0.023140,    0.130262,    0.082869,    0.050141,   -0.115058/
c *  20.9224 mb:
      data (coef(i,27,1),i=0,nx)/
     +  -10.552977,    1.475920,   -0.355611,   -0.209939,    0.047681,
     +   -0.006736,    0.119757,    0.072153,    0.024709,   -0.098524/
c *  23.4526 mb:
      data (coef(i,28,1),i=0,nx)/
     +  -11.404002,    1.427825,   -0.325822,   -0.150293,    0.020678,
     +   -0.002060,    0.100506,    0.098321,   -0.001260,   -0.095886/
c *  26.1829 mb:
      data (coef(i,29,1),i=0,nx)/
     +  -12.431881,    1.380402,   -0.295441,   -0.093749,   -0.004390,
     +    0.002642,    0.081415,    0.123695,   -0.026677,   -0.092560/
c *  29.1210 mb:
      data (coef(i,30,1),i=0,nx)/
     +  -13.700622,    1.333231,   -0.263907,   -0.040468,   -0.027290,
     +    0.007437,    0.062283,    0.148370,   -0.051712,   -0.088305/
c *  32.2744 mb:
      data (coef(i,31,1),i=0,nx)/
     +  -12.049180,    1.288188,   -0.233483,   -0.020884,   -0.028961,
     +    0.007528,    0.048806,    0.133084,   -0.050300,   -0.074679/
c *  35.6505 mb:
      data (coef(i,32,1),i=0,nx)/
     +   -9.317835,    1.244821,   -0.204065,   -0.014504,   -0.022521,
     +    0.005827,    0.037740,    0.102881,   -0.038848,   -0.057748/
c *  39.2566 mb:
      data (coef(i,33,1),i=0,nx)/
     +   -6.672400,    1.202820,   -0.175573,   -0.008325,   -0.016283,
     +    0.004179,    0.027023,    0.073628,   -0.027756,   -0.041350/
c *  43.1001 mb:
      data (coef(i,34,1),i=0,nx)/
     +   -4.107985,    1.162105,   -0.147954,   -0.002335,   -0.010237,
     +    0.002583,    0.016633,    0.045271,   -0.017004,   -0.025454/
c *  47.1882 mb:
      data (coef(i,35,1),i=0,nx)/
     +   -1.620239,    1.122606,   -0.121160,    0.003476,   -0.004369,
     +    0.001033,    0.006554,    0.017761,   -0.006572,   -0.010032/

c +++ Zone 2

c *   0.0050 mb:
      data (coef(i, 1,2),i=0,nx)/
     +  834.861145,   -0.395518,   -0.168292,   -0.522989,    0.135368,
     +   -0.648098,   -0.120294,   -0.432757,   -0.616573,   -0.049551/
c *   0.0161 mb:
      data (coef(i, 2,2),i=0,nx)/
     +  706.697998,   -0.564473,    0.061239,   -0.490616,    0.279786,
     +   -0.515674,   -0.074882,   -0.145787,   -0.887622,    0.150642/
c *   0.0384 mb:
      data (coef(i, 3,2),i=0,nx)/
     +  456.640137,   -0.515665,    0.238243,   -0.325252,    0.306180,
     +   -0.235547,   -0.028916,    0.203108,   -0.970324,    0.303589/
c *   0.0769 mb:
      data (coef(i, 4,2),i=0,nx)/
     +  196.650543,   -0.580368,    0.469164,   -0.167309,    0.389823,
     +    0.028119,    0.028069,    0.512286,   -0.967577,    0.455387/
c *   0.1370 mb:
      data (coef(i, 5,2),i=0,nx)/
     +   55.473385,   -0.577106,    0.510232,   -0.037092,    0.481069,
     +    0.095268,    0.059780,    0.602277,   -0.792924,    0.480604/
c *   0.2244 mb:
      data (coef(i, 6,2),i=0,nx)/
     +    3.620781,   -0.523341,    0.420839,    0.076334,    0.556607,
     +    0.033327,    0.063618,    0.596428,   -0.589134,    0.438588/
c *   0.3454 mb:
      data (coef(i, 7,2),i=0,nx)/
     +  -45.342590,   -0.570448,    0.387072,    0.196123,    0.643398,
     +   -0.025933,    0.053534,    0.719213,   -0.558915,    0.459968/
c *   0.5064 mb:
      data (coef(i, 8,2),i=0,nx)/
     +  -88.503006,   -0.608382,    0.353316,    0.301524,    0.715639,
     +   -0.073664,    0.042024,    0.832244,   -0.535806,    0.479966/
c *   0.7140 mb:
      data (coef(i, 9,2),i=0,nx)/
     + -120.005302,   -0.541255,    0.223263,    0.373417,    0.655807,
     +    0.010448,   -0.035588,    1.041446,   -0.612080,    0.524942/
c *   0.9753 mb:
      data (coef(i,10,2),i=0,nx)/
     + -148.589035,   -0.480346,    0.105258,    0.438649,    0.601519,
     +    0.086768,   -0.106009,    1.231265,   -0.681287,    0.565751/
c *   1.2972 mb:
      data (coef(i,11,2),i=0,nx)/
     + -161.973740,   -0.371849,    0.008552,    0.486533,    0.557450,
     +    0.112235,   -0.175632,    1.253408,   -0.627670,    0.561722/
c *   1.6872 mb:
      data (coef(i,12,2),i=0,nx)/
     + -173.125641,   -0.266075,   -0.083040,    0.532244,    0.515987,
     +    0.132649,   -0.241772,    1.260636,   -0.568181,    0.555640/
c *   2.1526 mb:
      data (coef(i,13,2),i=0,nx)/
     + -173.081696,   -0.094634,   -0.220959,    0.567489,    0.477990,
     +    0.120870,   -0.283719,    1.244042,   -0.533282,    0.546141/
c *   2.7009 mb:
      data (coef(i,14,2),i=0,nx)/
     + -150.789734,    0.220479,   -0.455112,    0.578471,    0.446871,
     +    0.041562,   -0.277308,    1.177748,   -0.543711,    0.525773/
c *   3.3398 mb:
      data (coef(i,15,2),i=0,nx)/
     + -130.673721,    0.435536,   -0.553035,    0.536432,    0.402999,
     +   -0.011230,   -0.249925,    1.061322,   -0.498618,    0.491557/
c *   4.0770 mb:
      data (coef(i,16,2),i=0,nx)/
     + -112.541351,    0.564179,   -0.533735,    0.448692,    0.348486,
     +   -0.040873,   -0.204408,    0.902095,   -0.405841,    0.445350/
c *   4.9204 mb:
      data (coef(i,17,2),i=0,nx)/
     +  -96.168854,    0.685395,   -0.516674,    0.364636,    0.299450,
     +   -0.066547,   -0.160329,    0.753376,   -0.318513,    0.400629/
c *   5.8776 mb:
      data (coef(i,18,2),i=0,nx)/
     + -121.107452,    0.789144,   -0.466781,    0.308981,    0.265531,
     +   -0.041420,   -0.096691,    0.627374,   -0.172447,    0.330783/
c *   6.9567 mb:
      data (coef(i,19,2),i=0,nx)/
     + -148.557831,    0.886490,   -0.416297,    0.258454,    0.234542,
     +   -0.012946,   -0.034285,    0.509281,   -0.027978,    0.261965/
c *   8.1655 mb:
      data (coef(i,20,2),i=0,nx)/
     + -146.946228,    0.924881,   -0.351371,    0.173360,    0.243558,
     +   -0.023024,   -0.043987,    0.511167,   -0.028737,    0.240058/
c *   9.5119 mb:
      data (coef(i,21,2),i=0,nx)/
     + -144.347504,    0.959377,   -0.288869,    0.090871,    0.253623,
     +   -0.034051,   -0.055879,    0.517345,   -0.034760,    0.220859/
c *  11.0038 mb:
      data (coef(i,22,2),i=0,nx)/
     + -136.725281,    0.967216,   -0.250049,    0.038883,    0.268749,
     +   -0.045579,   -0.034804,    0.465370,   -0.010830,    0.193404/
c *  12.6492 mb:
      data (coef(i,23,2),i=0,nx)/
     + -126.862816,    0.962160,   -0.223352,    0.002548,    0.285976,
     +   -0.057106,    0.001579,    0.386704,    0.026907,    0.162579/
c *  14.4559 mb:
      data (coef(i,24,2),i=0,nx)/
     + -117.413834,    0.957316,   -0.197774,   -0.032263,    0.302480,
     +   -0.068150,    0.036436,    0.311337,    0.063062,    0.133046/
c *  16.4318 mb:
      data (coef(i,25,2),i=0,nx)/
     + -108.387199,    0.951026,   -0.170067,   -0.067915,    0.319109,
     +   -0.079884,    0.071843,    0.243068,    0.098536,    0.099678/
c *  18.5847 mb:
      data (coef(i,26,2),i=0,nx)/
     +  -99.728241,    0.944342,   -0.142208,   -0.103052,    0.335396,
     +   -0.091603,    0.106632,    0.179040,    0.132930,    0.065653/
c *  20.9224 mb:
      data (coef(i,27,2),i=0,nx)/
     +  -90.187675,    0.953354,   -0.131290,   -0.116374,    0.326715,
     +   -0.094256,    0.128964,    0.132450,    0.141525,    0.044614/
c *  23.4526 mb:
      data (coef(i,28,2),i=0,nx)/
     +  -79.100502,    0.986260,   -0.145695,   -0.097073,    0.280154,
     +   -0.083285,    0.132998,    0.111127,    0.111379,    0.042697/
c *  26.1829 mb:
      data (coef(i,29,2),i=0,nx)/
     +  -68.699837,    1.017853,   -0.157891,   -0.083327,    0.241075,
     +   -0.074152,    0.137479,    0.089334,    0.082811,    0.041165/
c *  29.1210 mb:
      data (coef(i,30,2),i=0,nx)/
     +  -59.048897,    1.048161,   -0.167399,   -0.076558,    0.211126,
     +   -0.067266,    0.142594,    0.066655,    0.055907,    0.040111/
c *  32.2744 mb:
      data (coef(i,31,2),i=0,nx)/
     +  -48.262180,    1.062744,   -0.160901,   -0.062914,    0.173613,
     +   -0.056035,    0.123367,    0.051704,    0.041447,    0.034117/
c *  35.6505 mb:
      data (coef(i,32,2),i=0,nx)/
     +  -37.250641,    1.071060,   -0.148434,   -0.046917,    0.133942,
     +   -0.043364,    0.095243,    0.039982,    0.032004,    0.026359/
c *  39.2566 mb:
      data (coef(i,33,2),i=0,nx)/
     +  -26.585524,    1.079114,   -0.136360,   -0.031422,    0.095519,
     +   -0.031093,    0.068005,    0.028629,    0.022858,    0.018845/
c *  43.1001 mb:
      data (coef(i,34,2),i=0,nx)/
     +  -16.247211,    1.086921,   -0.124655,   -0.016402,    0.058273,
     +   -0.019197,    0.041600,    0.017624,    0.013991,    0.011561/
c *  47.1882 mb:
      data (coef(i,35,2),i=0,nx)/
     +   -6.217585,    1.094496,   -0.113301,   -0.001831,    0.022139,
     +   -0.007657,    0.015985,    0.006947,    0.005390,    0.004494/

c +++ Zone 3

c *   0.0050 mb:
      data (coef(i, 1,3),i=0,nx)/
     +  253.099335,   -0.059435,   -0.113692,   -0.076584,    0.018085,
     +    0.109696,   -0.245047,    0.207798,   -0.158821,   -0.005658/
c *   0.0161 mb:
      data (coef(i, 2,3),i=0,nx)/
     +  393.724976,   -0.159516,   -0.305137,   -0.205541,    0.048536,
     +    0.294412,   -0.657678,    0.557706,   -0.426257,   -0.015186/
c *   0.0384 mb:
      data (coef(i, 3,3),i=0,nx)/
     +  488.316467,   -0.220589,   -0.421964,   -0.284236,    0.067118,
     +    0.407132,   -0.909481,    0.771234,   -0.589456,   -0.021001/
c *   0.0769 mb:
      data (coef(i, 4,3),i=0,nx)/
     +  571.788940,   -0.274628,   -0.525335,   -0.353867,    0.083561,
     +    0.506871,   -1.132283,    0.960168,   -0.733859,   -0.026146/
c *   0.1370 mb:
      data (coef(i, 5,3),i=0,nx)/
     +  564.131775,   -0.280373,   -0.502551,   -0.299793,    0.065189,
     +    0.444018,   -1.027926,    0.888957,   -0.690901,   -0.011248/
c *   0.2244 mb:
      data (coef(i, 6,3),i=0,nx)/
     +  496.362030,   -0.249809,   -0.394284,   -0.179865,    0.027459,
     +    0.290057,   -0.717002,    0.634090,   -0.503432,    0.004950/
c *   0.3454 mb:
      data (coef(i, 7,3),i=0,nx)/
     +  432.359314,   -0.211118,   -0.286078,   -0.105599,   -0.000052,
     +    0.184573,   -0.442425,    0.356257,   -0.266806,   -0.011249/
c *   0.5064 mb:
      data (coef(i, 8,3),i=0,nx)/
     +  375.968506,   -0.173903,   -0.194453,   -0.041546,   -0.020229,
     +    0.090298,   -0.207287,    0.122826,   -0.062737,   -0.026678/
c *   0.7140 mb:
      data (coef(i, 9,3),i=0,nx)/
     +  335.498596,   -0.064613,   -0.226756,   -0.032027,    0.072557,
     +   -0.012652,   -0.217699,    0.255347,   -0.032882,   -0.068283/
c *   0.9753 mb:
      data (coef(i,10,3),i=0,nx)/
     +  298.778290,    0.034553,   -0.256067,   -0.023391,    0.156749,
     +   -0.106065,   -0.227148,    0.375591,   -0.005793,   -0.106034/
c *   1.2972 mb:
      data (coef(i,11,3),i=0,nx)/
     +  320.508484,    0.055737,   -0.215899,   -0.126142,    0.126488,
     +   -0.119170,   -0.163850,    0.312299,    0.057370,   -0.183508/
c *   1.6872 mb:
      data (coef(i,12,3),i=0,nx)/
     +  338.109192,    0.077756,   -0.167569,   -0.228182,    0.081726,
     +   -0.109376,   -0.114042,    0.249524,    0.124284,   -0.253638/
c *   2.1526 mb:
      data (coef(i,13,3),i=0,nx)/
     +  330.158142,    0.125628,   -0.127452,   -0.284897,    0.018188,
     +   -0.073173,   -0.058458,    0.164222,    0.202521,   -0.283837/
c *   2.7009 mb:
      data (coef(i,14,3),i=0,nx)/
     +  288.525238,    0.208034,   -0.113437,   -0.262148,   -0.070341,
     +   -0.019104,    0.050349,   -0.000575,    0.297252,   -0.249445/
c *   3.3398 mb:
      data (coef(i,15,3),i=0,nx)/
     +  247.402054,    0.288580,   -0.105926,   -0.234689,   -0.143902,
     +    0.044019,    0.059252,   -0.048267,    0.325860,   -0.189348/
c *   4.0770 mb:
      data (coef(i,16,3),i=0,nx)/
     +  207.174408,    0.365898,   -0.103597,   -0.202744,   -0.204937,
     +    0.116567,   -0.020110,    0.001814,    0.300148,   -0.106910/
c *   4.9204 mb:
      data (coef(i,17,3),i=0,nx)/
     +  172.771362,    0.425351,   -0.097642,   -0.168375,   -0.266588,
     +    0.200371,   -0.115398,    0.021846,    0.299176,   -0.026388/
c *   5.8776 mb:
      data (coef(i,18,3),i=0,nx)/
     +  170.819962,    0.377608,    0.003744,   -0.188331,   -0.277131,
     +    0.162590,   -0.022318,   -0.064583,    0.305203,   -0.023369/
c *   6.9567 mb:
      data (coef(i,19,3),i=0,nx)/
     +  171.848145,    0.322555,    0.108887,   -0.212190,   -0.282633,
     +    0.115752,    0.083180,   -0.156449,    0.311570,   -0.027389/
c *   8.1655 mb:
      data (coef(i,20,3),i=0,nx)/
     +  178.900009,    0.317759,    0.072287,   -0.147714,   -0.288088,
     +    0.079424,    0.109033,   -0.207073,    0.293082,   -0.003934/
c *   9.5119 mb:
      data (coef(i,21,3),i=0,nx)/
     +  185.850937,    0.315015,    0.032178,   -0.082944,   -0.293293,
     +    0.045129,    0.130805,   -0.253892,    0.274527,    0.019458/
c *  11.0038 mb:
      data (coef(i,22,3),i=0,nx)/
     +  181.033875,    0.342789,   -0.008198,   -0.050784,   -0.276437,
     +    0.028043,    0.139765,   -0.267162,    0.249229,    0.043197/
c *  12.6492 mb:
      data (coef(i,23,3),i=0,nx)/
     +  170.696304,    0.384559,   -0.047857,   -0.034871,   -0.249395,
     +    0.019533,    0.142417,   -0.264131,    0.221237,    0.066608/
c *  14.4559 mb:
      data (coef(i,24,3),i=0,nx)/
     +  160.792282,    0.424578,   -0.085853,   -0.019626,   -0.223487,
     +    0.011380,    0.144958,   -0.261227,    0.194418,    0.089037/
c *  16.4318 mb:
      data (coef(i,25,3),i=0,nx)/
     +  150.499008,    0.461093,   -0.122477,   -0.003129,   -0.196627,
     +    0.000506,    0.147210,   -0.251144,    0.167289,    0.109366/
c *  18.5847 mb:
      data (coef(i,26,3),i=0,nx)/
     +  140.299698,    0.495449,   -0.157735,    0.013451,   -0.170036,
     +   -0.011132,    0.149301,   -0.238611,    0.140674,    0.128437/
c *  20.9224 mb:
      data (coef(i,27,3),i=0,nx)/
     +  133.652267,    0.528523,   -0.174538,    0.028105,   -0.154952,
     +   -0.021779,    0.150992,   -0.227335,    0.124328,    0.123706/
c *  23.4526 mb:
      data (coef(i,28,3),i=0,nx)/
     +  132.214737,    0.560411,   -0.163866,    0.040182,   -0.156891,
     +   -0.031171,    0.152119,   -0.217697,    0.123112,    0.082945/
c *  26.1829 mb:
      data (coef(i,29,3),i=0,nx)/
     +  130.282120,    0.591975,   -0.157004,    0.055184,   -0.159140,
     +   -0.038202,    0.157241,   -0.213196,    0.124561,    0.042131/
c *  29.1210 mb:
      data (coef(i,30,3),i=0,nx)/
     +  127.686752,    0.623526,   -0.154959,    0.074146,   -0.161817,
     +   -0.042282,    0.167575,   -0.215254,    0.129460,    0.000725/
c *  32.2744 mb:
      data (coef(i,31,3),i=0,nx)/
     +  108.781189,    0.698895,   -0.147134,    0.068980,   -0.139383,
     +   -0.037117,    0.145935,   -0.184685,    0.111908,   -0.009194/
c *  35.6505 mb:
      data (coef(i,32,3),i=0,nx)/
     +   84.030518,    0.789494,   -0.137261,    0.054729,   -0.107821,
     +   -0.028531,    0.112540,   -0.142284,    0.086147,   -0.006936/
c *  39.2566 mb:
      data (coef(i,33,3),i=0,nx)/
     +   60.058582,    0.877243,   -0.127698,    0.040927,   -0.077253,
     +   -0.020216,    0.080196,   -0.101216,    0.061195,   -0.004749/
c *  43.1001 mb:
      data (coef(i,34,3),i=0,nx)/
     +   36.821239,    0.962303,   -0.118428,    0.027547,   -0.047621,
     +   -0.012155,    0.048844,   -0.061407,    0.037009,   -0.002629/
c *  47.1882 mb:
      data (coef(i,35,3),i=0,nx)/
     +   14.277436,    1.044823,   -0.109435,    0.014568,   -0.018874,
     +   -0.004335,    0.018428,   -0.022788,    0.013544,   -0.000572/

	end

	subroutine cofit3_64(xx1,xx2,xx3,yy1,yy2,yy3, c0,c1,c2)
c * Obtain coefficients for 3-point parabolic fit

        implicit real*8 (A-H,O-Z)

	x1=xx1
	x2=xx2
	x3=xx3

	y1=yy1
	y2=yy2
	y3=yy3

	x12=x1*x1
	x22=x2*x2
	x32=x3*x3

	t1=x2*x32-x3*x22
	t2=-(x1*x32-x3*x12)
	t3=x1*x22-x2*x12

	det=t1+t2+t3

	c0=(y1*t1+y2*t2+y3*t3)/det
	c1=((y2*x32-y3*x22)-(y1*x32-y3*x12)+(y1*x22-y2*x12))/det
	c2=((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))/det

	return
	end
