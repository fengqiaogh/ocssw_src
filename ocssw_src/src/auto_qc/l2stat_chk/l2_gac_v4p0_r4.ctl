# Control file for level 2 GAC data files
# Specifies what portions of the files should be verified
#  W. Robinson, GSC, 10 Dec 99  update for use with version 4.0
#  using the new L2 file from MSl12
# W. Robinson, SAIC, 27 Jun 2002 modify to check chlor_a as float,
#         to extend valid chlor range and to admit new flags
#
# % Gross value check
#              error thresh (%) Low thresh	High thresh
L2GCHK1       100.0		  0.0	    	 32.767003
L2GCHK2	      100.0               0.0            32.767003
L2GCHK3	      100.0               0.0            32.767003
L2GCHK4	      100.0               0.0            32.766003
L2GCHK5	      100.0               0.0            32.767003
L2GCHK6	      100.0               0.0            3.2767
L2GCHK7       100.0               0.0            6.4768     
L2GCHK8	      100.0               0.0            100.1
L2GCHK9       100.0               0.0            6.5534001
L2GCHK10      100.0               0.0            1.1600001
L2GCHK11      100.0               0.0            0.3250001
#
# % statistics check
#      
L2STAT1  
L2STAT2 
L2STAT3
L2STAT4  
L2STAT5 
L2STAT6
L2STAT7  
L2STAT8 
L2STAT9 
L2STAT10
L2STAT11
#
# Flag Percentage check
#             Low thresh	High thresh
L2FLGCK1  	0.0               100.0
L2FLGCK2  	0.0               100.0
L2FLGCK3  	0.0               100.0
L2FLGCK4  	0.0               100.0
L2FLGCK5  	0.0               100.0
L2FLGCK6  	0.0               100.0
L2FLGCK7  	0.0               100.0
L2FLGCK8  	0.0               100.0
L2FLGCK9  	0.0               100.0
L2FLGCK10 	0.0               100.0
L2FLGCK11 	0.0               100.0
L2FLGCK12 	0.0               100.0
L2FLGCK13 	0.0               100.0
L2FLGCK14 	0.0               100.0
L2FLGCK15 	0.0               100.0
L2FLGCK16 	0.0               100.0
L2FLGCK17       0.0               100.0
L2FLGCK18       0.0               100.0
L2FLGCK19       0.0               100.0
L2FLGCK20       0.0               100.0
L2FLGCK21       0.0               100.0
L2FLGCK22       0.0               100.0
L2FLGCK23       0.0               100.0
L2FLGCK24       0.0               100.0
L2FLGCK25       0.0               100.0
L2FLGCK26       0.0               100.0
L2FLGCK27       0.0               100.0
#
#  these are not in use yet
#L2FLGCK28       0.0               100.0
#L2FLGCK29       0.0               100.0
#L2FLGCK30       0.0               100.0
#L2FLGCK31       0.0               100.0
L2FLGCK32       0.0               100.0
#
# End of file
