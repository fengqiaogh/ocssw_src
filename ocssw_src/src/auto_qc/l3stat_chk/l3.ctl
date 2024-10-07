# Control file for level 3 Bin data files
# Specifies what portions of the files should be verified
#
# Data bin percentage check
#
#           Error Low thresh   Error High thresh      
L3DAT             20.00               80.0
# Weight check
#
L3WT
#
# % Gross value check (make % 30 if possible later)
#            error thresh (%) Low thresh     High thresh
L3GCHK1       50.0                0.6           2.8
L3GCHK2       50.0                0.3           3.0
L3GCHK3       50.0                0.2           2.0
L3GCHK4       50.0                0.25           1.00
L3GCHK5       50.0                0.1           0.70
L3GCHK6       50.0                0.0           2.0
L3GCHK7       50.0                0.0           2.0
L3GCHK8       50.0                0.0           2.0
L3GCHK9       50.0                0.01           1.0
L3GCHK10      50.0                0.0           20.0
L3GCHK11      50.0                0.5           1.5
L3GCHK12      50.0                0.0           0.3

# % Statistical check
#      
L3STAT1  
L3STAT2 
L3STAT3
L3STAT4  
L3STAT5 
L3STAT6
L3STAT7  
L3STAT8 
L3STAT9 
L3STAT10
L3STAT11
L3STAT12
#
# Climatology check
#             Thresh1L   Thresh1H  Thresh2L  Thresh1H  Thresh3L  Thresh3H 
L3CCHK1        100.0      100.0     100.0     100.0     100.0     100.0
L3CCHK2        100.0      100.0     100.0     100.0     100.0     100.0       
L3CCHK3        100.0      100.0     100.0     100.0     100.0     100.0      
L3CCHK4        100.0      100.0     100.0     100.0     100.0     100.0
L3CCHK5        100.0      100.0     100.0     100.0     100.0     100.0
L3CCHK6        100.0      100.0     100.0     100.0     100.0     100.0
L3CCHK7        100.0      100.0     100.0     100.0     100.0     100.0
L3CCHK8        100.0      100.0     100.0     100.0     100.0     100.0
L3CCHK9        100.0      100.0     100.0     100.0     100.0     100.0
L3CCHK10       100.0      100.0     100.0     100.0     100.0     100.0
L3CCHK11       100.0      100.0     100.0     100.0     100.0     100.0
L3CCHK12       100.0      100.0     100.0     100.0     100.0     100.0
#
#  for later, use these bounds
#L3CCHK1         25.0       25.0       8.0       8.0       1.5       1.5     
#
############### End of file
