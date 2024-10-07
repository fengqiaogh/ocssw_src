def jd(i,j,k):
# converted from: jd.f, Written by Frederick S. Patt, GSC, November 4, 1992
# Liang Hong, 2017/6/21
# Liang Hong, 2019/11/6:  made it python 3 compatible, all calculations are integer based
	
# This function converts a calendar date to the corresponding Julian
# day starting at noon on the calendar date.  The algorithm used is
# from Van Flandern and Pulkkinen, Ap. J. Supplement Series 41, 
# November 1979, p. 400.
	
#    Arguments
#  
#    Name    Type    I/O     Description
#    ----    ----    ---     -----------
#    i       I*4      I      Year - e.g. 1970
#    j       I*4      I      Month - (1-12)
#    k       I*4      I      Day  - (1-31)
#    jd      I*4      O      Julian day
	
	import numpy as np
	i = np.array(i)
	j = np.array(j)
	k = np.array(k)
	jd = 367*i - 7*(i+(j+9)//12)//4 + 275*j//9 + k + 1721014
	return jd
