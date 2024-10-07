        subroutine jddate(jd,i,j,k)
C
C       This routine computes the calendar date corresponding to
C       a given Julian day.  This code was brazenly copied from
C       a routine written by Myron Shear for CSC on Julian Day 1.
C       
C       ARGUMENT        TYPE    I/O     DESCRIPTION     
C       __________________________________________________________
C        JD             I*4      I      Julian Day (reference Jan 1, 4713 BC)
C        I              I*4      O      Year 
C        J              I*4      O      Month   
C        K              I*4      0      Day of Month
C
        l = jd + 68569
        n = 4*l/146097
        l = l - (146097*n + 3)/4
        i = 4000*(l+1)/1461001
        l = l - 1461*i/4 + 31
        j = 80*l/2447
        k = l - 2447*j/80
        l = j/11
        j = j + 2 - 12*l
        i = 100*(n-49) + i + l
        return
        end     
