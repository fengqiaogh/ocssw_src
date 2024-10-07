        subroutine matvec(xm,v1,v2)
c  Computes 3-vector v2 as product of 3x3 matrix xm and 3-vector v1

c       Program written by:     Frederick S. Patt
c                               General Sciences Corporation
c                               October 19, 1992

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  xm(3,3)      R*4      I      Input Matrix
c  v1(3)        R*4      I      Input Vector
c  v2(3)        R*4      I      Output Vector

        real*4 v1(3),v2(3),xm(3,3)
        do i=1,3
            v2(i) = 0.0
            do j=1,3
                v2(i) = v2(i) + xm(i,j)*v1(j) 
            end do
        end do
        return
        end
