c       data set invert (entered from hardcopy, 10 Jan 93, JWhiting)
      subroutine invert(a,b,n,l,c,ier)

c  invert -- inverts a matrix in place and solves a set of 
c            simultaneous linear equations.

c  input: 

c  a(l,l)    matrix to be inverted (inverse is returned in this 
c            same array).
c  b(n)      vector of right hand sides of linear equations
c  n         number of equations to be solved (1-20)
c  l         dimension of a, as declared in calling prog (l.ge.n)

c  output:

c  a(l,l)    inverse of original matrix
c  c(n)      vector of solution to linear equations
c  ier       return code
c            =0, normal return
c            =1, singular matrix (matrix a may be altered)
c            =2, invalid input (n.lt.1, or n.gt.20, or l.lt.n)
c                (matrix a unaltered)

c  external references: none

c    M. Shear, Computer Sciences Corp. 3-30-75
c    Minor revision to subroutine LINEAR. Name changed and error 
c    return ier=2 inserted.  January 1976

c  Method
c    Gauss-Jordan elimination with normalization and optimal pivoting.

c  Reference
c    Shan S. Kuo, Numerical Methods and Computers, 
c    Addison-Wesley, 1965. pp. 168-169.

c  Notes
c    1. If only matrix inversion is desired, the vector of right hand 
c       sides may be set to zeroes.
c    2. Execution time varies approx. as n**3. For n=20, execution 
c       time is about .1 seconds (IBM S/360-95)

      implicit real*8 (a-h,o-z)
      dimension a(l,l), b(1), c(1)

c Internal declarations
      dimension ipvot(20), index(20,2)

c Initialize
      if (n.lt.1 .or. n.gt.20 .or. l.lt.n) goto 902
      do 10 i=1,n
        c(i)=b(i)
   10 ipvot(i)=0

c Loop on eliminations
      do 500 npass=1,n

c Determine pivot element
        amax=0.d0
        do 30 j=1,n
          if (ipvot(j).ne.0) goto 30
          do 20 k=1,n
            if (ipvot(k).ne.0) goto 20
            if (dabs(amax).ge.dabs(a(j,k))) goto 20
            irow=j
            icol=k
            amax=a(j,k)
   20     continue
   30   continue
        if (amax.eq.0.d0) goto 900
        ipvot(icol)=1

c Normalization and exchange of rows
        a(irow,icol)=1.d0
        do 40 i=1,n
          t=a(irow,i)
          a(irow,i)=a(icol,i)
   40     a(icol,i)=t/amax
        t=c(irow)
        c(irow)=c(icol)
        c(icol)=t/amax
        index(npass,1)=irow
        index(npass,2)=icol

c Perform elimination
        do 100 i=1,n
          if (i.eq.icol) goto 100
          t=a(i,icol)
          a(i,icol)=0.d0
          c(i)=c(i)-c(icol)*t
          do 90 j=1,n
   90       a(i,j)=a(i,j)-a(icol,j)*t
  100   continue
  500 continue

c Reorder the inverted matrix by exchanging columns
      do 600 i=1,n
        k=n-i+1
        irow=index(k,1)
        icol=index(k,2)
        do 550 j=1,n
          t=a(j,irow)
          a(j,irow)=a(j,icol)
  550     a(j,icol)=t
  600 continue
      ier=0
      goto 999

  900 ier=1
      goto 999

  902 ier=2

  999 return
      end
