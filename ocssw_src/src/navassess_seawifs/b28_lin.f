      subroutine b28_lin(scan,npix,b2,b8)

c  This subroutine removes the effect of the bilinear gains from 
c   SeaWiFS bands 2 and 8

c  Modified to reflect the redefinition of the knees.  April 20, 2000.

      save
      real*4 b2(*),b8(*)
      real*4 b2tab(1024),b8tab(1024)
      real*4 b2k(5),b2s(5),b2r(5),b8k(5),b8s(5),b8r(5)
      integer*2 scan(8,*)
c      data b1k/21.,814.4,814.6,818.5,1023./
c      data b6k/23.,786.3,786.5,787.4,1023./
      data b2k/18.,808.27,809.16,811.31,1023.0/
      data b8k/21.,783.99,785.22,786.52,1023.0/
c      data b1s/0.013733,0.013733,0.020049,0.037170,0.240156/
c      data b6s/0.004254,0.004254,0.006344,0.012422,0.219265/
      data b2s/0.013426,0.013426,0.022472,0.037209,0.274/
      data b8s/0.0022241,0.0022241,0.0032520,0.0061538,0.139345/
c      data b1r/0.,10.899,10.903,11.049,60.159/
c      data b6r/0.,3.247,3.248,3.259,54.926/
      data b2r/0.,10.61,10.63,10.71,68.713/
      data b8r/0.,1.697,1.701,1.709,34.6613/
      logical first
      data first/.true./

c  If first time through, set up counts-to-radiances table
      if (first) then
        first = .false.
        i1 = 1
        i2 = b2k(2)+1
        do i=i1,i2
          b2tab(i) = b2s(2)*(i - b2k(1) - 1)
        end do
        do l=3,5
          i1 = i2 + 1
          i2 = b2k(l) + 1
          if (i2.ge.i1) then
            do i=i1,i2
              b2tab(i) = b2s(l)*(i - b2k(l-1) - 1) + b2r(l-1)
            end do
          end if
        end do

        i1 = 1
        i2 = b8k(2)+1
        do i=i1,i2
          b8tab(i) = b8s(2)*(i - b8k(1) - 1)
        end do
        do l=3,5
          i1 = i2 + 1
          i2 = b8k(l) + 1
          if (i2.ge.i1) then
            do i=i1,i2
              b8tab(i) = b8s(l)*(i - b8k(l-1) - 1) + b8r(l-1)
            end do
          end if
        end do
      end if
          
c  Loop through pixels

      do i=1,npix

c  Linearize Band 2

        b2(i) = b2tab(scan(2,i)+1)

c  Linearize Band 8

        b8(i) = b8tab(scan(8,i)+1)

      end do

      return
      end

