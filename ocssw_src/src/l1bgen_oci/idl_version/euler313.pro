	pro euler313,a,xm
;  Computes coordinate transformation matrix corresponding to Euler 
;   sequence; assumes order of rotations is 3,1,3
;
;  Reference:  Wertz, Appendix E

xma=fltarr(3,3,3)&xm=fltarr(3,3)&xmm=fltarr(3,3)
 
        if (n_elements(iord) eq 0) then iord = [1,2,3]
	c1=cos(a(0)/!radeg)
	s1=sin(a(0)/!radeg)
	c2=cos(a(1)/!radeg)
	s2=sin(a(1)/!radeg)
	c3=cos(a(2)/!radeg)
	s3=sin(a(2)/!radeg)
;  Convert individual rotations to matrices
	xma(2,2,0)=1.d0
	xma(1,1,0)=c1
	xma(0,0,0)=c1
	xma(0,1,0)=s1
	xma(1,0,0)=-s1
	xma(0,0,1)=1.d0
	xma(1,1,1)=c2
	xma(2,2,1)=c2
	xma(1,2,1)=s2
	xma(2,1,1)=-s2
	xma(2,2,2)=1.d0
	xma(1,1,2)=c3
	xma(0,0,2)=c3
	xma(0,1,2)=s3
	xma(1,0,2)=-s3
;  Compute total rotation as xm3*xm2*xm1
	xmm=xma(*,*,iord(1)-1)#xma(*,*,iord(0)-1)
	xm=xma(*,*,iord(2)-1)#xmm
	return
	end

