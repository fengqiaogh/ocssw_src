        pro mtoq,rm,q

;  Convert direction cosine matrix to equivalent quaternion

        q = dblarr(4)
        e = dblarr(3)

;  Compute Euler angle
        cphi = (rm(0,0)+rm(1,1)+rm(2,2)-1.d0)/2.d0
        if (abs(cphi) lt 0.98) then begin
            phi = acos(cphi)
        endif else begin
            ssphi = ((rm(0,1)-rm(1,0))^2 + (rm(2,0)-rm(0,2))^2 $
              + (rm(1,2)-rm(2,1))^2)/4.d0
            phi = asin(sqrt(ssphi))
            if (cphi lt 0) then phi = !pi - phi
        endelse

;  Compute Euler axis
        e(0) = (rm(1,2)-rm(2,1))/(sin(phi)*2.d0)
        e(1) = (rm(2,0)-rm(0,2))/(sin(phi)*2.d0)
        e(2) = (rm(0,1)-rm(1,0))/(sin(phi)*2.d0)
        e = e/sqrt(total(e*e))

;  Compute quaternion
        q(0) = e(0)*sin(phi/2.d0)
        q(1) = e(1)*sin(phi/2.d0)
        q(2) = e(2)*sin(phi/2.d0)
        q(3) = cos(phi/2.d0)

        return
        end
