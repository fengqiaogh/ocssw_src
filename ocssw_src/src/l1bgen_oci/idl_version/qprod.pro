pro qprod,q1,q2,q3

; Compute the product of two quaternions q3 = q1*q2

q3 = q1

q3(0,*) = q1(0,*)*q2(3,*) + q1(1,*)*q2(2,*) - q1(2,*)*q2(1,*) + q1(3,*)*q2(0,*)
q3(1,*) = -q1(0,*)*q2(2,*) + q1(1,*)*q2(3,*) + q1(2,*)*q2(0,*) + q1(3,*)*q2(1,*)
q3(2,*) = q1(0,*)*q2(1,*) - q1(1,*)*q2(0,*) + q1(2,*)*q2(3,*) + q1(3,*)*q2(2,*)
q3(3,*) = -q1(0,*)*q2(0,*) - q1(1,*)*q2(1,*) - q1(2,*)*q2(2,*) + q1(3,*)*q2(3,*)

return
end
