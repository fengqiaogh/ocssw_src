function iqr,dat

; IDL function to compute the inter-quartile range (IQR) of a data array

;  Name	Type 	I/O	Description
;
;  dat		Any	 I	Input data array

id = sort(dat)
nd = n_elements(dat)
iq1 = dat(id(nd/4))
iq2 = dat(id(3*nd/4))
iqr = iq2 - iq1

return,iqr
end
