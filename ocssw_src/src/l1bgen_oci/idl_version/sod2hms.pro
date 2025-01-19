pro sod2hms,secd,ih,mn,sec
ih = fix(secd/3600.d0)
mn = fix((secd mod 3600.d0)/60.d0)
sec = secd mod 60.d0
;print,ih,mn,sec,format='(2(i02,":"),i2)'
return
end
