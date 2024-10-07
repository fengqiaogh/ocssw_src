        subroutine sctime(mtime,iyr,iday,msec)
        integer*2 mtime(4)
        
c  Get Julian day from time tag
        jday = mtime(1)*8 + mtime(2)/128 + 2449001
        call jdate(jday,iyr,iday)

c  Get milliseconds of day and convert to UTC
        msec = 1048576*mod(mtime(2),128) + 1024*mtime(3) + mtime(4)

        return
        end

