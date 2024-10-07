        subroutine conv_soh(soh,scana,scdis)

c $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/conv_soh.f,v 1.2 1996/04/17 18:49:03 seawifsd Exp seawifsd $
c $Log: conv_soh.f,v $
c Revision 1.2  1996/04/17 18:49:03  seawifsd
c fixed a bug about the wrong dimension size.(ityp_acs(3) -> ityp_acs(40))
c
c Revision 1.1  1996/04/17 18:48:19  seawifsd
c Initial revision
c
c
c  Subroutine to convert SOH telemetry values.  Complete list of converted 
c  data is given in acs_block.f

c  February 2, 1994 by Frederick S. Patt

c  Modified to add multiple data types to analog specifications - 
c   December 7, 1994 by Frederick S. Patt

        real*4 scana(40)
        byte soh(776),scdis(40)

c  Location and conversions specifications are located in acs_comm
        integer*4 loc_acs(2,40)
        integer*4 dis_acs(3,40)
        integer*4 ityp_acs(40)
        real*4 con_acs(2,40)
c       common /acs_comm/ con_acs, loc_acs, dis_acs, ityp_acs
c MDM Oct. 15, 2004 making the block data routine acs_block inlined to
c avoid problems with the linker not picking it up...
c see acs_block.f for descriptions
c  Analog data start byte and length in SOH packet
        data loc_acs/   113,8,          !Orbit X position 
     *                  121,8,          !Orbit Y position
     *                  129,8,          !Orbit Z position
     *                  137,8,          !Orbit X velocity
     *                  145,8,          !Orbit Y velocity 
     *                  153,8,          !Orbit Z velocity 
     *                  265,4,          !Attitude yaw angle
     *                  273,4,          !Attitude roll angle
     *                  269,4,          !Attitude pitch angle
     *                  553,2,          !Sun sensor 1 angle 1
     *                  555,2,          !Sun sensor 1 angle 2
     *                  565,2,          !Sun sensor 2 angle 1
     *                  567,2,          !Sun sensor 2 angle 2
     *                  577,2,          !Sun sensor 3 angle 1
     *                  579,2,          !Sun sensor 3 angle 2
     *                  613,2,          !Earth scanner 1 phase  
     *                  615,2,          !Earth scanner 1 width  
     *                  625,2,          !Earth scanner 2 phase  
     *                  627,2,          !Earth scanner 2 width
     *                  57,2,           !GPS DOP value
     *                  31,4,           !GPS time tag fractional second
     *                  497,2,          !Momentum wheel 1 speed
     *                  503,2,          !Momentum wheel 1 current
     *                  509,2,          !Momentum wheel 2 speed
     *                  515,2,          !Momentum wheel 2 current
     *                  425,4,          !Torquer rod 1 x command levels 
     *                  429,4,          !Torquer rod 1 y command levels 
     *                  433,4,          !Torquer rod 1 z command levels 
     *                  437,4,          !Torquer rod 2 x command levels 
     *                  441,4,          !Torquer rod 2 y command levels 
     *                  445,4,          !Torquer rod 2 z command levels 
     *                  585,2,          !Magnetometer 1 x-axis
     *                  587,2,          !Magnetometer 1 y-axis
     *                  589,2,          !Magnetometer 1 z-axis
     *                  591,2,          !Magnetometer 1 r-axis
     *                  201,4,          !Attitude yaw rate
     *                  209,4,          !Attitude roll rate
     *                  205,4,          !Attitude pitch rate
     *                  4*0             /

c  Data types
        data ityp_acs/  5,              !Orbit X position 
     *                  5,              !Orbit Y position
     *                  5,              !Orbit Z position
     *                  5,              !Orbit X velocity
     *                  5,              !Orbit Y velocity 
     *                  5,              !Orbit Z velocity 
     *                  4,              !Attitude yaw angle
     *                  4,              !Attitude roll angle
     *                  4,              !Attitude pitch angle
     *                  1,              !Sun sensor 1 angle 1
     *                  1,              !Sun sensor 1 angle 2
     *                  1,              !Sun sensor 2 angle 1
     *                  1,              !Sun sensor 2 angle 2
     *                  1,              !Sun sensor 3 angle 1
     *                  1,              !Sun sensor 3 angle 2
     *                  2,              !Earth scanner 1 phase  
     *                  2,              !Earth scanner 1 width  
     *                  2,              !Earth scanner 2 phase  
     *                  2,              !Earth scanner 2 width
     *                  1,              !GPS DOP value
     *                  3,              !GPS time tag fractional second
     *                  2,              !Momentum wheel 1 speed
     *                  1,              !Momentum wheel 1 current
     *                  2,              !Momentum wheel 2 speed
     *                  1,              !Momentum wheel 2 current
     *                  4,              !Torquer rod 1 x command levels 
     *                  4,              !Torquer rod 1 y command levels 
     *                  4,              !Torquer rod 1 z command levels 
     *                  4,              !Torquer rod 2 x command levels 
     *                  4,              !Torquer rod 2 y command levels 
     *                  4,              !Torquer rod 2 z command levels 
     *                  2,              !Magnetometer 1 x-axis
     *                  2,              !Magnetometer 1 y-axis
     *                  2,              !Magnetometer 1 z-axis
     *                  2,              !Magnetometer 1 r-axis
     *                  4,              !Attitude yaw rate
     *                  4,              !Attitude roll rate
     *                  4,              !Attitude pitch rate
     *                  2*0             /

c  Linear analog conversion coefficients
        data con_acs/
     *          0.001, 0.,              !Orbit X position 
     *          0.001, 0.,              !Orbit Y position
     *          0.001, 0.,              !Orbit Z position
     *          0.001, 0.,              !Orbit X velocity
     *          0.001, 0.,              !Orbit Y velocity 
     *          0.001, 0.,              !Orbit Z velocity 
     *          1., 0.,                 !Attitude yaw angle
     *          1., 0.,                 !Attitude roll angle
     *          1., 0.,                 !Attitude pitch angle
     *          2.00225e-04, -2.0503,   !Sun sensor 1 angle 1
     *          2.00225e-04, -2.0503,   !Sun sensor 1 angle 2
     *          2.00225e-04, -2.0503,   !Sun sensor 2 angle 1
     *          2.00225e-04, -2.0503,   !Sun sensor 2 angle 2
     *          2.00225e-04, -2.0503,   !Sun sensor 3 angle 1
     *          2.00225e-04, -2.0503,   !Sun sensor 3 angle 2
     *          0.005493164, 0.0,       !Earth scanner 1 phase  
     *          0.005493164, 0.0,       !Earth scanner 1 width  
     *          0.005493164, 0.0,       !Earth scanner 2 phase  
     *          0.005493164, 0.0,       !Earth scanner 2 width
     *          0.1, 0.0,               !GPS DOP value
     *          1.0E-9, 0.0,            !GPS time tag fractional second
     *          0.2, 0.0,               !Momentum wheel 1 speed
     *          0.001, 0.0,             !Momentum wheel 1 current
     *          0.2, 0.0,               !Momentum wheel 2 speed
     *          0.001, 0.0,             !Momentum wheel 2 current
     *          1.0, 0.0,               !Torquer rod 1 x command levels 
     *          1.0, 0.0,               !Torquer rod 1 y command levels 
     *          1.0, 0.0,               !Torquer rod 1 z command levels 
     *          1.0, 0.0,               !Torquer rod 2 x command levels 
     *          1.0, 0.0,               !Torquer rod 2 y command levels 
     *          1.0, 0.0,               !Torquer rod 2 z command levels 
     *          2.0, 0.0,               !Magnetometer 1 x-axis
     *          2.0, 0.0,               !Magnetometer 1 y-axis
     *          2.0, 0.0,               !Magnetometer 1 z-axis
     *          2.0, 0.0,               !Magnetometer 1 r-axis
     *          1.0, 0.0,               !Attitude yaw rate
     *          1.0, 0.0,               !Attitude roll rate
     *          1.0, 0.0,               !Attitude pitch rate
     *                  4*0             /

c  Discrete data start byte, bit and bit length in SOH packet
        data dis_acs/
     *          551, 1, 2,              !Sun sensor 1 presence
     *          563, 1, 2,              !Sun sensor 2 presence
     *          575, 1, 2,              !Sun sensor 3 presence
     *          549, 1, 8,              !Sun sensor 1 status
     *          561, 1, 8,              !Sun sensor 2 status
     *          573, 1, 8,              !Sun sensor 3 status
     *          609, 1, 8,              !Earth scanner 1 status 1
     *          612, 1, 8,              !Earth scanner 1 status 2
     *          621, 1, 8,              !Earth scanner 2 status 1
     *          624, 1, 8,              !Earth scanner 2 status 2
     *          60, 1, 8,               !Number of GPS satellites visible
     *          61, 1, 8,               !Number of GPS satellites tracked
     *          111, 1, 8,              !GPS receiver status
     *          25, 1, 8,               !GPS time tag year (byte 1)
     *          26, 1, 8,               !GPS time tag year (byte 2)
     *          23, 1, 8,               !GPS time tag month
     *          24, 1, 8,               !GPS time tag day-of-month
     *          27, 1, 8,               !GPS time tag hour
     *          28, 1, 8,               !GPS time tag minute
     *          29, 1, 8,               !GPS time tag second
     *          557, 1, 8,              !Sun sensor 1 time tag (byte 1)
     *          558, 1, 8,              !Sun sensor 1 time tag (byte 2)
     *          559, 1, 8,              !Sun sensor 1 time tag (byte 3)
     *          560, 1, 8,              !Sun sensor 1 time tag (byte 4)
     *          569, 1, 8,              !Sun sensor 2 time tag (byte 1)
     *          570, 1, 8,              !Sun sensor 2 time tag (byte 2)
     *          571, 1, 8,              !Sun sensor 2 time tag (byte 3)
     *          572, 1, 8,              !Sun sensor 2 time tag (byte 4)
     *          581, 1, 8,              !Sun sensor 3 time tag (byte 1)
     *          582, 1, 8,              !Sun sensor 3 time tag (byte 2)
     *          583, 1, 8,              !Sun sensor 3 time tag (byte 3)
     *          584, 1, 8,              !Sun sensor 3 time tag (byte 4)
     *          617, 1, 8,              !Earth scanner 1 time tag (byte 1)
     *          618, 1, 8,              !Earth scanner 1 time tag (byte 2)
     *          619, 1, 8,              !Earth scanner 1 time tag (byte 3)
     *          620, 1, 8,              !Earth scanner 1 time tag (byte 4)
     *          629, 1, 8,              !Earth scanner 2 time tag (byte 1)
     *          630, 1, 8,              !Earth scanner 2 time tag (byte 2)
     *          631, 1, 8,              !Earth scanner 2 time tag (byte 3)
     *          632, 1, 8               !Earth scanner 2 time tag (byte 4)
     *                  /


c  Convert analog data 
        do i=1,40
          if (ityp_acs(i).eq.1) then
            call read_analog(scana(i),loc_acs(1,i),con_acs(1,i),soh)
          else if (ityp_acs(i).eq.2) then
            call read_short(scana(i),loc_acs(1,i),con_acs(1,i),soh)
          else if (ityp_acs(i).eq.3) then
            call read_long(scana(i),loc_acs(1,i),con_acs(1,i),soh)
          else if (ityp_acs(i).eq.4) then
            call read_float(scana(i),loc_acs(1,i),con_acs(1,i),soh)
          else if (ityp_acs(i).eq.5) then
            call read_double(scana(i),loc_acs(1,i),con_acs(1,i),soh)
          end if
        end do

c  Convert discrete data        
        do i=1,40
          if (dis_acs(1,i).gt.0) 
     *      call read_discrete(scdis(i),dis_acs(1,i),soh)
        end do
        
        return
        end
