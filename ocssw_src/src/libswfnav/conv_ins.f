        subroutine conv_ins(inst,insta,instd)
c  Subroutine to convert instrument telemetry values. 

c  February 2, 1994 by Frederick S. Patt
c
c Modification History
c
c  Changed integer*1 to byte for Sun OS compatibility, B. A. Franz,
c  GAC, November 14, 1997.


        real*4 insta(40)
        byte inst(88),instd(32)
        byte ilsb(44)

c  Location and conversions specifications are located in acs_con.fin
        integer*4 loc_ins(2,40)
        integer*4 dis_ins(3,32)
        real*4 con_ins(2,40)
c       common /ins_comm/ loc_ins, dis_ins, con_ins
c MDM Oct. 15, 2004 making the block data routine ins_block inlined to
c avoid problems with the linker not picking it up...
c see ins_block.f for descriptions
c  Analog data start word and length in instrument packet
        data loc_ins/
     *                  8, 1,           !Band 1/2 FPA Temperature
     *                  9, 1,           !Band 3/4 FPA Temperature
     *                  10, 1,          !Band 5/6 FPA Temperature
     *                  11, 1,          !Band 7/8 FPA Temperature
     *                  12, 1,          !Telescope Motor Temperature
     *                  13, 1,          !Tilt Base Temperature
     *                  14, 1,          !Tilt Platform Temperature
     *                  15, 1,          !Half Angle Motor Temperature
     *                  16, 1,          !Power Supply A Input Current
     *                  17, 1,          !Power Supply B Input Current
     *                  18, 1,          !+15 V Analog Power Voltage
     *                  19, 1,          !-15 V Analog Power Voltage
     *                  20, 1,          !+5 V Logic Power Voltage
     *                  21, 1,          !Power Supply Temperature
     *                  22, 1,          !B1/B2 Postamp Temperature
     *                  23, 1,          !Servo Drive Temperature
     *                  24, 1,          !+30 V Servo Power Voltage
     *                  25, 1,          !+21 V Servo Power Voltage
     *                  26, 1,          !-21 V Servo Power Voltage
     *                  27, 1,          !+5 V Servo Power Voltage
     *                  28, 1,          !Angular Comp. Phase Error
     *                  29, 1,          !Tilt Platform Position
     *                  30, 1,          !Tilt Base Position
     *                  31, 1,          !+28 V Heater Power
     *                  32, 1,          !Telescope A Motor Current
     *                  33, 1,          !Telescope B Motor Current
     *                  34, 1,          !Half Angle A Motor Current
     *                  35, 1,          !Half Angle B Motor Current
     *                  36, 1,          !Servo A Phase Error
     *                  37, 1,          !Servo B Phase Error
     *                  38, 1,          !Angular Comp. A Motor Current
     *                  39, 1,          !Angular Comp. B Motor Current
     *                  16*0/           !Spares

c  Linear analog conversion coefficients
        data con_ins/
     *          -0.2667,  66.667,       !Band 1/2 FPA Temperature
     *          -0.2667,  66.667,       !Band 3/4 FPA Temperature
     *          -0.2667,  66.667,       !Band 5/6 FPA Temperature
     *          -0.2667,  66.667,       !Band 7/8 FPA Temperature
     *          -0.2667,  66.667,       !Telescope Motor Temperature
     *          -0.2667,  66.667,       !Tilt Base Temperature
     *          -0.2667,  66.667,       !Tilt Platform Temperature
     *          -0.2667,  66.667,       !Half Angle Motor Temperature
     *          0.02, 0.26,             !Power Supply A Input Current
     *          0.02, 0.26,             !Power Supply B Input Current
     *          0.075, 0.0,             !+15 V Analog Power Voltage
     *          -0.075, 0.0,            !-15 V Analog Power Voltage
     *          0.025, 0.0,             !+5 V Logic Power Voltage
     *          -0.2667,  66.667,       !Power Supply Temperature
     *          -0.2667,  66.667,       !B1/B2 Postamp Temperature
     *          -0.2667,  66.667,       !Servo Drive Temperature
     *          0.15, 0.0,              !+30 V Servo Power Voltage
     *          0.1044, 0.0,            !+21 V Servo Power Voltage
     *          -0.1044, 0.0,           !-21 V Servo Power Voltage
     *          0.025, 0.0,             !+5 V Servo Power Voltage
     *          8.52, -377.,            !Angular Comp. Phase Error
     *          1.44, 0.0,              !Tilt Platform Position
     *          1.44, 0.0,              !Tilt Base Position
     *          0.14, 0.0,              !Heaters Current
     *          0.0024, 0.0,            !Telescope A Motor Current
     *          0.0024, 0.0,            !Telescope B Motor Current
     *          0.0024, 0.0,            !Half Angle A Motor Current
     *          0.0024, 0.0,            !Half Angle B Motor Current
     *          0.01, -1.25,            !Servo A Phase Error
     *          0.01, -1.25,            !Servo B Phase Error
     *          0.016, 0.0,     !Angular Comp. A Motor Current
     *          0.016, 0.0,     !Angular Comp. B Motor Current
     *          16*0.0/         !Spares

c  Discrete data start word, bit and bit length in SOH packet
        data dis_ins/
     *                  5, 1, 1,        !Servo A Selected (1=A)
     *                  5, 2, 1,        !Angular Comp. On (1=on)
     *                  5, 3, 1,        !Servo A Locked (1=on)
     *                  5, 4, 1,        !Servo B Locked (1=on)
     *                  5, 5, 1,        !Timing A Selected (1=A)
     *                  5, 6, 1,        !Tilt A On (1=on)
     *                  5, 7, 1,        !Tilt B On (1=on)
     *                  5, 8, 1,        !Tilt Telemetry On (1=on)
     *                  6, 1, 1,        !Stow On  (1=on)
     *                  6, 2, 1,        !Stow Aligned (1=yes)
     *                  6, 3, 1,        !Heaters Status (1=enables
     *                  6, 4, 1,        !Solar Door Open  (1=open)
     *                  6, 5, 1,        !Analog Power On (1=on)
     *                  6, 6, 1,        !Tilt Platform Limit (1=yes)
     *                  6, 7, 1,        !Tilt Base Limit (1=yes)
     *                  6, 8, 1,        !Tilt Nadir Aligned (1=yes)
     *                  7, 1, 1,        !Tilt Aft Aligned (1=yes)
     *                  7, 2, 1,        !Tilt Forward Aligned (1=yes)
     *                  7, 3, 1,        !Earth Mode Data On (1=Earth)
     *                  7, 4, 1,        !Half Angle Mirror Side  (1=side 2)
     *                  7, 5, 1,        !Image Data Sync  (1=yes)
     *                  7, 6, 1,        !Angular Comp. at Speed  (1=yes)
     *                  30*0    /       !Spares


c  Only need even bytes from input array
        do i=1,44
          ilsb(i) = inst(2*i)
        end do

c  Convert analog data 
        do i=1,32
          call read_analog(insta(i),loc_ins(1,i),con_ins(1,i),ilsb)
        end do
        
        do i=1,22
          call read_discrete(instd(i),dis_ins(1,i),ilsb)
        end do
        
        return
        end
