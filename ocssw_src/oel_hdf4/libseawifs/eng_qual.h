/*
 $Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.2/L012_Util/util/osc/eng_qual.h,v 4.11 1995/01/25 17:24:26 seawifsd Exp seawifsd $
 $Log: eng_qual.h,v $
 Revision 4.11  1995/01/25 17:24:26  seawifsd
 1. added macro definition of USE_POWER_A/USE_POWER_B to indicate the
 choice of power supply.
 2. updated macro definition FILLLIMIT to include the extra fields that was
 used when analog power is off.
 3. defined several macros(LIMITTYPE,LIMITHARDLOW,...) to indicate the
 order of the index in defining limits for each telemetry field.
 4. defined a structure eng_qualbitStruct for the engineer quality flag field
 and a union eng_qual.
 5. included the necessary include header files.

 Revision 4.10  1995/01/18 14:38:23  seawifsd
 defined constants for engineer quality flags.

 */

#ifndef ENG_QUAL_H_
#define ENG_QUAL_H_



#define LFLAG0 0
#define LFLAG1 1
#define LFLAG2 2
#define LFLAG3 4
#define LFLAG4 8

#define USE_POWER_A 1
#define USE_POWER_B 0

#define GRANULE 0.0001
#define FILLLIMIT LFLAG0,-GRANULE,GRANULE,-GRANULE,GRANULE,-GRANULE,GRANULE
#define LIMITTYPE  0
#define LIMITHARDLOW  1
#define LIMITHARDHIGH  2
#define LIMITSOFTLOW  3
#define LIMITSOFTHIGH  4
#define LIMITSOFTOFFLOW  5
#define LIMITSOFTOFFHIGH 6

typedef struct eng_qualbitStruct {
    /* first byte */
    unsigned int BAND_1_2_FPA_TEMPERATURE_FLAG : 1;
    unsigned int BAND_3_4_FPA_TEMPERATURE_FLAG : 1;
    unsigned int BAND_5_6_FPA_TEMPERATURE_FLAG : 1;
    unsigned int BAND_7_8_FPA_TEMPERATURE_FLAG : 1;
    unsigned int TELESCOPE_MOTOR_TEMPERATURE_FLAG : 1;
    unsigned int TILT_BASE_TEMPERATURE_FLAG : 1;
    unsigned int TILT_PLATFORM_TEMPERATURE_FLAG : 1;
    unsigned int HALF_ANG_MOTOR_TEMPERATURE_FLAG : 1;
    /* second byte */
    unsigned int POWER_SUPPLY_A_INPUT_CURRENT_FLAG : 1;
    unsigned int POWER_SUPPLY_B_INPUT_CURRENT_FLAG : 1;
    unsigned int POWER_SUPPLY_POS_15_VOLT_ANALOG_FLAG : 1;
    unsigned int POWER_SUPPLY_NEG_15_VOLT_ANALOG_FLAG : 1;
    unsigned int POWER_SUPPLY_POS_5_VOLT_LOGIC_FLAG : 1;
    unsigned int POWER_SUPPLY_TEMPERATURE_FLAG : 1;
    unsigned int BAND_1_2_POST_AMP_TEMPERATURE_FLAG : 1;
    unsigned int SERVO_DRIVER_TEMPERATURE_FLAG : 1;
    /* thrid byte */
    unsigned int POWER_SUPPLY_POS_30_VOLT_SERVO_FLAG : 1;
    unsigned int POWER_SUPPLY_POS_21_VOLT_SERVO_FLAG : 1;
    unsigned int POWER_SUPPLY_NEG_21_VOLT_SERVO_FLAG : 1;
    unsigned int POWER_SUPPLY_POS_5_VOLT_SERVO_FLAG : 1;
    unsigned int ANG_MOM_COMP_PHASE_ERROR_FLAG : 1;
    unsigned int TILT_PLATFORM_POSITION_FLAG : 1;
    unsigned int TILT_BASE_POSITION_FLAG : 1;
    unsigned int HEATERS_CURRENT_FLAG : 1;
    /* fourth byte */
    unsigned int TELESCOPE_A_MOTOR_CURRENT_FLAG : 1;
    unsigned int TELESCOPE_B_MOTOR_CURRENT_FLAG : 1;
    unsigned int HALF_ANG_MIR_A_MOTOR_CURRENT_FLAG : 1;
    unsigned int HALF_ANG_MIR_B_MOTOR_CURRENT_FLAG : 1;
    unsigned int SERVO_A_PHASE_ERROR_FLAG : 1;
    unsigned int SERVO_B_PHASE_ERROR_FLAG : 1;
    unsigned int ANG_MOM_COMP_A_MOTOR_CURRENT_FLAG : 1;
    unsigned int ANG_MOM_COMP_B_MOTOR_CURRENT_FLAG : 1;
} eng_qualbitType;

#ifndef byte
#define byte unsigned char
#endif /* byte */

typedef union eng_qual {
    byte bt[4];
    short st[2];
    int in;
    eng_qualbitType bits;
} eng_qualType;

#include "stdio.h"
//#include	"usrhdr.h"
#include "usrmac.h"
#include "eng_qual_proto.h"

#endif/* ENG_QUAL_H_ */
