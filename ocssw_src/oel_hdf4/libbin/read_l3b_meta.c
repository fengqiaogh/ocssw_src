/*------------------------------------------------------------------------------
 File:  getl3b.c

 Contents:
 get_l3b_open   - initializes and retrieves product-level metadata
 and information that is required by getl3b_record
 get_l3b_record - reads Level-3 binned data from nrec bins at a time
 get_l3b_close  - closes the product files; performs clean-up

 Notes:
 o These routines must be called in order as listed above.
 o get_l3b_record() may be called multiple times.
 o I/O specifications document used - v3.2 (17 June 1994)

 Other relevant files:
 seabin.h    - various #defined constants for level 3 binned data
 also #includes hdf.h
 seaproto.h  - prototypes for public level 3 output functions
 seaprotoi.h - prototypes for low-layer level 3 output functions
 l3b_misc.c  - a lower layer of level 3 output functions

 Modification history:
 Programmer     Organization      Date      Description of change
 --------------   ------------    --------    ---------------------
 Lakshmi Kumar    Hughes STX      11/22/93    Original development
 Lakshmi Kumar    Hughes STX      14/02/94    Added arguments bin_syear
 bin_sday, bin_eyear,
 bin_eday to get_l3b_open.
 Lakshmi Kumar    Hughes STX      04/12/94    Made some changes to
 agree with specs sec 7.3.1
 Lakshmi Kumar    Hughes STX      07/01/94    Modified to reflect specs
 v3.0 changes
 Lakshmi Kumar    Hughes STX      11/17/94    Added meta_l3b struct &
 removed station (I/O Specs
 v4.0)
 Lakshmi Kumar    Hughes STX      02/09/95    Added flag_use & eng_q_use
 flds to get_l3b_open call
 (Ref to I/O specs v4.2)
 Lakshmi Kumar    Hughes STX      05/18/95    Changed datatype of global
 attribute "flag_use"
 (Ref. V4.3 I/O Specs)
 Added required code to
 read start_num fld and
 flags_set fld
 Lakshmi Kumar    HITC            05/25/95    Added code to read global
 attribute "flag_names"
 Lakshmi Kumar    HSTX            09/25/95    Added "orbit" output argu-
 ment to get_l3b_open call
 (V4.4 I/O & V2.8 product
 specs)
 Lakshmi Kumar    HSTX            10/27/95    Read "Start Orbit" and
 "End Orbit" global attrs.
 Lakshmi Kumar    HSTX            03/18/96    Added "start_orbit" and
 "infiles" output argumetns
 Lakshmi Kumar    HSTX            03/29/96    Added "meta_l3b" argument
 to 'get_l3b_close' call
 Lakshmi Kumar    HSTX            04/30/96    Removed "meta_l3b" arg
 from get_l3b_close and
 changed all meta_l3b char
 pointers to char arrays.
 This change was done in
 order to avoid Miami code
 change, as Miami's code
 does not use meta_l3b
 struct.
 Lakshmi Kumar    Hughes STX      06/18/96    Changed defn. of MAX to
 MAXVAL inorder to remove
 compile time warning
 Lakshmi Kumar    HSTX            07/26/96    Changed meta_l3b arrays
 sizes to handle larger
 strings.
 Lakshmi Kumar    Hughes STX      10/11/96    Removed arg "orbit" &
 added "end_orb" to
 get_l3b_open call.
 (Ref. V5.0 I/O specs)

 Joel Gales       Futuretech      11/24/99    Modify to handle variable
 products.
------------------------------------------------------------------------------*/

#include <math.h>
#include "netcdf.h"
#include "meta_l3b.h"
#include <sensorInfo.h>
#include <timeutils.h>
#include <genutils.h>

#include <hdf.h>
#include <mfhdf.h>

/*-----------------------------------------------------------------------------
 Function:  read_attrs

 Returns:   int32_t (status)
 The return code is a negative value if any error occurs, otherwise,
 returns 0.

 Description:
 The function read_attrs reads all the global attributes from the given
 Bin file

 Parameters: (in calling order)
 Type      Name        I/O   Description
 ----      ----        ---   -----------
 int32_t   sdfid        I    ID req to access HDF SDS interface
 meta_l3bType *meta_l3m O    Product meta data passed to the calling
 function
 Modification history:
 Programmer     Organization   Date      Description of change
 -------------- ------------   --------  ---------------------
 Lakshmi Kumar  Hughes STX     11/23/94  Original development

 Lakshmi Kumar  Hughes STX     02/09/95  Added code to read global attrs -
 flag_use and eng_q_use
 (Ref to I/O specs v4.2)
 Lakshmi Kumar  Hughes STX        05/09/95  Changed datatype of global attr.
 "flag_use" (ref. V4.3 I/O Specs)
 Lakshmi Kumar  Hughes STX   07/26/96  Added code to truncate the strings,
 if the strings read are larger than
 the space allocated in meta_l3b
 structure
 Joel Gales     Futuretech    03/14/00 Read "Units" metadata attribute
 Joel Gales     Futuretech    06/27/11 Remove "Data Center", "Station",
 "Station Lon/Lat"
----------------------------------------------------------------------------*/
int32_t read_l3b_meta_hdf4(int32_t sdfid, meta_l3bType *meta_l3b) {
    // first zero out the metadata
    bzero(meta_l3b, sizeof (meta_l3bType));

    int32_t nt, count;
    char *attr_buf;
    if ((attr_buf = (char *) calloc(LG_ATTRSZ + 1, sizeof (char))) == NULL) {
        printf("%s -Error: Cannot allocate memory for attribute buffer\n",
                __FILE__);
        exit(EXIT_FAILURE);
    }

    char attr_name[128];
    int32_t status;

    // fill in missing param values
    strcpy(meta_l3b->data_center, "Unknown");
    strcpy(meta_l3b->mission_char, "Unknown");
    strcpy(meta_l3b->sensor_char, "Unknown");
    strcpy(meta_l3b->station, "Unknown");
    meta_l3b->station_lat = 0.0;
    meta_l3b->station_lon = 0.0;
    strcpy(meta_l3b->replace, "Unknown");
    meta_l3b->orbit = -1;
    meta_l3b->resolution = BAD_FLT;
    strcpy(meta_l3b->binning_scheme, "Integerized Sinusoidal Grid");

    /****  read global attribute "product_name"  */

    status = getattrsz(sdfid, "Product Name", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "product_name", &nt, &count);
        if (status < 0)
            return FAIL;
        else
            strcpy(attr_name, "product_name");
    } else {
        strcpy(attr_name, "Product Name");
    }

    if (count > SM_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b ");
        printf("\nstructure.  Attr size read = %d\n", count);
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        memcpy(meta_l3b->product_name, attr_buf, SM_ATTRSZ - 1);
        meta_l3b->product_name[SM_ATTRSZ - 1] = 0;
    } else {
        if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->product_name)) < 0)
            return FAIL;
        meta_l3b->product_name[count] = 0;
    }

    /**** read global attribute "title" */
    status = getattrsz(sdfid, "Title", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "title", &nt, &count);
        if (status < 0)
            return FAIL;
        else
            strcpy(attr_name, "title");
    } else {
        strcpy(attr_name, "Title");
    }

    if (count > SM_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b ");
        printf("\nstructure.  Attr size read = %d\n", count);
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        memcpy(meta_l3b->title, attr_buf, SM_ATTRSZ - 1);
        meta_l3b->title[SM_ATTRSZ - 1] = 0;
    } else {
        if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->title)) < 0)
            return FAIL;
        meta_l3b->title[count] = 0;
    }

    /**** read global attribute "instrument" */
    status = getattrsz(sdfid, "Sensor Name", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "instrument", &nt, &count);
        if (status < 0)
            return FAIL;
        else
            strcpy(attr_name, "instrument");
    } else {
        strcpy(attr_name, "Sensor Name");
    }

    if (count > SM_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b ");
        printf("\nstructure.  Attr size read = %d\n", count);
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        memcpy(meta_l3b->sensor_name, attr_buf, SM_ATTRSZ - 1);
        meta_l3b->sensor_name[SM_ATTRSZ - 1] = 0;
    } else {
        if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->sensor_name)) < 0)
            return FAIL;
        meta_l3b->sensor_name[count] = 0;
    }

    // need to fill in the sensorID and sensor
    meta_l3b->sensorID = sensorName2SensorId(meta_l3b->sensor_name);
    if (meta_l3b->sensorID != -1)
        strcpy(meta_l3b->sensor, sensorId2InstrumentName(meta_l3b->sensorID));
    else
        strcpy(meta_l3b->sensor, "Unknown");

    /**** read global attribute "units" */
    status = getattrsz(sdfid, "Units", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "units", &nt, &count);
        if (status < 0)
            return FAIL;
        else
            strcpy(attr_name, "units");
    } else {
        strcpy(attr_name, "Units");
    }

    if (count > MD_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b ");
        printf("\nstructure.  Attr size read = %d\n", count);
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        memcpy(meta_l3b->units, attr_buf, MD_ATTRSZ - 1);
        meta_l3b->units[MD_ATTRSZ - 1] = 0;
    } else {
        if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->units)) < 0)
            return FAIL;
        meta_l3b->units[count] = 0;
    }

    /**** read global attribute "platform" */
    status = getattrsz(sdfid, "Mission", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "platform", &nt, &count);
        if (status < 0) {
            count = -1;
        } else {
            strcpy(attr_name, "platform");
        }
    } else {
        strcpy(attr_name, "Mission");
    }

    if (count > SM_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b ");
        printf("\nstructure.  Attr size read = %d\n", count);
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        memcpy(meta_l3b->mission, attr_buf, SM_ATTRSZ - 1);
        meta_l3b->mission[SM_ATTRSZ - 1] = 0;
    } else {
        if (count == -1) {
            strcpy(meta_l3b->mission, "Unknown");
        } else {
            if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->mission)) < 0)
                return FAIL;
            meta_l3b->mission[count] = 0;
        }
    }

    /**** read global attribute "temporal_range" */
    status = getattrsz(sdfid, "Product Type", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "temporal_range", &nt, &count);
        if (status < 0)
            return FAIL;
        else
            strcpy(attr_name, "temporal_range");
    } else {
        strcpy(attr_name, "Product Type");
    }

    if (count > SM_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b ");
        printf("\nstructure.  Attr size read = %d\n", count);
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        memcpy(meta_l3b->prod_type, attr_buf, SM_ATTRSZ - 1);
        meta_l3b->prod_type[SM_ATTRSZ - 1] = 0;
    } else {
        if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->prod_type)) < 0)
            return FAIL;
        meta_l3b->prod_type[count] = 0;
    }

    /**** read global attribute "processing_version" */
    status = getattrsz(sdfid, "Processing Version", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "processing_version", &nt, &count);
        if (status < 0) {
            count = -1;
        } else {
            strcpy(attr_name, "processing_version");
        }

    } else {
        strcpy(attr_name, "Processing Version");
    }

    if (count > SM_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b ");
        printf("\nstructure.  Attr size read = %d\n", count);
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        memcpy(meta_l3b->pversion, attr_buf, SM_ATTRSZ - 1);
        meta_l3b->pversion[SM_ATTRSZ - 1] = 0;
    } else {
        if (count == -1) {
            strcpy(meta_l3b->pversion, "Unknown");
        } else {
            if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->pversion)) < 0)
                return FAIL;
            meta_l3b->pversion[count] = 0;
        }
    }

    /**** read global attribute "date_created" */
    status = getattrsz(sdfid, "Processing Time", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "date_created", &nt, &count);
        if (status < 0)
            return FAIL;
        else
            strcpy(attr_name, "date_created");
    } else {
        strcpy(attr_name, "Processing Time");
    }

    if (count > SM_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b ");
        printf("\nstructure.  Attr size read = %d\n", count);
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        memcpy(meta_l3b->ptime, attr_buf, SM_ATTRSZ - 1);
        meta_l3b->ptime[SM_ATTRSZ - 1] = 0;
    } else {
        if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->ptime)) < 0)
            return FAIL;
        meta_l3b->ptime[count] = 0;
    }

    /**** read global attribute "history" */
    status = getattrsz(sdfid, "Processing Control", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "history", &nt, &count);
        if (status < 0)
            return FAIL;
        else
            strcpy(attr_name, "history");
    } else {
        strcpy(attr_name, "Processing Control");
    }

    if (count > MD_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b ");
        printf("\nstructure.  Attr size read = %d\n", count);
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        memcpy(meta_l3b->proc_con, attr_buf, MD_ATTRSZ - 1);
        meta_l3b->proc_con[MD_ATTRSZ - 1] = 0;
    } else {
        if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->proc_con)) < 0)
            return FAIL;
        meta_l3b->proc_con[count] = 0;
    }

    /**** read global attribute "l2_flag_names" */
    status = getattrsz(sdfid, "L2 Flag Names", &nt, &count);
    if (status >= 0) {
        strcpy(attr_name, "L2 Flag Names");

        if (count > SM_ATTRSZ) {
            printf("\n****read_attrs: Attribute '%s'", attr_name);
            printf(" size is greater than it is defined in meta_l3b ");
            printf("\nstructure.  Attr size read = %d\n", count);
            if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
                return FAIL;
            memcpy(meta_l3b->flag_names, attr_buf, SM_ATTRSZ - 1);
            meta_l3b->flag_names[SM_ATTRSZ - 1] = 0;
        } else {
            if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->flag_names)) < 0)
                return FAIL;
            meta_l3b->flag_names[count] = 0;
        }
    }

    /**** read global attribute "time_coverage_start" */
    status = getattrsz(sdfid, "Start Time", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "time_coverage_start", &nt, &count);
        if (status < 0)
            return FAIL;
        else
            strcpy(attr_name, "time_coverage_start");
    } else {
        strcpy(attr_name, "Start Time");
    }

    if (count > SM_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b structure.");
        return FAIL;
    } else {
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        if (strstr(attr_buf, "T") != NULL) {
            meta_l3b->startTime = isodate2unix(attr_buf);
        } else {
            meta_l3b->startTime = zulu2unix(attr_buf);
        }
    }

    /**** read global attribute "time_coverage_end" */
    status = getattrsz(sdfid, "End Time", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "time_coverage_end", &nt, &count);
        if (status < 0)
            return FAIL;
        else
            strcpy(attr_name, "time_coverage_end");
    } else {
        strcpy(attr_name, "End Time");
    }

    if (count > SM_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b structure.");
        return FAIL;
    } else {
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        if (strstr(attr_buf, "T") != NULL) {
            meta_l3b->endTime = isodate2unix(attr_buf);
        } else {
            meta_l3b->endTime = zulu2unix(attr_buf);
        }
    }

    /**** read global attribute "start_orbit_number" */
    if ((rdattr(sdfid, "Start Orbit", (VOIDP *) & meta_l3b->start_orb)) < 0)
        if ((rdattr(sdfid, "start_orbit_number", (VOIDP *) & meta_l3b->start_orb))
                < 0)
            return FAIL;

    /**** read global attribute "end_orbit_number" */
    if ((rdattr(sdfid, "End Orbit", (VOIDP *) & meta_l3b->end_orb)) < 0)
        if ((rdattr(sdfid, "end_orbit_number", (VOIDP *) & meta_l3b->end_orb))
                < 0)
            return FAIL;

    /**** read global attribute "geospatial_lat_units" */
    status = getattrsz(sdfid, "Latitude Units", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "geospatial_lat_units", &nt, &count);
        if (status < 0)
            return FAIL;
        else
            strcpy(attr_name, "geospatial_lat_units");
    } else {
        strcpy(attr_name, "Latitude Units");
    }

    if (count > SM_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b ");
        printf("\nstructure.  Attr size read = %d\n", count);
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        memcpy(meta_l3b->lat_units, attr_buf, SM_ATTRSZ - 1);
        meta_l3b->lat_units[SM_ATTRSZ - 1] = 0;
    } else {
        if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->lat_units)) < 0)
            return FAIL;
        meta_l3b->lat_units[count] = 0;
    }

    /**** read global attribute "geospatial_lon_units" */
    status = getattrsz(sdfid, "Longitude Units", &nt, &count);
    if (status < 0) {
        status = getattrsz(sdfid, "geospatial_lon_units", &nt, &count);
        if (status < 0)
            return FAIL;
        else
            strcpy(attr_name, "geospatial_lon_units");
    } else {
        strcpy(attr_name, "Longitude Units");
    }

    if (count > SM_ATTRSZ) {
        printf("\n****read_attrs: Attribute '%s'", attr_name);
        printf(" size is greater than it is defined in meta_l3b ");
        printf("\nstructure.  Attr size read = %d\n", count);
        if ((rdattr(sdfid, attr_name, (VOIDP *) attr_buf)) < 0)
            return FAIL;
        memcpy(meta_l3b->lon_units, attr_buf, SM_ATTRSZ - 1);
        meta_l3b->lon_units[SM_ATTRSZ - 1] = 0;
    } else {
        if ((rdattr(sdfid, attr_name, (VOIDP *) meta_l3b->lon_units)) < 0)
            return FAIL;
        meta_l3b->lon_units[count] = 0;
    }

    /**** read global attribute "data_bins" */
    int32 tmpInt;
    if ((rdattr(sdfid, "Data Bins", (VOIDP *) &tmpInt)) < 0)
        if ((rdattr(sdfid, "data_bins", (VOIDP *) &tmpInt)) < 0)
            return FAIL;
    meta_l3b->data_bins = tmpInt;
    
    /**** read global attribute "percent_data_bins" */
    if ((rdattr(sdfid, "Percent Data Bins", (VOIDP *) & meta_l3b->pct_databins))
            < 0)
        if ((rdattr(sdfid, "percent_data_bins",
                (VOIDP *) & meta_l3b->pct_databins)) < 0)
            return FAIL;

    /**** read global attribute "software_name" */
    if ((getattrsz(sdfid, "Software Name", &nt, &count)) >= 0) {
        if (count > SM_ATTRSZ) {
            printf("\n****read_attrs: Attribute '%s'", "Software Name");
            printf(" size is greater than it is defined in meta_l3b ");
            printf("\nstructure.  Attr size read = %d\n", count);
            if ((rdattr(sdfid, "Software Name", (VOIDP *) attr_buf)) < 0)
                return FAIL;
            memcpy(meta_l3b->soft_name, attr_buf, SM_ATTRSZ - 1);
            meta_l3b->soft_name[SM_ATTRSZ - 1] = 0;
        } else {
            if ((rdattr(sdfid, "Software Name", (VOIDP *) meta_l3b->soft_name))
                    < 0)
                return FAIL;
            meta_l3b->soft_name[count] = 0;
        }
    }

    /**** read global attribute "software_version" */
    if ((getattrsz(sdfid, "Software Version", &nt, &count)) >= 0) {
        if (count > SM_ATTRSZ) {
            printf("\n****read_attrs: Attribute '%s'", "Software Version");
            printf(" size is greater than it is defined in meta_l3b ");
            printf("\nstructure.  Attr size read = %d\n", count);
            if ((rdattr(sdfid, "Software Version", (VOIDP *) attr_buf)) < 0)
                return FAIL;
            memcpy(meta_l3b->soft_ver, attr_buf, SM_ATTRSZ - 1);
            meta_l3b->soft_ver[SM_ATTRSZ - 1] = 0;
        } else {
            if ((rdattr(sdfid, "Software Version", (VOIDP *) meta_l3b->soft_ver))
                    < 0)
                return FAIL;
            meta_l3b->soft_ver[count] = 0;
        }
    }

    /**** read global attribute "Input Files" */
    if ((getattrsz(sdfid, "Input Files", &nt, &count)) >= 0) {
        if (count > LG_ATTRSZ) {
            printf("\n****read_attrs: Attribute '%s'", "Input Files");
            printf(" size is greater than it is defined in meta_l3b ");
            printf("\nstructure.  Attr size read = %d\n", count);
            if ((rdattr(sdfid, "Input Files", (VOIDP *) attr_buf)) < 0)
                return FAIL;
            memcpy(meta_l3b->infiles, attr_buf, LG_ATTRSZ - 1);
            meta_l3b->infiles[LG_ATTRSZ - 1] = 0;
        } else {
            if ((rdattr(sdfid, "Input Files", (VOIDP *) meta_l3b->infiles)) < 0)
                return FAIL;
            meta_l3b->infiles[count] = 0;
        }
    }

    /**** read global attribute "input_parameters" */
    if ((getattrsz(sdfid, "Input Parameters", &nt, &count)) >= 0) {
        if (count > LG_ATTRSZ) {
            printf("\n****read_attrs: Attribute '%s'", "Input Parameters");
            printf(" size is greater than it is defined in meta_l3b ");
            printf("\nstructure.  Attr size read = %d\n", count);
            if ((rdattr(sdfid, "Input Parameters", (VOIDP *) attr_buf)) < 0)
                return FAIL;
            memcpy(meta_l3b->input_parms, attr_buf, LG_ATTRSZ - 1);
            meta_l3b->input_parms[LG_ATTRSZ - 1] = 0;
        } else {
            if ((rdattr(sdfid, "Input Parameters",
                    (VOIDP *) meta_l3b->input_parms)) < 0)
                return FAIL;
            meta_l3b->input_parms[count] = 0;
        }
    }

    /**** read global attribute "Northernmost Latitude" */
    if ((getattrsz(sdfid, "Northernmost Latitude", &nt, &count) == 0)) {
        if ((rdattr(sdfid, "Northernmost Latitude", (VOIDP *) & meta_l3b->north)) < 0)
            return FAIL;
    }
    /**** read global attribute "Southernmost Latitud" */
    if ((getattrsz(sdfid, "Southernmost Latitude", &nt, &count) == 0)) {
        if ((rdattr(sdfid, "Southernmost Latitude", (VOIDP *) & meta_l3b->south)) < 0)
            return FAIL;
    }
    /**** read global attribute "Easternmost Longitude" */
    if ((getattrsz(sdfid, "Easternmost Longitude", &nt, &count) == 0)) {
        if ((rdattr(sdfid, "Easternmost Longitude", (VOIDP *) & meta_l3b->east)) < 0)
            return FAIL;
    }
    /**** read global attribute "Westernmost Longitude" */
    if ((getattrsz(sdfid, "Westernmost Longitude", &nt, &count) == 0)) {
        if ((rdattr(sdfid, "Westernmost Longitude", (VOIDP *) & meta_l3b->west)) < 0)
            return FAIL;
    }

    /**** read global attribute "Bin Resolution" */
    if ((getattrsz(sdfid, "Bin Resolution", &nt, &count)) >= 0) {
            if ((rdattr(sdfid, "Bin Resolution",
                (VOIDP *) attr_buf)) < 0)
                return FAIL;
        meta_l3b->resolution = str2resolution(attr_buf);
    }

    /**** read global attribute "Bin Scheme" */
    if ((getattrsz(sdfid, "Binning Scheme", &nt, &count)) >= 0) {
        if (count > SM_ATTRSZ) {
            printf("\n****read_attrs: Attribute '%s'", "Binning Scheme");
            printf(" size is greater than it is defined in meta_l3b ");
            printf("\nstructure.  Attr size read = %d\n", count);
            if ((rdattr(sdfid, "Binning Scheme", (VOIDP *) attr_buf)) < 0)
                return FAIL;
            memcpy(meta_l3b->binning_scheme, attr_buf, SM_ATTRSZ - 1);
            meta_l3b->binning_scheme[SM_ATTRSZ - 1] = 0;
        } else {
            if ((rdattr(sdfid, "Binning Scheme",
                    (VOIDP *) meta_l3b->binning_scheme)) < 0)
                return FAIL;
            meta_l3b->binning_scheme[count] = 0;
        }
    }

    free(attr_buf);
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
 Function:  rdattr

 Returns:   int32_t (status)
 The return code is a negative value if any error occurs, otherwise,
 returns 0.

 Description:
 The function rdattr reads the requested global attribute

 Parameters: (in calling order)
 Type      Name        I/O   Description
 ----      ----        ---   -----------
 int32_t   sdfid        I    ID req to access HDF SDS interface
 char  *   attr_name    I    attribute name
 void  *   buf         I/O   pointer to data buffer

 Modification history:
 Programmer     Organization   Date      Description of change
 -------------- ------------   --------  ---------------------
 Lakshmi Kumar  Hughes STX     11/07/94  Original development
 Joel Gales     Futuretech     06/24/13  Add support for NETCDF4

----------------------------------------------------------------------------*/

int32_t rdattr(int32_t sdfid, char *attr_name, void *buf) {
    int32_t attrnum;
    int status;

    if (sdfid > 0) {
        attrnum = SDfindattr(sdfid, attr_name);
        if ((SDreadattr(sdfid, attrnum, buf)) < 0)
            return FAIL;
    } else {
        status = nc_get_att(-sdfid, NC_GLOBAL, attr_name, buf);
        if (status != NC_NOERR)
            return FAIL;
    }

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
 Function:  getattrsz

 Returns:   int32_t (status)
 The return code is a negative value if any error occurs, otherwise,
 returns 0.

 Description:
 The function getattrsz passes the requested global attribute's
 number type (data type) and the number of values

 Parameters: (in calling order)
 Type      Name        I/O   Description
 ----      ----        ---   -----------
 int32_t   sdfid        I    ID req to access HDF SDS interface
 char  *   attr_name    I    attribute name
 int32_t * nt           O    HDF data type
 int32_t * count        O    number of values in the specified attribute

 Modification history:
 Programmer     Organization   Date      Description of change
 -------------- ------------   --------  ---------------------
 Lakshmi Kumar  Hughes STX     11/07/94  Original development
 Joel Gales     Futuretech     06/24/13  Add support for NETCDF4

----------------------------------------------------------------------------*/

int32_t getattrsz(int32_t id, char *attr_name, int32_t *nt, int32_t *count) {
    int32_t attrnum;
    char name[SM_ATTRSZ];
    int status;
    size_t cnt;

    if (id > 0) {
        attrnum = SDfindattr(id, attr_name);
        if ((SDattrinfo(id, attrnum, name, nt, count)) < 0)
            return FAIL;
    } else {

        status = nc_inq_attlen(-id, NC_GLOBAL, attr_name, &cnt);
        *count = cnt;
        if (status != NC_NOERR)
            return FAIL;
    }
    return SUCCEED;
}

int32_t read_l3b_meta_hdf5(hid_t grp0, meta_l3bType *meta_l3b) {
    // first zero out the metadata
    bzero(meta_l3b, sizeof (meta_l3bType));

    hid_t attr;
    hid_t atype = H5Tcopy(H5T_C_S1);
    char *attr_buf;
    if ((attr_buf = (char *) calloc(MAX_ORDER + 1, sizeof (char))) == NULL) {
        printf("%s -Error: Cannot allocate memory for attribute buffer\n",
                __FILE__);
        exit(EXIT_FAILURE);
    }

    // fill in missing values
    strcpy(meta_l3b->data_center, "Unknown");
    strcpy(meta_l3b->sensor_char, "Unknown");
    strcpy(meta_l3b->station, "Unknown");
    meta_l3b->station_lat = 0.0;
    meta_l3b->station_lon = 0.0;
    strcpy(meta_l3b->units, "Unknown");
    strcpy(meta_l3b->prod_type, "Unknown");
    strcpy(meta_l3b->replace, "Unknown");
    strcpy(meta_l3b->input_parms, "Unknown");
    strcpy(meta_l3b->flag_names, "Unknown");
    strcpy(meta_l3b->infiles, "Unknown");
    meta_l3b->orbit = -1;
    meta_l3b->start_orb = -1;
    meta_l3b->end_orb = -1;
    strcpy(meta_l3b->binning_scheme, "Unknown");
    meta_l3b->resolution = BAD_FLT;

    attr = H5Aopen_name(grp0, "Product Name");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->product_name);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Title");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->title);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Sensor");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->sensor_name);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Sensor Characteristics");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->sensor_char);
    H5Aclose(attr);

    // need to fill in the sensorID and sensor
    meta_l3b->sensorID = sensorName2SensorId(meta_l3b->sensor_name);
    if (meta_l3b->sensorID != -1)
        strcpy(meta_l3b->sensor, sensorId2InstrumentName(meta_l3b->sensorID));
    else
        strcpy(meta_l3b->sensor, "Unknown");

    attr = H5Aopen_name(grp0, "Mission");
    if (attr == -1)
        attr = H5Aopen_name(grp0, "mission");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->mission);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Mission Characteristics");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->mission_char);
    H5Aclose(attr);

    meta_l3b->sensorID = sensorName2SensorId(meta_l3b->sensor);

    attr = H5Aopen_name(grp0, "Processing Version");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->pversion);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Software Name");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->soft_name);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Software ID");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->soft_ver);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Processing Time");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->ptime);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Processing Control");
    H5Tset_size(atype, MD_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->proc_con);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "L2 Flag Names");
    H5Tset_size(atype, MD_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->flag_names);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Start Orbit");
    H5Aread(attr, H5T_STD_U32LE, &meta_l3b->start_orb);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "End Orbit");
    H5Aread(attr, H5T_STD_U32LE, &meta_l3b->end_orb);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Start Time");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, attr_buf);
    H5Aclose(attr);
    meta_l3b->startTime = zulu2unix(attr_buf);

    attr = H5Aopen_name(grp0, "End Time");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, attr_buf);
    H5Aclose(attr);
    meta_l3b->endTime = zulu2unix(attr_buf);

    attr = H5Aopen_name(grp0, "Latitude Units");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->lat_units);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Longitude Units");
    H5Tset_size(atype, SM_ATTRSZ);
    H5Aread(attr, atype, &meta_l3b->lon_units);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Data Bins");
    int32 tmpInt;
    H5Aread(attr, H5T_STD_U32LE, &tmpInt);
    H5Aclose(attr);
    meta_l3b->data_bins = tmpInt;
    
    attr = H5Aopen_name(grp0, "Percent Data Bins");
    H5Aread(attr, H5T_NATIVE_FLOAT, &meta_l3b->pct_databins);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Northernmost Latitude");
    H5Aread(attr, H5T_NATIVE_FLOAT, &meta_l3b->north);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Southernmost Latitude");
    H5Aread(attr, H5T_NATIVE_FLOAT, &meta_l3b->south);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Easternmost Longitude");
    H5Aread(attr, H5T_NATIVE_FLOAT, &meta_l3b->east);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Westernmost Longitude");
    H5Aread(attr, H5T_NATIVE_FLOAT, &meta_l3b->west);
    H5Aclose(attr);

    H5Tclose(atype);
    free(attr_buf);
    return SUCCEED;
}

int32_t read_l3b_meta_netcdf4(int ncid, meta_l3bType *meta_l3b) {
    // first zero out the metadata
    bzero(meta_l3b, sizeof (meta_l3bType));

    char *attr_buf;
    if ((attr_buf = (char *) calloc(MAX_ORDER + 1, sizeof (char))) == NULL) {
        printf("%s -Error: Cannot allocate memory for attribute buffer\n",
                __FILE__);
        exit(EXIT_FAILURE);
    }

    // fill in missing values
    strcpy(meta_l3b->data_center, "Unknown");
    strcpy(meta_l3b->mission_char, "Unknown");
    strcpy(meta_l3b->sensor_char, "Unknown");
    strcpy(meta_l3b->station, "Unknown");
    meta_l3b->station_lat = 0.0;
    meta_l3b->station_lon = 0.0;
    strcpy(meta_l3b->replace, "Unknown");
    meta_l3b->orbit = -1;
    meta_l3b->resolution = BAD_FLT;

    nc_get_att(ncid, NC_GLOBAL, "product_name", &meta_l3b->product_name);
    nc_get_att(ncid, NC_GLOBAL, "title", &meta_l3b->title);
    nc_get_att(ncid, NC_GLOBAL, "instrument", &meta_l3b->sensor);
    nc_get_att(ncid, NC_GLOBAL, "platform", &meta_l3b->mission);

    // create sensor_name
    meta_l3b->sensorID = instrumentPlatform2SensorId(meta_l3b->sensor,
            meta_l3b->mission);
    if (meta_l3b->sensorID == -1)
        meta_l3b->sensorID = OCRVC;
    strcpy(meta_l3b->sensor_name, sensorId2SensorName(meta_l3b->sensorID));

    nc_get_att(ncid, NC_GLOBAL, "units", &meta_l3b->units);
    nc_get_att(ncid, NC_GLOBAL, "temporal_range", &meta_l3b->prod_type);
    nc_get_att(ncid, NC_GLOBAL, "processing_version", &meta_l3b->pversion);
    nc_get_att(ncid, NC_GLOBAL, "date_created", &meta_l3b->ptime);
    nc_get_att(ncid, NC_GLOBAL, "history", &meta_l3b->proc_con);
    nc_get_att(ncid, NC_GLOBAL, "time_coverage_start", attr_buf);
    meta_l3b->startTime = isodate2unix(attr_buf);
    nc_get_att(ncid, NC_GLOBAL, "time_coverage_end", attr_buf);
    meta_l3b->endTime = isodate2unix(attr_buf);

    nc_get_att(ncid, NC_GLOBAL, "start_orbit_number", &meta_l3b->start_orb);
    nc_get_att(ncid, NC_GLOBAL, "end_orbit_number", &meta_l3b->end_orb);

    nc_get_att(ncid, NC_GLOBAL, "northernmost_latitude", &meta_l3b->north);
    nc_get_att(ncid, NC_GLOBAL, "southernmost_latitude", &meta_l3b->south);
    nc_get_att(ncid, NC_GLOBAL, "easternmost_longitude", &meta_l3b->east);
    nc_get_att(ncid, NC_GLOBAL, "westernmost_longitude", &meta_l3b->west);
    nc_get_att_long(ncid, NC_GLOBAL, "data_bins", (long*)&meta_l3b->data_bins);
    nc_get_att(ncid, NC_GLOBAL, "percent_data_bins", &meta_l3b->pct_databins);
    nc_get_att(ncid, NC_GLOBAL, "binning_scheme", &meta_l3b->binning_scheme);

    nc_type attType;
    nc_inq_atttype(ncid, NC_GLOBAL, "geospatial_lon_resolution", &attType);
    if(attType==NC_DOUBLE) {
        double val;
        nc_get_att(ncid, NC_GLOBAL, "geospatial_lon_resolution", &val);
        meta_l3b->resolution = degrees2resolution(val);
    } else {
        nc_get_att(ncid, NC_GLOBAL, "geospatial_lon_resolution", attr_buf);
        meta_l3b->resolution = str2resolution(attr_buf);
    }

    nc_get_att(ncid, NC_GLOBAL, "geospatial_lat_units", &meta_l3b->lat_units);
    nc_get_att(ncid, NC_GLOBAL, "geospatial_lon_units", &meta_l3b->lon_units);

    nc_get_att(ncid, NC_GLOBAL, "identifier_product_doi", &meta_l3b->doi);
    nc_get_att(ncid, NC_GLOBAL, "keywords", &meta_l3b->keywords);
    
    int grp_id;
    nc_inq_grp_ncid(ncid, "processing_control", &grp_id);
    nc_get_att(grp_id, NC_GLOBAL, "software_name", &meta_l3b->soft_name);
    nc_get_att(grp_id, NC_GLOBAL, "software_version", &meta_l3b->soft_ver);
    nc_get_att(grp_id, NC_GLOBAL, "source", &meta_l3b->infiles);
    nc_get_att(grp_id, NC_GLOBAL, "l2_flag_names", meta_l3b->flag_names);

    // get input parameters
    int grp_id2;
    int i;
    int numParams = 0;
    char attName[1024];
    meta_l3b->input_parms[0] = 0;
    i = nc_inq_grp_ncid(grp_id, "input_parameters", &grp_id2);
    nc_inq_natts(grp_id2, &numParams);
    for (i = 0; i < numParams; i++) {
        nc_inq_attname(grp_id2, NC_GLOBAL, i, attName);
        nc_get_att(grp_id2, NC_GLOBAL, attName, attr_buf);
        strcat(meta_l3b->input_parms, attName);
        strcat(meta_l3b->input_parms, " = ");
        strcat(meta_l3b->input_parms, attr_buf);
        strcat(meta_l3b->input_parms, "\n");
    }

    free(attr_buf);
    return SUCCEED;
}
