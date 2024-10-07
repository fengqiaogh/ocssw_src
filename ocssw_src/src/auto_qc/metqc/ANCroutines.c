#include <mfhdf.h>
#include "ancil.h"
#include "l1io.h"


/*********************************************************************
 * NAME:     startHDF
 *
 * PURPOSE:  open HDF file structure using SDstart and Hopen routines
 *
 * ARGS:
 *  outfile - output filename
 *  sdfid   - SD file id
 *  fid     - HDF file id
 *  mode    - DFACC_RDONLY, DFACC_CREATE, DFACC_RDWR
 *
 * OUTPUTS:  returns a SD structure ID and a File ID.
 *
 * EFFECTS:  opens HDF file structures
 *
 * RETURNS:  status flag
 * Mods: BDS 2/10/94 ver 3.3, (pg: SD32) does not use DFAAC_WRITE so
 *            its removed here as an option arg.
 * 8/21/96 BDS - renamed 'perror' to 'pexit' to avoid HDF4.0 conflict.
 *********************************************************************/

int startHDF(outfile, sdfid, fid, mode)
char *outfile;
int32_t *sdfid, *fid, mode;
{
    int32 lsdfid, lfid;

    if ((lsdfid = SDstart(outfile, mode)) < 0) return FAIL;
    if ((lfid = Hopen(outfile, DFACC_RDONLY, 0)) < 0) return FAIL;
    Vstart(lfid);

    *sdfid = lsdfid;
    *fid = lfid;

    return SUCCESS;

} /* startHDF */

/*********************************************************************
 * setupGrid 
 *
 * create and HDF Vgroup grid and name it with passed arguments
 *
 *********************************************************************/

int32_t setupGrid(fid, grpname)
int32_t fid;
char *grpname;
{
    int32 gridid;

    gridid = Vattach(fid, -1, "w");
    /* Vsetclass(gridid, vgroupclass); */
    Vsetname(gridid, grpname);

    return gridid;
}

/*********************************************************************
 * gridToGrid 
 *
 * attach inner grid to outer grid: nested Vgroups
 *
 *********************************************************************/

int32_t gridToGrid(outergridid, innergridid)
int32_t outergridid, innergridid;
{
    int32 groupref;

    groupref = VQueryref(innergridid);
    if (groupref < 0) return -1;

    Vaddtagref(outergridid, DFTAG_VG, groupref);

    return 0;
}

/*********************************************************************
 * writeGeom
 *
 * write HDF geometry Vdata using SeaWiFS defined geometry structure
 *
 * return geometry ID
 *
 *********************************************************************/

int32_t writeGeom(int32_t fid, int32_t gridid, char *geomname, int32_t bin_meth, int32_t registration,
        float vsize, float hsize, float max_north, float max_south, float max_west, float max_east)
{
    int32 geomid, geomref, *i32;
    int result;
    float *f32;
    unsigned char arrhold[GEOMSIZE];

    /*
     * --- Write geometry struct for data array and QC data array ----
     */

    geomid = VSattach(fid, -1, "w");
    result = VSfdefine(geomid, "bin_meth", DFNT_INT32, 1);
    result = VSfdefine(geomid, "registration", DFNT_INT32, 1);
    result = VSfdefine(geomid, "vsize", DFNT_FLOAT32, 1);
    result = VSfdefine(geomid, "hsize", DFNT_FLOAT32, 1);
    result = VSfdefine(geomid, "max_north", DFNT_FLOAT32, 1);
    result = VSfdefine(geomid, "max_south", DFNT_FLOAT32, 1);
    result = VSfdefine(geomid, "max_west", DFNT_FLOAT32, 1);
    result = VSfdefine(geomid, "max_east", DFNT_FLOAT32, 1);

    VSsetclass(geomid, GEOMCLASS);
    VSsetname(geomid, geomname);

    result = VSsetfields(geomid,
            "bin_meth,registration,vsize,hsize,max_north,max_south,max_west,max_east");

    i32 = (void *) &arrhold[0];
    *i32 = bin_meth;
    i32 = (void *) &arrhold[4];
    *i32 = registration;

    f32 = (void *) &arrhold[8];
    *f32 = vsize;
    f32 = (void *) &arrhold[12];
    *f32 = hsize;
    f32 = (void *) &arrhold[16];
    *f32 = max_north;
    f32 = (void *) &arrhold[20];
    *f32 = max_south;
    f32 = (void *) &arrhold[24];
    *f32 = max_west;
    f32 = (void *) &arrhold[28];
    *f32 = max_east;

    result = VSwrite(geomid, arrhold, 1, FULL_INTERLACE);
    VSdetach(geomid);

    /*
     * get geometry reference and add to vgroup
     */

    geomref = VSQueryref(geomid);
    Vaddtagref(gridid, DFTAG_VH, geomref);

    if (result < 0)
        return ERROR;
    else
        return geomid;

} /* writeGeom */

/*********************************************************************
 * findGeomId
 *
 * find existing HDF geometry ID 
 *
 *********************************************************************/

int32 findGeomId(fid, geomname)
int32 fid;
char *geomname;
{
    int32 geomid, vsid;

    if ((vsid = VSfind(fid, geomname)) < 0)
        pexit("getting geom vsid");

    if ((geomid = VSattach(fid, vsid, "r")) < 0)
        pexit("attaching geom vsid");

    return geomid;

} /* findGeomId */

/*********************************************************************
 * linkGeom
 *
 * link to existing HDF geometry Vdata 
 * get geometry reference and add to vgroup
 *
 *********************************************************************/

int32 linkGeom(gridid, geomid)
int32 gridid, geomid;
{
    int32 geomref, result;

    geomref = VSQueryref(geomid);
    result = Vaddtagref(gridid, DFTAG_VH, geomref);

    return result;

} /* linkGeom */

/*********************************************************************
 * detachGeom
 *
 * detach Geometry ID 
 *
 *********************************************************************/

int32 detachGeom(geomid)
int32 geomid;
{

    VSdetach(geomid);

    return 0;

} /* detachGeom */

/*********************************************************************
 * addAttr
 *
 * add an attribute to an existing SDS
 *
 *********************************************************************/

int addAttr(sdsid, dataattr, datatype, dataunit)
int32 sdsid;
char *dataattr;
int32 datatype;
char *dataunit;
{
    int result;

    if ((result = SDsetattr(sdsid, dataattr, datatype, strlen(dataunit) + 1,
            dataunit)) != 0) return (-1);

    return result;

} /* addAttr */

/*********************************************************************
 * setSDSref
 *
 * get and SDS reference number and add it to an existing Vgroup
 * end SDS access
 *
 *********************************************************************/

int setSDSref(sdsid, gridid)
int32 sdsid, gridid;
{
    int32 sdsref;
    int result;

    if ((sdsref = (int32) SDidtoref(sdsid)) < 0) return (-1);

    Vaddtagref(gridid, DFTAG_NDG, sdsref);
    if ((result = SDendaccess(sdsid)) < 0) return result;

    return 0;

} /* setSDSref */

/*********************************************************************
 * deattachHDFgrid
 *
 * make HDF grid structure inactive
 *
 *********************************************************************/

int deattachHDFgrid(gridid)
int32 gridid;
{
    return Vdetach(gridid);

} /* deattachHDFgrid */

/*********************************************************************
 * closeHDFstructs
 *
 * close HDF file structures
 *
 *********************************************************************/

int closeHDFstructs(sdfid, fid)
int32 sdfid, fid;
{
    int result;

    if ((result = Vend(fid)) < 0) {
        printf("Vend: result: %d\n", result);
        return (-1);
    }

    if ((result = SDend(sdfid)) < 0) {
        printf("SDend: result: %d\n", result);
        return (-1);
    }

    if ((result = Hclose(fid)) < 0) {
        printf("Error on Hclose(); stack follows:\n");
        HEprint(stderr, 0);
        printf("Stack complete\n");
        return (-1);
    }

    return result;

} /* closeHDFstructs */

/***************************************************************** 
 * File:				wrtsds
 *
 * Purpose:			write an SDS file
 *
 * Description:	writes an SDS given labels, formats, units, and a
 *              	data block.	
 *
 * Input parms: 		
 *  char *sdsfid  	- file id
 *  int rank		   - number of dimensions in SDS
 *  int32 shape[] 	- dimension sizes for each dimension
 *  int32 datatype	- HDF datatype of data block
 *  char *datalabel 	- label descriptor for SDS
 *  char *dataunit	- unit descriptor for SDS
 *  char *datafmt 	- format descriptor for SDS
 *  void *data	   	- array of data values
 *
 * Output parms: 	none		
 *
 * Returns:     	0 on success, non-zero on failure		
 *
 * Local vars:	
 *  int result 	- HDF routine return value
 *
 * Subs called:	HDF library routines SDcreate, SDwritedata
 *
 * History:     	none		
 *
 * Note:				uses "ancil.h"
 *
 * Author:	Brian D. Schieber, GSC, 2/93
 *
 * Modification history:	
 *  BDS, 9/3/93	- modified for new HDF interface
 *
 *****************************************************************/

int32 wrtsds(sdfid, rank, shape, datatype, datalabel, data)

int32 sdfid, shape[], datatype;
int rank;
char *datalabel;
void *data;
{
    int result = 0;
    int32 sdsid;
    static int32 start[2] = {0, 0};

    if ((sdsid = SDcreate(sdfid, datalabel, datatype, rank, shape)) < 0)
        return sdsid;

    if ((result = SDwritedata(sdsid, start, NULL, shape, data)) < 0)
        return result;

    return sdsid;

} /* wrtsds */

/***************************************************************** 
 * File:				rewrtsds
 *
 * Purpose:			overwrite an SDS file
 *
 * Description:	overwrite an SDS given prior sdsid, data shape, and
 *              	data block.	
 *
 * Input parms: 		
 *  char *sdsid  		- sds id
 *  int32 shape[] 	- dimension sizes for each dimension
 *  void *data	   	- array of data values
 *
 * Output parms: 	none		
 *
 * Returns:     	0 on success, non-zero on failure		
 *
 * Local vars:	
 *  int result 	- HDF routine return value
 *
 * Subs called:	HDF library routines SDwritedata
 *
 * History:     	none		
 *
 * Note:				uses "ancil.h"
 *
 * Author:	Brian D. Schieber, GSC, 2/93
 *
 * Modification history:	
 *  BDS, 9/3/93	- modified for new HDF interface
 *
 *****************************************************************/

int32 rewrtsds(sdsid, shape, data)
int32 sdsid, shape[];
void *data;
{
    int result = 0;
    static int32 start[2] = {0, 0};

    printf("rewrtsds sdsid %d, shape[0] %d, shape[1] %d\n",
            sdsid, shape[0], shape[1]);

    result = SDwritedata(sdsid, start, NULL, shape, data);
    printf("rewrtsds SDwritedata result: %d\n", result);

    result += SDendaccess(sdsid);
    printf("rewrtsds SDendaccess result: %d\n", result);

    return result;

} /* rewrtsds */

/***************************************************************** 
 * File:	rdsds
 *
 * Purpose:	read HDF Scientific Data Set (SDS) 
 *
 * Description:	read HDF Scientific Data Set (SDS) and return 
 *                descriptors and data values.  Reads one SDS at a time.
 *
 * Input parameters: 		
 *   char    *filename 		- HDF base name (no ext)
 *   char    *vgname 		- name of vgroup where SDS exists
 *   char    *sdsname 		- name of SDS data grid
 *
 * Output parameters: 		
 *   int32    dimsizes[]	- size of each dimension
 *   void    *inData 		- array of data
 *
 * Returns:	-1 on failure, 0 on success	
 *
 * Local variables:	
 *   int result		      - counters, return value
 *
 * Subroutines called:	HDF library routines
 *
 * History:     none
 *
 * Note:	arrays returned from this routine must 'fit' the variables
 *       allocated space in the calling routine.  Any dimension, any size
 *       SDS can be returned.  The 'dimsizes' returned should be verified
 *       to determine the SDS size retreived.
 *
 * Author:	Brian D. Schieber, GSC, 4/93
 *
 * Modification history:	
 *   modified for the new HDF formats 10/8/93, BDS, GSC/SAIC
 *  BDS, 8/14/95 Mods per OAPS spec 2.7
 *****************************************************************/

int rdsds(filename, vgname, sdsname, dimsizes, inData)
char *filename, *sdsname;
char *vgname;
int32 dimsizes[];
void *inData;
{
    int result;
    int j;
    int32 sdfid, fid, vid, vg;
    int32 sdsindex, sdsid;
    int32 rank, numbertype;
    int32 start[2];
    int32 tagarray1[255], refarray1[255], nrefs1, nattrs1;
    char lname[MAXNAMELNG];

    /* 
     * Default settings for SDreaddata
     */

    start[0] = 0;
    start[1] = 0;

    /*
     * ---- Open data file ------------
     */

    if ((result = startHDF(filename, &sdfid, &fid, DFACC_RDONLY)) != 0) {
        pwarning("Fatal error starting HDF file");
        return (ERROR);
    }

    /*
     * Traverse down to the correct SDS
     */

    /* Open the Vgroup by name and fileid, attach to it and get the number 
     * and tag type (SDG, VH, etc) in it, loop over each item and find 
     * the ones that match the desired SDS and then use the SDS reference 
     * and extract the SDS array to the calling procedure. 
     */

    if ((vid = Vfind(fid, vgname)) <= 0) return (ERROR);
    vg = Vattach(fid, vid, "r");
    nrefs1 = Vntagrefs(vg);
    nattrs1 = nrefs1;
    Vgettagrefs(vg, tagarray1, refarray1, nattrs1);

    for (j = 0; j < nrefs1; j++) {
        if ((tagarray1[j] == DFTAG_NDG) ||
                (tagarray1[j] == DFTAG_SD)) { /* SDS found */
            sdsindex = SDreftoindex(sdfid, refarray1[j]);
            sdsid = SDselect(sdfid, sdsindex);
            SDgetinfo(sdsid, lname, &rank, dimsizes, &numbertype, &nattrs1);
            if (!strcmp(sdsname, lname)) {
                if ((SDreaddata(sdsid, start, (int32 *) NULL, dimsizes, inData)) != 0) {
                    pwarning("Fatal error reading SDS data array");
                    return (ERROR);
                }
                Vdetach(vg);
                return (0); /* SUCCESS */
            } /* if (!strcmp(vgname, lname)) */
        } /* tagarray1[k] */
    } /* for j */
    Vdetach(vg);
    pwarning("Fatal error: Case 1 SDS data array not found");
    return (ERROR);

    /*
     * ---- Close data file ------------
     */

    if ((closeHDFstructs(sdfid, fid)) < 0) {
        pwarning("Fatal error closing SDS: closeHDFstructs");
        return (ERROR);
    }

    return (result);

} /* rdsds */

/********************************************************************** 
 * SUBROUTINE: SDSinFile
 *
 * PURPOSE: assign attributes, add to Vgroup and deattach SDS 
 *
 * AUTHOR: B.D.Schieber, GSC, 10/93
 * Mods: BDS 5/19/95 Support SeaWiFS Specs 2.7
 **********************************************************************/

int32 SDSinFile(sdsname, longname, units, datafmt,
        datatype, sdfid, rank, shape, data, gridid)
char *sdsname, *longname, *units, *datafmt;
int32 sdfid, datatype, rank, shape[], gridid;
void *data;
{
    int result;
    int32 sdsid;

    /*
     * write SDS array
     */

    if ((sdsid = wrtsds(sdfid, rank, shape, datatype, sdsname, data)) < 0) {
        pwarning("wrtsds");
        return (ERROR);
    }

    /*
     * set SDS attributes
     */

    if ((result = addAttr(sdsid, "long_name", DFNT_CHAR, longname)) != 0) {
        pwarning("addAttr");
        return (ERROR);
    }

    if (strlen(units) > 0) /* don't write "" attributes */
        if ((result = addAttr(sdsid, "units", DFNT_CHAR, units)) != 0) {
            pwarning("addAttr");
            return (ERROR);
        }

    /*
     * add SDS to Vgroup and deattach SDS
     */

    if ((result = setSDSref(sdsid, gridid)) != 0) {
        pwarning("setSDSref");
        return (ERROR);
    }

    return 0;

} /* SDSinFile */

/***************************************************************** 
 * File:				wrtattr	
 *
 * Purpose:			write DAAC style global attributes to an HDF file
 *
 * Description:	string arrays of label/value pairs written to HDF 
 *              	header file.
 *
 * Input parms: 
 *	dfile         - HDF file id 
 *      num      - number of strings
 *
 * Output parms:none
 *
 * Returns: 	    -1 on failure, 0 on success
 *
 * Local vars:   
 *   int i		                   - counter
 *   char outlabel[MAXLABLEN+1]   - local label string
 *   char outdescr[MAXDESCLEN+1]  - local description string
 *   int32 dfile	                - file id number
 *   int32 result	                - return value for HDF calls 
 *
 * Subs called: HDF library routines
 *
 * Note:	uses header (ancil.h) definitions for label length, 
 *       description length, and name length
 *       uses argument setting for numannarr
 *
 * History:	based on DAAC "Metadata Submission Guide" 2/93
 *
 * Author:	Brian D. Schieber, GSC, 2/93
 *
 * Mod history:
 *  BDS, 8/30/93 	- modified to work with new HDF design.
 *
 *****************************************************************/

int wrtattr(dfile, annot, numannarr)
int32 dfile;
struct annotation *annot;
int numannarr;
{

    /* 
     * local variables 
     */

    char outlabel[MAXLABLEN + 1], outchar[MAXDESCLEN + 1];
    int16 outint16[1];
    int32 outint32[1];
    float outfloat32[1];
    int i = 0;
    int32 result = 0;

    /* 
     * Place each label/description pair into HDF header 
     */

    for (i = 0; i < numannarr; i++) {

        /* printf ("label [%d] [%s]\n",i, annot[i].label); */
        /* printf ("descr [%d] [%s]\n",i, annot[i].descr); */
        strcpy(outlabel, annot[i].label);

        switch (annot[i].type) {
        case (DFNT_CHAR8):
            strcpy(outchar, annot[i].descr);
            result = SDsetattr(dfile, outlabel, DFNT_CHAR8,
                    strlen(outchar) + 1, outchar);
            break;

        case (DFNT_INT16):
            /* strcpy (outchar, annot[i].descr); */
            /* sscanf(&outchar[0], "%hd", &outint16[0]); */
            sscanf(annot[i].descr, "%hd", &outint16[0]);
            /* printf("outint16 [%hd]\n", outint16[0]); */

            result = SDsetattr(dfile, outlabel, DFNT_INT16, 1, outint16);
            break;

        case (DFNT_INT32):
            sscanf(annot[i].descr, "%d", &outint32[0]);
            result = SDsetattr(dfile, outlabel, DFNT_INT32, 1, outint32);
            break;

        case (DFNT_FLOAT32):
            /* strcpy (outchar, annot[i].descr); */
            /* sscanf(&outchar[0], "%f", &outfloat32[0]); */
            sscanf(annot[i].descr, "%f", &outfloat32[0]);
            /* printf("outfloat32 [%f]\n", outfloat32[0]); */
            result = SDsetattr(dfile, outlabel, DFNT_FLOAT32, 1, outfloat32);
            break;

        default:
            printf("Error in function WRTATTR, default CASE encountered\n");
            return (-1);
            break;

        } /* end switch */
        if (result) printf("Error in writing header annotations\n");
    }
    return (result);
} /* wrtattr */

