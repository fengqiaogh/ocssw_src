/* --------------------------------------------------------------- */
/* get_depth.c - water depth classification function for MSl12.    */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/* Outputs:                                                        */
/*     depth - water depth (m)                                     */
/*                                                                 */
/* Algorithm Provided By: R. A. Stumph, NOAA                       */
/* Written By: B. A. Franz, SAIC GSC, SIMBIOS Project, 9 July 1999 */
/* Updated:    29 July 2002, BAF                                   */
/*                                                                 */
/* --------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

void get_depth_classification(l2str *l2rec, float depth[]) {
    static float *Fo;
    static int firstCall = 1;

    int32_t ip, ipb;
    int32_t ib, iib;
    float mu0, mu;
    float Rrs555;
    float Rrs490;
    float diff;

    l1str *l1rec = l2rec->l1rec;
    int32_t nbands = l1rec->l1file->nbands;

    /* Need to reconstruct mean Fo or use nominal*/
    if (firstCall) {
        firstCall = 0;
        if ((Fo = (float *) calloc(nbands, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to Fo in get_depth\n");
            exit(FATAL_ERROR);
        }

        for (ib = 0; ib < nbands; ib++) {
            iib = l2rec->bindx[ib];
            if (l1_input->outband_opt >= 2)
                Fo[iib] = l1rec->l1file->Fonom[iib];
            else
                Fo[iib] = l1rec->l1file->Fobar[iib];
        }
    }

    for (ip = 0; ip < l1rec->npix; ip++) {

        ipb = nbands*ip;

        /*                                                         */
        /* Is it land?                                             */
        /* We can't rely on the standard landmask, since reefs may */
        /* get flagged as land if the navigation is off by a pixel.*/
        /* Algorithm provided by F. S. Patt, SAIC GSC.             */
        /*                                                         */
        mu0 = l1rec->csolz[ip];
        mu = l1rec->csenz[ip];
        diff = (l2rec->l1rec->Lt[ipb + 7] - 0.7 * l2rec->l1rec->Lt[ipb + 1]) / (1.0 + mu0 / mu);
        if (diff > -1.5) {
            depth[ip] = -1.0;
            continue;
        }

        /*                                                         */
        /* Is it cloud, or other conditions which result in no     */
        /* Rrs retrieval?                                          */
        /*                                                         */
        if (l2rec->nLw[ipb + 2] <= 0.0 || l2rec->nLw[ipb + 4] <= 0.0) {
            depth[ip] = -2.0;
            continue;
        }

        Rrs490 = l2rec->nLw[ipb + 2] / Fo[2];
        Rrs555 = l2rec->nLw[ipb + 4] / Fo[4];

        depth[ip] = 43.0 * (-0.95 + log(1000.0 * Rrs490) / log(1000.0 * Rrs555));
        if (depth[ip] < 0.0) depth[ip] = 0.0;

    }
}



/*

SeaWiFS Corals.  Here is the preliminary depth algorithm.
Replace the current depth algorithm with the following: 

Calculate the depth as:
depth = 0.43 * {-0.95 +  ln[1000 Rrs(490)] / ln[1000 Rrs(555) }

I am guessing that the coefficients of 1.4 and -0.95 are acceptable
for SeaWiFS.

Use the same binning strategy.
Then

shallow     = NOT land  and depth <=5 m  
mod. shallow = 5- 10 m
mod. deep    = 10-20 m
deep         = 20-30 m
very deep    = > 30 m

The algorithm is insensitive to bottom cover. 

                             --rick

_______________________________________________________________________
Richard P. Stumpf                              phone:  301-713-3028 x173
NOAA National Ocean Service                      fax:  301-713-4388
Center for Coastal Monitoring and Assessment
1305 East-West Highway, N/SCI1 rm 9115          
Silver Spring, MD 20910                        email:  
richard.stumpf@noaa.gov
_______________________________________________________________________

 */
/*
From: richard.stumpf@noaa.gov
Date: Thu, 08 Jul 99 10:15:34 -0500
To: <gene>, <norman@tursiops>
Subject: SeaWiFS Coral Reef algorithm


Norm and Gene:

Algorithm for using SeaWiFS to find potential coral reef environments.
We should keep in mind a paper in Eos or Oceanography and perhaps a RS journal.

Settings are for 8 bit, all the information warranted at present.
If you prefer 16-bit output, set the intertidal, clouds, and land
to large numbers.  I have tried to anticipate hypothetical future
improvements to the algorithm.

Settings                               suggested  color  rgb

0 = missing data                                  black  0,0,0
5 = shallow water (<~5m)                          cyan   0,255,255
20 = medium water (~5m - ~20m)                          blue   0, 0, 255
100 = deep water                                  purple 128,0, 128

201 = intertidal                                   not  used yet
240 = clouds                                      gray   200,200,200
250 = land                                        dk.green  0,180,0

For high pigment, either a similar color to corresponding shallow,
   or the same color:
105 = pigmented bottom, shallow (prelim.)   cyan/green   0,255,127 or 0,255,255
120 = pigmented bottom, mod. shallow (prelim.)dark cyan   0,127,127 or
0,255,255

Here is the basic algorithm (in fake code form):
!------------------------------------------------------------------------------------
image = 0
if (Rrs(555) le 0.0035 and Rrs(555) gt 0.0)) then image = 100  !deep water
if (Rrs(555) gt 0.0035) then image = 20                     !mod. shallow < ~20
m
if (Rrs(670) gt 0.0015) then image = 5                      !shallow water < ~5
m


!At this 3-depth resolution, pigment should have a slight effect on the
results.
!The chlorophyll/pigment estimate inherently increases with shallowing water,
!with a bare sand bottom, this will reach a semi-analytical value of between
1.5
!and 2.0 ug/L.  The settings for the pigmented bottom allow us to easily merge
!or separate from the non-pigmented bottom

if (chl gt 1.5 and image gt 0)  image = image+100  !allows easy management
                                                   !of rudimentary
                                                   !pigmented-bottom algorithm


! use standard cloudmask,  this will probably need to be changed if we need
! to identify intertidal waters, as some pixels will flag out under both
! the cloud mask and land mask.

if (cloudmask) then image = 240           !if cloudmask is set change values

! land mask, note that this may nail some clouds,
!  it must be done after the cloud mask
!  In thinking about this, I believe that Rrs(670) will work, but See NOTE 1
below:

if (Rrs(670) le 0.) then image = 250      ! land


!---------------------------------------------------------------------------------


NOTES:
1. For the landmask, the algorithm I have used for years for AVHRR and TM is:
        Llm = (Lt(670) - Lr(670))/F0(670)  -  (Lt(865) - Lr(865))/F0(865)
        if (Llm le 0) image = 250  !land
   Using (Rrs(670) le 0) should give the same results.

2. For compositing, take composite of all values
where(image gt 0 and image le 200).
Keep in mind that wind resuspension will be a problem, and could cause some
peculiar results. This is the best for now.

My solution to resuspension has been to take averages, then find the minimum
value
of the averages.  Taking minimum values of individual
scenes makes the result highly sensitive to one faulty scene. But if the
SeaWiFS
output is reliable enough, minimum non-cloud values will work.

3.  I have not resolved appropriate values to permit ready identification of
intertidal shallow water.  The obvious method is to identify areas that shift
between
land and shallow water in different scenes.  The problem is that the land mask
and
cloud mask do not discriminate well against each other.  I will think about
this.

Look forward to seeing some results.
--
                                        --rick
_______________________________________________________________________________
Richard P. Stumpf                                  phone:  301-713-3028 x173
NOAA National Ocean Service                           fax:  301-713-4388
Center for Coastal Monitoring and Assessment
1305 East-West Highway, N/SCI1 rm 9109
Silver Spring, MD 20910                                   email:
 richard.stumpf@noaa.gov
_______________________________________________________________________________



Subject: 
       Re: depth index, Florida test image
   Date: 
       Wed, 14 Jul 99 13:36:42 -0500
  From: 
       richard.stumpf@noaa.gov
    To: 
       <franz>
    CC: 
       <norman@tursiops>, <gene@tutuila>




Bryan:

To compensate for the residual atmosphere, Rrs(555) > 0.004 will work. 

However, to compensate for the resuspension over west Florida, including the 
entrainment of sediment into the Florida Current, Rrs(555) > 0.008, for < 20m.  
This isn't too bad for most areas, but some of the Bahamas high chlorophyll areas 
get flagged as deep water.

That value will correspond to 20 m on the Florida shelf. 

This means that the shallow water value will also need to be boosted. 
Rrs(555) is ok for shallow water,  Rrs(555) > 0.015 would describe the shallow water
adequately for this scene, although over-compensating for the Bahamas Banks. 

(For shallow water with Rrs(670), my guess is to double the value to Rrs(670) > 0.003,
but we might as well use Rrs555 for this scene.)

As a bonus, I put a gif version of the nautical chart for the region in. 
Blue is <10 fathoms (almost 20m). 

orca.gsfc.nasa.gov/pub/incoming/chart_11013_1.gif




-- 
                                        --rick
_______________________________________________________________________________
Richard P. Stumpf                                  phone:  301-713-3028 x173
NOAA National Ocean Service                           fax:  301-713-4388
Center for Coastal Monitoring and Assessment
1305 East-West Highway, N/SCI1 rm 9109                
Silver Spring, MD 20910                                   email:
richard.stumpf@noaa.gov
_______________________________________________________________________________

 */
