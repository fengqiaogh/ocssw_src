/*
 *  Implement the polarization correction for Hawkeye
 */

#include <netcdf>
#include "polcor_hawkeye.h"

#include <math.h>

#include <allocate2d.h>

#include "l12_proto.h"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

// global variables
static bool firstRun = true;
static float** m12;
static float** m13;


extern "C" void polcor_hawkeye(l1str *l1rec, int32_t ip) {
    int numBands = l1rec->l1file->nbands;
    
    if(firstRun) {
        firstRun = false;
        
        int numPixels = l1rec->npix;
        try {
            NcFile ncFile(input->polfile, NcFile::read);
            NcVar ncVar = ncFile.getVar("m12");
            int lutNumBands = ncVar.getDim(0).getSize();
            if(lutNumBands != numBands) {
                printf("-E- %s:%d - Problem reading polfile %s\n", __FILE__, __LINE__, input->polfile);
                printf("            L1 file #bands=%d, LUT m12 #bands=%d\n", numBands, lutNumBands);
                exit(EXIT_FAILURE);
            }
            int lutNumPixels = ncVar.getDim(1).getSize();
            if(lutNumPixels != numPixels) {
                printf("-E- %s:%d - Problem reading polfile %s\n", __FILE__, __LINE__, input->polfile);
                printf("            L1 file #pix=%d, LUT m12 #pix=%d\n", numPixels, lutNumPixels);
                exit(EXIT_FAILURE);
            }
            m12 = allocate2d_float(numBands, numPixels);
            ncVar.getVar(m12[0]);
            
            ncVar = ncFile.getVar("m13");
            lutNumBands = ncVar.getDim(0).getSize();
            if(lutNumBands != numBands) {
                printf("-E- %s:%d - Problem reading polfile %s\n", __FILE__, __LINE__, input->polfile);
                printf("            L1 file #bands=%d, LUT m13 #bands=%d\n", numBands, lutNumBands);
                exit(EXIT_FAILURE);
            }
            lutNumPixels = ncVar.getDim(1).getSize();
            if(lutNumPixels != numPixels) {
                printf("-E- %s:%d - Problem reading polfile %s\n", __FILE__, __LINE__, input->polfile);
                printf("            L1 file #pix=%d, LUT m13 #pix=%d\n", numPixels, lutNumPixels);
                exit(EXIT_FAILURE);
            }
            m13 = allocate2d_float(numBands, numPixels);
            ncVar.getVar(m13[0]);
            
        } catch   (NcException& e) {
            printf("-E- %s:%d - Problem reading polfile %s\n", __FILE__, __LINE__, input->polfile);
            printf("            %s\n", e.what());
            exit(EXIT_FAILURE);
        }  
    } // firstRun

    for(int band=0; band<numBands; band++) {
        int ipb = ip * numBands + band;
        
        double alpha = l1rec->alpha[ip] / OEL_RADEG;
        double L_x = l1rec->Lt[ipb] / l1rec->tg_sol[ipb] / l1rec->tg_sen[ipb];
        double L_qp = l1rec->L_q[ipb] * cos(2 * alpha) + l1rec->L_u[ipb] * sin(2 * alpha);
        double L_up = l1rec->L_u[ipb] * cos(2 * alpha) - l1rec->L_q[ipb] * sin(2 * alpha);

        l1rec->polcor[ipb] = 1.0 / (1.0 - m12[band][ip] * L_qp / L_x - m13[band][ip] * L_up / L_x);
        l1rec->dpol[ipb] = sqrt(pow(l1rec->L_q[ipb], 2.0) + pow(l1rec->L_u[ipb], 2.0)) / L_x;
    } //for bands
    
    
}