/*
 * get_avw.c
 *
 * This algorithm returns the weighted harmonic mean of the visible-range (400 – 700 nm) remote sensing reflectance (Rrs) wavelengths,
 * outputting the Apparent Visible Wavelength (AVW) in units of nanometers.
 *  The AVW is an optical water classification index, representing a one-dimensional geophysical metric that is inherently correlated
 *  to Rrs spectral shape (Vandermeulen et al. 2020).
 *
 *
 *Vandermeulen, R. A., Mannino, A., Craig, S.E., Werdell, P.J., 2020: “150 shades of green: Using the full spectrum of remote sensing reflectance
 *   to elucidate color shifts in the ocean,” Remote Sensing of Environment, 247, 111900, https://doi.org/10.1016/j.rse.2020.111900
 *
 *  Created on: Sep 30, 2020
 *      Author: M. Zhang
 */

#include "l12_proto.h"


float avw_cal_hypspectral(float *Rrs, float *wave, int nwave){

    float avw=BAD_FLT;
    int32_t ib;

    float dfirst=0,dlast=0;
    float *Secondderiv=(float *)malloc(nwave*sizeof(float));

    int32_t nwave_out=700-400+1;
    float *wave_out=(float *)malloc(nwave_out*sizeof(float));
    float *Rrs_out =(float *)malloc(nwave_out*sizeof(float));
    double nominator=0., denominator=0.;

    for(ib=0;ib<nwave_out;ib++)
        wave_out[ib]=400+ib;

    dfirst=first_deriv(wave,Rrs,0);
    dlast=first_deriv(wave,Rrs,nwave);

    spline(wave,Rrs,nwave,dfirst,dlast,Secondderiv);

    for(ib=0;ib<nwave_out;ib++)
        splint(wave,Rrs,Secondderiv,nwave,wave_out[ib],&Rrs_out[ib]);

    for(ib=0;ib<nwave_out;ib++){
        nominator+=Rrs_out[ib];
        denominator+=(Rrs_out[ib]/wave_out[ib]);
    }
    avw=nominator/denominator;

    free(wave_out);free(Rrs_out);
    free(Secondderiv);
    return avw;

}


float avw_cal_multispectral(float *Rrs, float *wave, int nwave){

    float avw=BAD_FLT;
    int32_t ib;
    float *avw_coef=input->avw_coef;

    double nominator=0., denominator=0.,temp=0.;

    for(ib=0;ib<nwave; ib++){
        nominator+=Rrs[ib];
        denominator+=(Rrs[ib]/wave[ib]);
    }
    avw=nominator/denominator;

    for(ib=0;ib<6;ib++)
        temp+=avw_coef[ib]*pow(avw,ib);

    avw=temp;

    return avw;

}

void get_avw(l2str *l2rec, float avw[]){
    filehandle* l1file = l2rec->l1rec->l1file;

    int32_t ip,ib;
    int32_t ipb;
    
    int32_t sensorID = l1file->sensorID;
    int32_t nbands = l1file->nbands;
    int32_t npix = l2rec->l1rec->npix;
    static float fsol;
    float *wave = l1file->fwave;
    float *Rrs;
    static float *Rrs_avw, *wave_avw;
    static int firstcall=1;
    int32_t nwave_avw;
    int32_t negative=0;
    static int ifhyper=0;
    static  int ib400=0,ib700=0;

    if(firstcall){
        firstcall=0;
        Rrs_avw =(float *)malloc(nbands*sizeof(float));
        wave_avw=(float *)malloc(nbands*sizeof(float));
        switch (sensorID){
        case HICO:
        case PRISM:
        case OCI:
            ifhyper=1;
            break;
        }
        if(ifhyper){
            ib400=windex(400, wave,nbands);
            ib700=windex(700, wave,nbands);

            if(ib400>0 && wave[ib400]>=400)
                ib400--;
            if(ib700<nbands-1 && wave[ib700]<700)
                ib700++;
        }
        fsol=l2rec->l1rec->fsol;
    }

    /*                                                      */
    /* Compute desired products at each pixel               */
    /*                                                      */
    for (ip = 0; ip < npix; ip++) {
        avw[ip]=BAD_FLT;

        if(l2rec->l1rec->mask[ip] || l2rec->l1rec->solz[ip] >= SOLZNIGHT)
            continue;

        ipb=ip*nbands;
        Rrs=&l2rec->Rrs[ipb];

        nwave_avw=0;
        negative=0;
        if(ifhyper){
            for(ib=ib400;ib<=ib700;ib++){
                if(Rrs[ib]!=BAD_FLT){
                    Rrs_avw [nwave_avw]=Rrs[ib] * (l1file->Fonom[ib]*fsol/l2rec->l1rec->Fo[ib])/l2rec->outband_correction[ipb+ib];
                    wave_avw[nwave_avw]=wave[ib];
                    nwave_avw++;
                    if(Rrs[ib]<0)
                        negative=1;
                }
            }
        }
        else{
            for(ib=0;ib<nbands;ib++){
                if(wave[ib]>=400 && wave[ib]<=700){
                    if(Rrs[ib]!=BAD_FLT){
                        Rrs_avw [nwave_avw]=Rrs[ib] * (l1file->Fonom[ib]*fsol/l2rec->l1rec->Fo[ib])/l2rec->outband_correction[ipb+ib];
                        wave_avw[nwave_avw]=wave[ib];
                        nwave_avw++;
                        if(Rrs[ib]<0)
                            negative=1;
                    }
                }
            }
        }

        if(negative)
            l2rec->l1rec->flags[ip] |=PRODWARN;

        // ensure enough bands to run a spline
        if(nwave_avw < 4) {
            l2rec->l1rec->flags[ip] |=PRODFAIL;
        } else {
            if(ifhyper)
                avw[ip]=avw_cal_hypspectral(Rrs_avw,wave_avw,nwave_avw);
            else
                avw[ip]=avw_cal_multispectral(Rrs_avw,wave_avw,nwave_avw);
        }

    }
}
void get_Rrs_brightness(l2str *l2rec, float Rrs_brightness[]){

    int32_t ip,ib;
    int32_t ipb;

    int32_t nbands = l2rec->l1rec->l1file->nbands;
    int32_t npix = l2rec->l1rec->npix;
    float *wave=l2rec->l1rec->l1file->fwave;
    float *Rrs;
    static float *Rrs_avw, *wave_avw;
    static int firstcall=1;
    int32_t nwave_avw;
    static  int ib400=0,ib700=0;

    if(firstcall){
        firstcall=0;
        Rrs_avw =(float *)malloc(nbands*sizeof(float));
        wave_avw=(float *)malloc(nbands*sizeof(float));

        ib400=windex(400, wave,nbands);
        ib700=windex(700, wave,nbands);
    }

    /*                                                      */
    /* Compute desired products at each pixel               */
    /*                                                      */
    for (ip = 0; ip < npix; ip++) {
        Rrs_brightness[ip]=BAD_FLT;

        if(l2rec->l1rec->mask[ip] || l2rec->l1rec->solz[ip] >= SOLZNIGHT)
            continue;

        ipb=ip*nbands;
        Rrs=&l2rec->Rrs[ipb];
        nwave_avw=0;

        for(ib=ib400;ib<=ib700;ib++){
                if(Rrs[ib]!=BAD_FLT){
                    Rrs_avw [nwave_avw]=Rrs[ib];
                    wave_avw[nwave_avw]=wave[ib];
                    nwave_avw++;
                }
                else if (Rrs[ib]<0)
                    l2rec->l1rec->flags[ip]|=PRODWARN;
        }

        if(nwave_avw<2)
            continue;

        Rrs_brightness[ip]=0;
        for(ib=0;ib<nwave_avw-1;ib++){
            Rrs_brightness[ip]+=(Rrs_avw[ib]+Rrs_avw[ib+1])/2.0*(wave_avw[ib+1]-wave_avw[ib]);
        }
    }
}

void get_lambda_max(l2str *l2rec, float lambda_max[]){

    int32_t ip,ib;
    int32_t ipb;

    int32_t nbands = l2rec->l1rec->l1file->nbands;
    int32_t npix = l2rec->l1rec->npix;
    float *wave=l2rec->l1rec->l1file->fwave;
    float *Rrs;
    float Rrs_temp;
    static int firstcall=1;
    static int ib400, ib700;
    int negative=0;

    if(firstcall){
        ib400=windex(400, wave,nbands);
        ib700=windex(700, wave,nbands);
        firstcall=0;
    }

    /*                                                      */
    /* Compute desired products at each pixel               */
    /*                                                      */
    for (ip = 0; ip < npix; ip++) {
        negative=0;
        lambda_max[ip]=BAD_FLT;

        if(l2rec->l1rec->mask[ip] || l2rec->l1rec->solz[ip] >= SOLZNIGHT)
            continue;

        ipb=ip*nbands;
        Rrs=&l2rec->Rrs[ipb];
        Rrs_temp=Rrs[ib400];

        if(Rrs_temp>=0)
            lambda_max[ip]=wave[ib400];
        else
            negative=1;

        for(ib=ib400+1;ib<=ib700;ib++){
            if(Rrs[ib]<0)
                negative=1;
            if(Rrs[ib]>Rrs_temp){
                Rrs_temp=Rrs[ib];
                lambda_max[ip]=wave[ib];
            }
        }
        if(negative)
            l2rec->l1rec->flags[ip] |=PRODWARN;
    }
}
