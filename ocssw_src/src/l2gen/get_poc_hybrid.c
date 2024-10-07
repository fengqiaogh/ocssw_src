/*
 * get_poc_hybrid.c
 *
 *  Created on: May 30, 2023
 *      Author: mzhang11
 */

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include <jansson.h>


float range_brdi (float x)
{
    return 1-log10(0.9*x-12.5);
}
float range_mbr (float x)
{
    return log10(0.9*x-12.5);
}
float poc_stramski_hybrid(float *Rrs, int32_t sensorID) {
    static int firstCall = 1;
    static int ib[4]={-1,-1,-1,-1};
    static int nwave;
    int i;
    float ratio[3];
    static float virtual_coef_A[2],virtual_coef_B[2],virtual_rat;
    static int ifvirtual;
    int32_t wave[4];

    float poc = BAD_FLT;
    float Rrs510A,Rrs510B,mbr,brdi,poc_mbr,poc_brdi,w_mbr,w_brdi;
    float Rrs_temp[4];
    static float mbr_coef[4],brdi_coef[6];


    if (firstCall) {

        firstCall = 0;
        ifvirtual=0;

        char filename[FILENAME_MAX];
        char *filedir;
        const char *subsensorDir;

        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(EXIT_FAILURE);
        }
        strcpy(filename, filedir);
        strcat(filename, "/");
        strcat(filename, sensorId2SensorDir(sensorID));
        strcat(filename, "/");

        subsensorDir = subsensorId2SubsensorDir(sensorId2SubsensorId(sensorID));
        if(subsensorDir) {
            strcat(filename, subsensorDir);
            strcat(filename, "/");
        }
        strcat(filename, "poc_stramski_hybrid.json");
        json_t *json,*wv;
        json_error_t error;
        const char *line;

        json = json_load_file(filename, 0, &error);
        if (!json){
            printf("the file %s not exist\n",filename);
            exit(1);
        }

      /*  size_t size=json_dumpb(json,NULL,0,0);

        buffer=(char *)malloc(size);
        size=json_dumpb(json,buffer,size,0);
        line=json_dumps(json,0);*/

        wv=json_object_get(json,"poc_wave");
        line=json_string_value(wv);

        i=0;
        while(line){
        wave[i++]=atof(line);
        line=strchr(line,',');
        if(line)
            line++;
        }
        nwave=i;


        wv=json_object_get(json,"poc_mbr_coef");
        line=json_string_value(wv);
        i=0;
        while(line){
            mbr_coef[i++]=atof(line);
            line=strchr(line,',');
            if(line)
                line++;
        }

        wv=json_object_get(json,"poc_brdi_coef");
        line=json_string_value(wv);
        i=0;
        while(line){
            brdi_coef[i++]=atof(line);
            line=strchr(line,',');
            if(line)
                line++;
        }

        for(i=0;i<nwave;i++)
            ib[i]=bindex_get(wave[i]);

        switch(sensorID){
        case MODISA: case MODIST: case VIIRSJ1:case VIIRSN:
            ifvirtual=1;
            break;
        default:
            break;
        }
        if(ifvirtual){
            switch(sensorID){
            case MODISA: case MODIST:
                virtual_coef_A[0]=-0.00008;
                virtual_coef_A[1]= 1.085;
                virtual_coef_B[0]=-0.00041 ;
                virtual_coef_B[1]= 1.104;
                virtual_rat=0.5;
                break;
            case VIIRSJ1:
                virtual_coef_A[0]=-0.0000004 ;
                virtual_coef_A[1]= 1.068;
                virtual_coef_B[0]=-0.00130;
                virtual_coef_B[1]= 1.291;
                virtual_rat=0.69;
                break;
            case VIIRSN:
                virtual_coef_A[0]=-0.000070;
                virtual_coef_A[1]= 1.096 ;
                virtual_coef_B[0]=-0.00094;
                virtual_coef_B[1]= 1.221;
                virtual_rat=0.63;
                break;
            default:
                printf("-E- %s line %d: coefficients to calculate virtual Rrs(510) are not defined for Stramski hybrid POC\n",
                        __FILE__, __LINE__);
                exit(1);
                break;
            }
        }


        if(nwave==3){
            ib[3]=ib[2];
            ib[2]=-1;
        }
        if (ib[0] < 0 || ib[1] < 0|| ib[3] < 0) {
            printf("-E- %s line %d: required bands not available for Stramski hybrid POC\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    for(i=0;i<nwave;i++){
        if(ib[i]>0 && Rrs[ib[i]]<0)
            return BAD_FLT;
    }

    if(nwave==4){
        for(i=0;i<4;i++)
            Rrs_temp[i]=Rrs[ib[i]];
    }else if(nwave==3){
        for(i=0;i<2;i++)
            Rrs_temp[i]=Rrs[ib[i]];
        Rrs_temp[3]=Rrs[ib[3]];
    }

    if(ifvirtual){
        Rrs510A=virtual_coef_A[0]+virtual_coef_A[1]*Rrs_temp[1];
        if(nwave==4)
            Rrs510B=virtual_coef_B[0]+virtual_coef_B[1]*Rrs_temp[2];
        else if(nwave==3)
            Rrs510B=virtual_coef_B[0]+virtual_coef_B[1]*Rrs_temp[3];

        Rrs_temp[2]=virtual_rat*Rrs510A+(1-virtual_rat)*Rrs510B;
    }

    for(i=0;i<3;i++)
        ratio[i]=Rrs_temp[i]/Rrs_temp[3];

    if(ifvirtual){
        if(isnan(ratio[0])||isnan(ratio[1])||isnan(ratio[2]))
            mbr=MAX(MAX(ratio[0],ratio[1]),ratio[2]);
        else if(ratio[2]<1.2 && Rrs_temp[2]>Rrs_temp[1]&&Rrs_temp[2]>Rrs_temp[0])
            mbr=MAX(MAX(ratio[0],ratio[1]),ratio[2]);
        else
            mbr=MAX(ratio[0],ratio[1]);
    }
    else
        mbr=MAX(MAX(ratio[0],ratio[1]),ratio[2]);

    mbr=log10(mbr);
    poc_mbr=0.;
    for(i=0;i<4;i++)
        poc_mbr+=mbr_coef[i]*pow(mbr,i);
    poc_mbr=pow(10.,poc_mbr);

    brdi=(Rrs_temp[0]-Rrs_temp[3])/Rrs_temp[1];
    poc_brdi=0.;
    for(i=0;i<6;i++)
        poc_brdi+=brdi_coef[i]*pow(brdi,i);
    poc_brdi=pow(10.,poc_brdi);

    if(isnan(poc_brdi))
        w_brdi=-1;
    else if(poc_brdi<15.)
        w_brdi=1.0;
    else if(poc_brdi>25.)
        w_brdi=0.0;
    else
        w_brdi=range_brdi(poc_brdi);

    if(isnan(poc_mbr))
        w_mbr=-1;
    else if(poc_mbr<15.)
        w_mbr=0.0;
    else if(poc_mbr>25.)
        w_mbr=1.0;
    else
        w_mbr=range_mbr(poc_mbr);


    w_mbr=(w_mbr+(1-w_brdi))*0.5;
    w_brdi=1-w_mbr;

    if(isnan(brdi))
        poc=poc_mbr;
    else if(brdi<1.)
        poc=poc_mbr;
    else
        poc=poc_brdi*w_brdi+poc_mbr*w_mbr;

    return (poc);
}