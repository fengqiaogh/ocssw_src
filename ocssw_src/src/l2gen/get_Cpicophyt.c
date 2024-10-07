/*
 * get_Cpicophyt.c
 *
 *This algorithm calculates cell abundances of Prochlorococcus, Synechococcus, and autotrophic picoeukaryotes in surface waters using principal
 *component analysis (PCA) of hyperspectral and multispectral Rrs.
 *
 *
 *Reference:
 *Priscila Kienteca Lange, P. Jeremy Werdell, Zachary K. Erickson, Giorgio Dall’Olmo, Robert J. W. Brewin, Mikhail V. Zubkov,
 *Glen A. Tarran, Heather A. Bouman, Wayne H. Slade, Susanne E. Craig, Nicole J. Poulton, Astrid Bracher, Michael W. Lomas,
 *and Ivona Cetinić, "Radiometric approach for the detection of picophytoplankton assemblages across oceanic fronts," Opt. Express 28, 25682-25705 (2020)
 *
 * Created on: Sep 29, 2023
 *      Author: Minwei Zhang
 */

#include "l12_proto.h"
#include <hdf5.h>
#include <jansson.h>

void rrs_standardize(float *rrs,int nwave)
{
    int i;
    float mean=0., std=0.;

    for(i=0;i<nwave;i++)
        mean+=rrs[i];
    mean/=nwave;

    for(i=0;i<nwave;i++)
        std+=pow(rrs[i]-mean,2);
    std=sqrt(std/(nwave-1));

    for(i=0;i<nwave;i++)
        rrs[i]=(rrs[i]-mean)/std;
}

void get_Cpicophyt(l2str *l2rec, l2prodstr *p, float *cell_abundance)
{
    int32_t i,j,ip,ipb;
    static int firstcall=1;
    static int nbands, npix;
    static float32 *pca_table;
    static double *pro_coef, *syn_coef, *apeuk_coef;
    static int nwave_pca, npc,npc_used;
    static int *wave_pca;
    static double *pc;
    float *sst;
    static int *bindx;
    float *Rrs=l2rec->Rrs;
    static float *interp_rrs;
    static float *wave;
    static float *pcc_all;

    if(firstcall){
        firstcall=0;
        nbands=l2rec->l1rec->l1file->nbands;
        wave=l2rec->l1rec->l1file->fwave;
        npix=l2rec->l1rec->npix;

        //read the pca components from the table
        char filename[FILENAME_MAX];
        char *filedir;

        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(EXIT_FAILURE);
        }
        strcpy(filename, filedir);
        strcat(filename, "/common/");

        strcat(filename, "pca_picophyto.h5");

        printf("Reading pca table for picophyto from: %s\n", filename);

        hid_t file_id,set_id, space_id;
        hsize_t dims[2], maxdims[2];
        hid_t  mem_space_id;

        hsize_t start[2]= {(hsize_t) 0, (hsize_t) 0};
        hsize_t stride[2]={(hsize_t) 1, (hsize_t) 1};
        hsize_t count[2] ={(hsize_t) 1, (hsize_t) 1};


        if( (file_id=H5Fopen(filename,H5F_ACC_RDONLY, H5P_DEFAULT))==-1){
            printf("error in opening -%s-\n", filename);
            exit(FATAL_ERROR);
        }

        set_id=H5Dopen(file_id,"component", H5P_DEFAULT);
        space_id=H5Dget_space(set_id);
        H5Sget_simple_extent_dims(space_id, dims, maxdims);

        nwave_pca=dims[0];
        npc=dims[1];

        pca_table=(float32 *)malloc(nwave_pca*npc*sizeof(float32));
        pc=(double *)malloc(npc*sizeof(double));
        wave_pca=(int *)malloc(nwave_pca*sizeof(int));
        bindx=(int *)malloc(nwave_pca*sizeof(int));
        pcc_all=(float *)malloc(3*npix*sizeof(float));

        count[0]=nwave_pca;
        count[1]=npc;

        mem_space_id=H5Screate_simple(2, count, NULL);
        H5Sselect_hyperslab(space_id, H5S_SELECT_SET, start, stride, count, NULL);
        H5Dread(set_id,H5T_NATIVE_FLOAT, mem_space_id, space_id, H5P_DEFAULT, (void *)pca_table);

        H5Sclose(space_id);
        H5Dclose(set_id);

        set_id=H5Dopen(file_id,"wavelength", H5P_DEFAULT);
        space_id=H5Dget_space(set_id);
        mem_space_id=H5Screate_simple(1, count, NULL);
        H5Sselect_hyperslab(space_id, H5S_SELECT_SET, start, stride, count, NULL);
        H5Dread(set_id,H5T_NATIVE_INT, mem_space_id, space_id, H5P_DEFAULT, (void *)wave_pca);

        H5Sclose(space_id);
        H5Dclose(set_id);

        H5Fclose(file_id);

        for(i=0;i<nwave_pca;i++)
            bindx[i]=windex(1.0*wave_pca[i],wave,nbands);

        //read the coefficients from msl12_default.par
        filedir=strrchr(filename,'/');
        filename[filedir-filename]='\0';
        strcat(filename, "/picophyt.json");

        json_t *json,*wv;
        json_error_t error;
        const char *line;

        json = json_load_file(filename, 0, &error);
        if (!json){
            printf("the file %s not exist\n",filename);
            exit(1);
        }

        wv=json_object_get(json,"npc");
        line=json_string_value(wv);
        npc_used=atof(line);

        pro_coef=(double *)malloc((npc+2)*sizeof(double));
        syn_coef=(double *)malloc((npc+1)*sizeof(double));
        apeuk_coef=(double *)malloc((npc+1)*sizeof(double));
        interp_rrs=(float *)malloc(nwave_pca*sizeof(float));

        wv=json_object_get(json,"pro_coef");
        line=json_string_value(wv);
        i=0;
        while(line){
            pro_coef[i++]=atof(line);
            line=strchr(line,',');
            if(line)
                line++;
        }

        wv=json_object_get(json,"syn_coef");
        line=json_string_value(wv);
        i=0;
        while(line){
            syn_coef[i++]=atof(line);
            line=strchr(line,',');
            if(line)
                line++;
        }

        wv=json_object_get(json,"apeuk_coef");
        line=json_string_value(wv);
        i=0;
        while(line){
            apeuk_coef[i++]=atof(line);
            line=strchr(line,',');
            if(line)
                line++;
        }

    }

    if (input->proc_sst)
        sst = l2rec->sst;
    else
        sst=l2rec->l1rec->sstref;

    for(ip=0;ip<npix;ip++){

        ipb=ip*nbands;

        cell_abundance[ip]=BAD_FLT;
        for(i=0;i<3;i++)
            pcc_all[ip*3+i]=BAD_FLT;

         //interpolate input Rrs to the pca wavelength
        for(i=0;i<nwave_pca;i++){
            for(j=0;j<nbands;j++){
                if(wave_pca[i]< wave[j])
                    break;
            }
            if(j==0)
                j+=1;
            else if(j==nbands)
                j=nbands-1;

            //interpolate Rrs in pca wavelength using input Rrs at wavelengths j and j-1

            interp_rrs[i]=Rrs[ipb+j-1]+(Rrs[ipb+j]-Rrs[ipb+j-1])*(wave_pca[i]-wave[j-1])/(wave[j]-wave[j-1]);
        }

        rrs_standardize(interp_rrs,nwave_pca);

        for(i=0;i<npc;i++){
            pc[i]=0.;
            for(j=0;j<nwave_pca;j++){
                pc[i]+=pca_table[j*npc+i]*interp_rrs[j];
            }
        }

        pcc_all[ip*3]=pro_coef[0]+pro_coef[1]*log10(sst[ip]);
        pcc_all[ip*3+1]=syn_coef[0];
        pcc_all[ip*3+2]=apeuk_coef[0];
        for(i=0;i<npc_used;i++){
            pcc_all[ip*3]+=pro_coef[i+2]*pc[i];
            pcc_all[ip*3+1]+=syn_coef[i+1]*pc[i];
            pcc_all[ip*3+2]+=apeuk_coef[i+1]*pc[i];
        }
        pcc_all[ip*3+1]=pow(10.,pcc_all[ip*3+1]);
        pcc_all[ip*3+2]=pow(10.,pcc_all[ip*3+2]);
        if(pcc_all[ip*3+1]<0.)
            pcc_all[ip*3+1]=BAD_FLT;
        if(pcc_all[ip*3+2]<0.)
            pcc_all[ip*3+2]=BAD_FLT;

        if(pcc_all[ip*3]<0 && pcc_all[ip*3+1]>0. && pcc_all[ip*3+2]>0.)
            pcc_all[ip*3]=0.;

        switch(p->cat_ix){
        case CAT_prochlorococcus:
            cell_abundance[ip]=pcc_all[ip*3];
            break;
        case CAT_synechococcus:
            cell_abundance[ip]=pcc_all[ip*3+1];
            break;
        case CAT_autotrophic_picoeukaryotes:
            cell_abundance[ip]=pcc_all[ip*3+2];
            break;
        default:
            break;
        }

    }

}

