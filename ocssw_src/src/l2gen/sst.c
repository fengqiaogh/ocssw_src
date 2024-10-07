#include "sst.h"


/**
 * @brief Retrive the dust SST correction for MODIS and product SST * 
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_sst_dust_correction(l2str *l2rec)
{
    float *dust = NULL;
    extern l1qstr l1que;
    call_dust_correction(l2rec->l1rec, &l1que, input, "SST", &dust); 
    return (dust);
}

/**
 * @brief Retrive the bias SSES for SST (both VIIRS and MODIS)
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_bias_sst(l2str *l2rec)
{
    float *sses_var = NULL;
    extern l1qstr l1que;
    call_sses_bias(l2rec->l1rec, &l1que, input, "SST", &sses_var); 
    return (sses_var);
}

/**
 * @brief Retrive the bias SSES for SST4 (MODIS)
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_bias_sst4(l2str *l2rec)
{
    float *sses_var = NULL;
    extern l1qstr l1que;
    call_sses_bias(l2rec->l1rec, &l1que, input, "SST4", &sses_var); 
    return (sses_var);
}

/**
 * @brief Retrive the bias SSES for SST3 (VIIRS)
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_bias_sst_triple(l2str *l2rec)
{
    float *sses_var = NULL;
    extern l1qstr l1que;
    call_sses_bias(l2rec->l1rec, &l1que, input, "SST3", &sses_var); 
    return (sses_var);
}

/**
 * @brief Retrive the STD bias SSES for SST
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_stdv_sst(l2str *l2rec)
{
    float *sses_var = NULL;
    extern l1qstr l1que;
    call_sses_std(l2rec->l1rec, &l1que, input, "SST", &sses_var);
    return (sses_var);
}

/**
 * @brief Retrive Counts SSES for SST
 * @param l2rec - l2 record
 * @return int16*  - pointer to array
 */
int16 *get_counts_sst(l2str *l2rec)
{
    int16 *counts = NULL;
    extern l1qstr l1que;
    call_sses_counts(l2rec->l1rec, &l1que, input, "SST", &counts);
    return counts;
}

/**
 * @brief Retrive Counts SSES for SST4
 * @param l2rec - l2 record
 * @return int16*  - pointer to array
 */
int16 *get_counts_sst4(l2str *l2rec)
{
    int16 *counts = NULL;
    extern l1qstr l1que;
    call_sses_counts(l2rec->l1rec, &l1que, input, "SST4", &counts);
    return counts;
}

/**
 * @brief Retrive Counts SSES for SST3
 * @param l2rec - l2 record
 * @return int16*  - pointer to array
 */
int16 *get_counts_sst_triple(l2str *l2rec)
{
    int16 *counts = NULL;
    extern l1qstr l1que;
    call_sses_counts(l2rec->l1rec, &l1que, input, "SST3", &counts);
    return counts;
}

/**
 * @brief Retrive the STD bias SSES for SST4
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_stdv_sst4(l2str *l2rec)
{
    float *sses_var = NULL;
    extern l1qstr l1que;
    call_sses_std(l2rec->l1rec, &l1que, input, "SST4", &sses_var); // we pass here the SST
    return (sses_var);
}

/**
 * @brief Retrive the STD bias SSES for SST3
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_stdv_sst_triple(l2str *l2rec)
{
    float *sses_var = NULL;
    extern l1qstr l1que;
    call_sses_std(l2rec->l1rec, &l1que, input, "SST3", &sses_var); // we pass here the SST
    return (sses_var);
}

/**
 * @brief Retrive Qual level (0-5) for SST
 * @param l2rec - l2 record
 * @return int8*  - pointer to array
 */
int8 *get_qual_sst(l2str *l2rec)
{
    int8 *qual_flags = NULL;
    extern l1qstr l1que;
    call_qual_flags(l2rec->l1rec, &l1que, input, "SST", &qual_flags); // we pass here the SST
    return qual_flags;

}

/**
 * @brief Retrive Qual level for SST4
 * @param l2rec - l2 record
 * @return int8*  - pointer to array
 */
int8 *get_qual_sst4(l2str *l2rec)
{
    int8 *qual_flags = NULL;
    extern l1qstr l1que;
    call_qual_flags(l2rec->l1rec, &l1que, input, "SST4", &qual_flags); // we pass here the SST
    return qual_flags;
}

/**
 * @brief Retrive Qual level for SST3
 * @param l2rec - l2 record
 * @return int8*  - pointer to array
 */
int8 *get_qual_sst_triple(l2str *l2rec)
{
    int8 *qual_flags = NULL;
    extern l1qstr l1que;
    call_qual_flags(l2rec->l1rec, &l1que, input, "SST3", &qual_flags); // we pass here the SST
    return qual_flags;
}

/**
 * @brief Retrive Cloud Mask Flags for SST
 * @param l2rec - l2 record
 * @return int16*  - pointer to array
 */
int16 *get_flags_sst(l2str *l2rec)
{
    int16 *flags_sst = NULL;
    extern l1qstr l1que;
    call_sst_flags(l2rec->l1rec, &l1que, input, "SST", &flags_sst);
    return flags_sst;

}

/**
 * @brief Retrive Cloud Mask Flags for SST4
 * @param l2rec - l2 record
 * @return int16*  - pointer to array
 */
int16 *get_flags_sst4(l2str *l2rec)
{
    int16 *flags_sst = NULL;
    extern l1qstr l1que;
    call_sst_flags(l2rec->l1rec, &l1que, input, "SST4", &flags_sst);
    return flags_sst;
}

/**
 * @brief Retrive Cloud Mask Flags for SST3
 * @param l2rec - l2 record
 * @return int16*  - pointer to array
 */
int16 *get_flags_sst_triple(l2str *l2rec)
{
    int16 *flags_sst = NULL;
    extern l1qstr l1que;
    call_sst_flags(l2rec->l1rec, &l1que, input, "SST3", &flags_sst);
    return flags_sst;
}

/**
 * @brief Get the sst for  SST
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_sst(l2str *l2rec)
{
    float *sst_cpp = NULL;
    extern l1qstr l1que;

    call_sst(l2rec->l1rec, &l1que, input, "SST", &sst_cpp);    
    return sst_cpp;
}

/**
 * @brief Get the sst for  SST4
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_sst4(l2str *l2rec)
{
    float *sst_cpp = NULL;
    extern l1qstr l1que;
    call_sst(l2rec->l1rec, &l1que, input, "SST4", &sst_cpp);   
    return sst_cpp;
}

/**
 * @brief Get the sst for  SST3
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_sst_triple(l2str *l2rec)
{
    float *sst_cpp = NULL;
    extern l1qstr l1que;
    call_sst(l2rec->l1rec, &l1que, input, "SST3", &sst_cpp);   
    return sst_cpp;
}

/**
 * @brief Get the means SSES bias for  SST
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_bias_mean_sst(l2str *l2rec)
{
    float *sses_var = NULL;
    extern l1qstr l1que;
    call_sses_bias_mean(l2rec->l1rec, &l1que, input, "SST", &sses_var); 
    return (sses_var);
}

/**
 * @brief Get the means SSES bias for  SST4
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_bias_mean_sst4(l2str *l2rec)
{
    float *sses_var = NULL;
    extern l1qstr l1que;
    call_sses_bias_mean(l2rec->l1rec, &l1que, input, "SST4", &sses_var); 
    return (sses_var);
}

/**
 * @brief Get the means SSES bias for  SST3
 * @param l2rec - l2 record
 * @return float*  - pointer to array
 */
float *get_bias_mean_sst_triple(l2str *l2rec)
{
    float *sses_var = NULL;
    extern l1qstr l1que;
    call_sses_bias_mean(l2rec->l1rec, &l1que, input, "SST3", &sses_var); 
    return (sses_var);
}


