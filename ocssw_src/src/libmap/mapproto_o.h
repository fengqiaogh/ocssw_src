#ifndef MAPPROTO_O_H
#define MAPPROTO_O_H

#include <meta_l3b.h>

int32 put_l3m(char *l3m_path, char *replaces, int16 bin_syear, int16 bin_sday,
        int16 bin_eyear, int16 bin_eday, int16 syear, int16 sday,
        int32 smsec, int16 eyear, int16 eday, int32 emsec,
        float32 *lat_range, float32 *lon_range, int32 lines,
        int32 columns, char *flag_names, char *flag_use,
        uint8 *eng_q_use, char *ptime, char *infiles, char *prod_type,
        int32 nbins, char *l3m_name, void *l3m_data, char *measure,
        char *proc_con, char *proc_log, meta_l3bType *meta_l3b);

#endif /* MAPPROTO_O_H */
