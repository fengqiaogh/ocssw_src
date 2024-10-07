#ifndef MAPPROTO_H
#define MAPPROTO_H

#include <stdint.h>
#include <map.h>

int32_t get_l3m(char *l3m_path, int16_t *syear, int16_t *sday, int32_t *smsec,
        int16_t *eyear, int16_t *eday, int32_t *emsec,
        char *prod_type, char *l3m_name, uint8_t *l3m_data,
        unsigned char *palette, meta_struct *meta_l3m);

int32_t read_meta(int32_t sdfid, meta_struct *meta_l3m);

int32_t getattrsz(int32_t id, char *attr_name, int32_t *nt, int32_t *count);

int32_t rdattr(int32_t sdfid, char *attr_name, void *buf);

#endif /* MAPPROTO_H */
