#ifndef REGEN_PROTO_H
#define REGEN_PROTO_H

extern int32 regen(char *infile, int32 *spix, int32 *epix, int32 *sscan,
        int32 *escan, int32 pix_sub, int32 sc_sub, char *parm_list,
        char *outfile);

extern void set_subsc_params(int32 nsamp, int32 nrec, int32 xsub, int32 ysub,
        int32 *ssamp, int32 *esamp, int32 *srec, int32 *erec,
        int32 *subsc_samps, int32 *subsc_recs, int32 *max_samp_used,
        int32 *max_rec_used);

extern int32 rdattr(int32 sdfid, char *attr_name, void *buf);

extern int32 getset_gattrs(int32 isdfid, int32 osdfid);

extern int32 dupHDF(int32 ifid, int32 ofid, int32 isdfid, int32 osdfid,
        int32 *ssamp, int32 *esamp, int32 *srec, int32 *erec, int32 xsub,
        int32 ysub, int32 subsc_samps, int32 subsc_recs, int32 nvgps,
        char *parm_list);

extern int32 set_sds(int32 isdfid, int32 osdfid, char *vgname, int32 ovid,
        int32 tag, int32 ref, int32 subsc_samps, int32 subsc_recs,
        char *parm_list, int32 *in_sdsid, int32 *out_sdsid);

extern int32 write_data(int32 isdfid, int32 osdfid, int32 isdsid, int32 osdsid,
        int32 ssamp, int32 srec, int32 esamp, int32 erec, int32 xsub,
        int32 ysub);

extern int32 rdslice(int32 sdfid, char *name, int32 *start, int32 *edge,
        void *buf);

extern int32 subsample(int32 rank, int32 *idims, int32 *odims,
        int32 num_type, int32 srec, int32 erec, int32 ssamp,
        int32 esamp, int32 xsub, int32 ysub, void *ibuf, void *obuf);

extern int32 update_sds(int32 sdfid, char *name, int32 *start, int32 *edge,
        void *buf);

extern int32 subsamp_rec(int32 srec, int32 maxrec, int32 rank, int32 *dims,
        int32 ysub, int32 nt, void *obuf, void *nbuf);

extern void subsamp_2D(int32 spix, int32 sscan, int32 epix, int32 escan,
        int32 xdim, int32 xsub, int32 ysub, int32 nt, void *inbuf,
        void *outbuf);

extern int32 alloc_nav_buffs(int32 nrec, NavType *nav_rec);

extern int32 get_navdata(int32 sdfid, int32 nrec, NavType *nav_rec);

extern int32 alloc_geonav_buffs(int32 nrec, GeoType *geo_rec);

extern int32 get_geodata(int32 pix_start, int32 pix_sub, int32 srec,
        int32 ssamp, int32 subsc_recs,
        int32 subsc_samps, int32 ysub, int32 xsub, char *dtype,
        NavType *nav_rec, GeoType *geo_rec);

extern int32 write_coords(int32 sdfid, int32 nrec, GeoType *geo_rec);

extern int32 get_tiltdata(int32 sdfid, tilt_Type *tiltrec);

extern int32 set_tiltdata(int32 sdfid, int32 srec, int32 erec, int32 ysub,
        float32 *slat, float32 *slon, float32 *elat, float32 *elon,
        tilt_Type *old_tiltrec, tilt_Type *new_tiltrec);

extern int32 set_l1adata(int32 isdfid, int32 osdfid, int32 num_recs,
        int32 num_samps, int32 subsc_recs, int32 subsc_samps,
        int32 srec, int32 ssamp, int32 erec, int32 esamp, int32 xsub,
        int32 ysub, FilemetricsType *fm_rec);

extern int32 subsamp_l1adata(int32 srec, int32 erec, int32 ssamp,
        int32 esamp, int32 xsub, int32 ysub, int32 nsamp,
        int16 *gain, int16 *dark_rest, FilemetricsType *fm_rec,
        int16 *inbuf, int16 *outbuf, int16 *s_satp, int16 *s_zerop,
        int32 *new_buf_recs);

extern int32 set_flags(int32 sdfid, int32 nrec, int32 nsamp, int32 *flags);

extern int32 set_globalattrs(char *outfile, int32 sdfid, int32 l1aflag,
        int32 nrec, int32 nsamp, int32 ssamp, int32 xsub, int32 *flags,
        GeoType *geo_rec, FilemetricsType *fm_rec);

extern void free_geonav_buffs(GeoType *geo_rec);

extern void free_nav_buffs(NavType *nav_rec);

#endif /* REGEN_PROTO_H */

