#ifndef DATATYPE_PROTO_H_
#define DATATYPE_PROTO_H_

extern char *get_data_type_str_from_name
PROTO((char *fname));

extern int get_data_type_from_name
PROTO((char *fname));

extern int get_data_type_from_ffm
PROTO((byte rec[21504]));

extern char *get_data_type_str
PROTO((int datatype));


#endif /* DATATYPE_PROTO_H_ */
