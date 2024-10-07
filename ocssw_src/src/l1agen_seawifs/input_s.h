#ifndef INPUT_S_H_
#define INPUT_S_H_

#ifndef byte
#define byte unsigned char
#endif

typedef struct input_struct {
    int flag;
    short int sc_id[2];
    short int iyear;
    short int iday;
    int msec;
    float sc_ana[40];
    byte sc_dis[40];
    float inst_ana[5][40];
    byte inst_dis[5][32];
    int nflag[8];
} input_sType;

#endif /* INPUT_S_H_ */
