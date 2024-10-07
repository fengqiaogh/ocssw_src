#include "l12_proto.h"

int32 prodlist(int32 sensorID, int32 evalmask, const char *inprod, const char *defprod, char outprod[L1_MAXPROD][32]) {
    static int32 nwave = -1;
    static int32 nwaveIR = -1;
    static int32 nwaveVIS = -1;
    static int32 *wave = NULL;

    int32 tot_prod;
    char wavestr[6];
    char prefix[256];
    char suffix[256];
    char *p1, *p2;
    int32 i, j, k, iw;
    char *tmp_ptr;
    char tmp_str[2048];

    if (wave == NULL) {
        nwave = rdsensorinfo(sensorID, evalmask, "Lambda", (void **) &wave);
        nwaveIR = rdsensorinfo(sensorID, evalmask, "NbandsIR", NULL);
        nwaveVIS = rdsensorinfo(sensorID, evalmask, "NbandsVIS", NULL);
    }

    // insert default products if requested or no prods specified

    if (inprod[0] == '\0') {
        strcpy(tmp_str, defprod);
    } else {
        strcpy(tmp_str, inprod);
        tmp_ptr = tmp_str;
        while (*tmp_ptr != '\0') {
            if (isupper(*tmp_ptr)) *tmp_ptr = tolower(*tmp_ptr);
            tmp_ptr++;
        }
        if ((tmp_ptr = strstr(tmp_str, "default")) != NULL) {
            i = strlen(tmp_ptr);
            k = strlen(inprod);
            strncpy(tmp_str, inprod, k - i);
            tmp_str[k - i] = '\0';
            strcat(tmp_str, defprod);
            strcat(tmp_str, inprod + k - i + 7);
        } else
            strcpy(tmp_str, inprod);
    }

    // count products and separate into array of strings

    tot_prod = 0;
    if (tmp_str[0] != '\0') {
        i = 0;
        while (1) {
            while (tmp_str[i] == ' ' || tmp_str[i] == ',')
                i++;
            if (tmp_str[i] == '\0')
                break;
            k = i;
            while (tmp_str[i] != ' ' &&
                    tmp_str[i] != ',' &&
                    tmp_str[i] != '\0')
                i++;

            // exclude l2_flags, it will be explicitly added if needed
            if (strstr(&tmp_str[k], "l2_flags") == &tmp_str[k])
                continue;
            strncpy(outprod[tot_prod], &tmp_str[k], i - k);
            outprod[tot_prod][i - k] = '\0';
            tot_prod++;
            if (tmp_str[i] == '\0')
                break;
        }
    }
    if (tot_prod == 0) {
        return (tot_prod);
    }

    // expand any nnn wavelength placeholders with sensor wavelengths
    for (i = 0; i < tot_prod; i++) {
        if ((tmp_ptr = strstr(outprod[i], "_nnn")) != NULL) {

            // make room
            for (k = tot_prod - 1; k > i; k--)
                strcpy(outprod[k + nwave - 1], outprod[k]);
            tot_prod = tot_prod + nwave - 1;

            // replace
            p1 = outprod[i];
            p2 = tmp_ptr;
            prefix[0] = '\0';
            strncpy(prefix, p1, p2 - p1 + 1);
            prefix[p2 - p1 + 1] = '\0';
            suffix[0] = '\0';
            strcpy(suffix, p2 + 4);
            outprod[i][0] = '\0';

            for (iw = 0; iw < nwave; iw++) {
                wavestr[0] = '\0';
                sprintf(wavestr, "%d", (int) wave[iw]);
                outprod[i + iw][0] = '\0';
                strcat(outprod[i + iw], prefix);
                strcat(outprod[i + iw], wavestr);
                strcat(outprod[i + iw], suffix);
            }
        }
    }

    // expand any vvv wavelength placeholders with visible sensor wavelengths

    for (i = 0; i < tot_prod; i++) {
        if ((tmp_ptr = strstr(outprod[i], "_vvv")) != NULL) {

            // make room
            for (k = tot_prod - 1; k > i; k--)
                strcpy(outprod[k + nwaveVIS - 1], outprod[k]);
            tot_prod = tot_prod + nwaveVIS - 1;

            // replace
            p1 = outprod[i];
            p2 = tmp_ptr;
            prefix[0] = '\0';
            strncpy(prefix, p1, p2 - p1 + 1);
            prefix[p2 - p1 + 1] = '\0';
            suffix[0] = '\0';
            strcpy(suffix, p2 + 4);
            outprod[i][0] = '\0';

            for (iw = 0; iw < nwaveVIS; iw++) {
                wavestr[0] = '\0';
                sprintf(wavestr, "%d", (int) wave[iw]);
                outprod[i + iw][0] = '\0';
                strcat(outprod[i + iw], prefix);
                strcat(outprod[i + iw], wavestr);
                strcat(outprod[i + iw], suffix);
            }
        }
    }

    // expand any iii wavelength placeholders with IR sensor wavelengths
    for (i = 0; i < tot_prod; i++) {
        if ((tmp_ptr = strstr(outprod[i], "_iii")) != NULL) {

            // make room
            for (k = tot_prod - 1; k > i; k--)
                strcpy(outprod[k + nwaveIR - 1], outprod[k]);
            tot_prod = tot_prod + nwaveIR - 1;

            // replace
            p1 = outprod[i];
            p2 = tmp_ptr;
            prefix[0] = '\0';
            strncpy(prefix, p1, p2 - p1 + 1);
            prefix[p2 - p1 + 1] = '\0';
            suffix[0] = '\0';
            strcpy(suffix, p2 + 4);
            outprod[i][0] = '\0';

            for (iw = 0; iw < nwaveIR; iw++) {
                wavestr[0] = '\0';
                sprintf(wavestr, "%d", (int) wave[iw + nwave]);
                outprod[i + iw][0] = '\0';
                strcat(outprod[i + iw], prefix);
                strcat(outprod[i + iw], wavestr);
                strcat(outprod[i + iw], suffix);
                //printf("%d %d %d %s %d\n",iw,iw+nwave,wave[iw+nwave],outprod[i+iw],tot_prod);
            }
        }
    }

    // get rid of duplicate products from output file
    for (i = 0; i < tot_prod; i++) {
        for (j = i + 1; j < tot_prod; j++) {
            if (strcmp(outprod[i], outprod[j]) == 0) {
                for (k = j + 1; k < tot_prod; k++) {
                    strcpy(outprod[k - 1], outprod[k]);
                }
                tot_prod--;
                j--;
            }
        }
    }

    return (tot_prod);
}
