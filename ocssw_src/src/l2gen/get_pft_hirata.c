#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

float calc_microplankton_hirata(float chl) {
    float a0 = 0.9117;
    float a1 = -2.7330;
    float a2 = 0.4003;
    float x = log10(chl);

    float ans = 1.0 / (a0 + exp(a1 * x + a2));
    if (ans > 1)
        ans = 1;
    if (ans < 0)
        ans = 0;
    return ans;
}

float calc_diatoms_hirata(float chl) {
    float a0 = 1.3272;
    float a1 = -3.9828;
    float a2 = 0.1953;
    float x = log10(chl);

    float ans = 1.0 / (a0 + exp(a1 * x + a2));
    if (ans > 1)
        ans = 1;
    if (ans < 0)
        ans = 0;

    return ans;
}

float calc_greenalgae_hirata(float chl) {
    float a0 = 0.2490;
    float a1 = -1.2621;
    float a2 = -0.5523;
    float x = log10(chl);
    float y = chl;

    float ans = (a0 / y) * exp(a1 * pow(x + a2, 2));
    if (ans > 1)
        ans = 1;
    if (ans < 0)
        ans = 0;
    return ans;
}

float calc_picoplankton_hirata(float chl) {
    float a0 = 0.1529;
    float a1 = 1.0306;
    float a2 = -1.5576;
    float a3 = -1.859;
    float a4 = 2.9954;
    float x = log10(chl);

    float ans = -(1.0 / (a0 + exp(a1 * x + a2))) + a3 * x + a4;
    if (ans > 1)
        ans = 1;
    if (ans < 0)
        ans = 0;
    return ans;
}

float calc_prokaryotes_hirata(float chl) {
    float a0 = 0.0067;
    float a1 = 0.6154;
    float a2 = -19.5190;
    float a3 = 0.9643;
    float a4 = 0.1027;
    float a5 = -0.1189;
    float a6 = 0.0626;
    float x = log10(chl);
    float y = chl;

    float ans = (a0 / a1 / y) * exp(((a2 * pow(x + a3, 2)) / pow(a1, 2))) + a4 * pow(x, 2) + a5 * x + a6;
    if (ans > 1)
        ans = 1;
    if (ans < 0)
        ans = 0;
    return ans;
}

float calc_prochlorococcus_hirata(float chl) {
    float a0 = .0099;
    float a1 = 0.6808;
    //float a2 = -8.6276;
    float a3 = 0.9668;
    float a4 = 0.0074;
    float a5 = -0.1621;
    float a6 = 0.0436;
    float x = log10(chl);
    float y = chl;

    float ans = (a0 / a1 / y) * exp(((a3 * pow(x + a4, 2))) / pow(a1, 2)) + a4 * pow(x, 2) + a5 * x + a6;
    if (ans > 1)
        ans = 1;
    if (ans < 0)
        ans = 0;

    return ans;
}

float calc_dinoflagellates_hirata(float chl) {
    float ans = calc_microplankton_hirata(chl) - calc_diatoms_hirata(chl);
    if (ans > 1)
        ans = 1;
    if (ans < 0)
        ans = 0;

    return ans;

}

float calc_nanoplankton_hirata(float chl) {
    float ans = 1.0 - calc_microplankton_hirata(chl) - calc_picoplankton_hirata(chl);
    if (ans > 1)
        ans = 1;
    if (ans < 0)
        ans = 0;

    return ans;
}

float calc_picoeukaryotes_hirata(float chl) {
    float ans = calc_picoplankton_hirata(chl) - calc_prokaryotes_hirata(chl);
    if (ans > 1)
        ans = 1;
    if (ans < 0)
        ans = 0;

    return ans;
}

float calc_prymnesiophytes_hirata(float chl) {
    float ans = calc_nanoplankton_hirata(chl) - calc_greenalgae_hirata(chl);
    if (ans > 1)
        ans = 1;
    if (ans < 0)
        ans = 0;

    return ans;
}

void get_pft_hirata(l2str *l2rec, l2prodstr *p, float prod[]) {
    int i;
    l1str *l1rec = l2rec->l1rec;

    for (i = 0; i < l1rec->npix; i++) {
        if (l2rec->chl[i] == BAD_FLT) {
            prod[i] = BAD_FLT;
            continue;
        }

        switch (p->cat_ix) {
        case CAT_microplankton_hirata:
            prod[i] = calc_microplankton_hirata(l2rec->chl[i]);
            break;
        case CAT_diatoms_hirata:
            prod[i] = calc_diatoms_hirata(l2rec->chl[i]);
            break;
        case CAT_greenalgae_hirata:
            prod[i] = calc_greenalgae_hirata(l2rec->chl[i]);
            break;
        case CAT_picoplankton_hirata:
            prod[i] = calc_picoplankton_hirata(l2rec->chl[i]);
            break;
        case CAT_prokaryotes_hirata:
            prod[i] = calc_prokaryotes_hirata(l2rec->chl[i]);
            break;
        case CAT_prochlorococcus_hirata:
            prod[i] = calc_prochlorococcus_hirata(l2rec->chl[i]);
            break;
        case CAT_dinoflagellates_hirata:
            prod[i] = calc_dinoflagellates_hirata(l2rec->chl[i]);
            break;
        case CAT_nanoplankton_hirata:
            prod[i] = calc_nanoplankton_hirata(l2rec->chl[i]);
            break;
        case CAT_picoeukaryotes_hirata:
            prod[i] = calc_picoeukaryotes_hirata(l2rec->chl[i]);
            break;
        case CAT_prymnesiophytes_hirata:
            prod[i] = calc_prymnesiophytes_hirata(l2rec->chl[i]);
            break;

        default:
            printf("get_pft_hirata can not produce product %s\n", p->algorithm_id);
            exit(1);
        }
        if (isnan(prod[i])) {
            prod[i] = BAD_FLT;
            l1rec->flags[i] |= PRODFAIL;
        }
    }


}

/*
x = log(Chl-a)
y = Chl-a
Microplankton = (0.9117+exp(-2.7330x+0.4003))^-1
Diatoms = [1.3272 + exp (-3.9828*x + 0.1953)]^−1
Dinoflagellates = Microplankton-Diatoms
Nanplankton = 1-Microplanton-Picoplankton
Greenalgae = (0.2490/y) exp[-1.2621(x +-.05523)^2]
Prymnesiophytes = Nanoplankton-Green Algae #(Haptophytes)
Picoplankton = – [0.1529+ exp (1.0306 *x + -1.5576)]^−1 + -1.8597*x + 2.9954
Prokaryotes = (.0067 /.6154 /y) exp [-19.5190 (x + 0.9643)^2 /0.0067^2 ] + 0.1027*x^2 + -0.1189*x + 0.0626
Picoeukaryotes = Picoplankton-Prokaryotes
Prochlorococcus = (0.0099 /0.6808/y) exp [0.9668 (x + 0.0074 )^2 /0.0099^2 ]+ 0.0074*x^2 + -0.1621*x + 0.0436

 */
