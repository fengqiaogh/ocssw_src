#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include "get_nitrate.h"

/*
 *  Calculating Global and Basin Scale Maps of Sea surface Nitrate from Chlorophyll and SST using Empirical Algorithms.
 *  Modified by Joaquim Goes - 24 Aug 2017
 *  Implemented into l2gen by Robert Lossing May 2019
 */

float calc_nitrate(float chl, float sst, float lat, float lon) {

// Initiating Variables

	float no3;
	float nitrate = BAD_FLT;

	if (chl < 0) {
		return nitrate;
	}

	if (chl == BAD_FLT) {
		return nitrate;
	}

	no3 = 25.68 - 1.97 * sst + 0.04 * pow(sst, 2) - 1.63 * chl
			+ 0.012 * pow(chl, 2);

// If Equatorial (Latitude: 30 degrees N and 18.33 degrees S):

	if (lat <= 30 && lat >= -18.333 && sst > 20 && chl > 0.24 && chl < 10) {
		no3 = 345.41 - 23.4 * sst + 0.4 * pow(sst, 2) + 0.012 * chl
				- pow(chl, 2) * 0.00048;
	}

// If Tropical Gyres:

	if (sst > 28 && chl >= 0.3 && chl <= 2) {
		no3 = 90.12 - 6.02 * sst + 0.09 * pow(sst, 2) - .05 * chl;
	}

// If Antartic (latitude) * and Arctic (latitude):

	if (sst < 2 && chl <= 2.2 && chl >= 0.1) {
		no3 = 12.98 - 13.61 * sst - 4.16 * pow(sst, 2) - 4.97 * chl
				+ 0.267 * pow(chl, 2);
	}

// Condition 4 (Arabian Sea 45 degrees E and 86.666 E )

	if (lon >= 45 && lon <= 86.666 && lat <= 30 && lat >= 5) {
		no3 = 119.14 - 8.38 * sst + 0.15 * pow(sst, 2) - .2 * chl;
	}

// All other regions:

	nitrate = no3;

	if (nitrate < 0) {
		nitrate = 0;
	}
	if (nitrate > 40) {
		nitrate = 40;
	}
	return nitrate;

}

void get_nitrate(l2str *l2rec, l2prodstr *p, float prod[]) {
	int i;
	int16_t year, day;
	double sec;
	unix2yds(l2rec->l1rec->scantime, &year, &day, &sec);

	l1str *l1rec = l2rec->l1rec;
	l2rec->sst = get_sst(l2rec);

	for (i = 0; i < l1rec->npix; i++) {
		if (l2rec->chl[i] == BAD_FLT) {
			prod[i] = BAD_FLT;
			continue;
		}

		switch (p->cat_ix) {

		case CAT_nitrate:
			prod[i] = calc_nitrate(l2rec->chl[i], l2rec->sst[i], l1rec->lat[i], l1rec->lon[i]);
			break;

		default:
			printf("get_nitrate.c can not produce product %s\n",
					p->algorithm_id);
			exit(1);
		}

		if (isnan(prod[i])) {
			prod[i] = BAD_FLT;
			l1rec->flags[i] |= PRODFAIL;

		}
	}
}
