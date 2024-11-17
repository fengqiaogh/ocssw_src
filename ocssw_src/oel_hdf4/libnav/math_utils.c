#include "math_utils.h"

void transpose3d(const float inp[3][3], float out[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i][j] = inp[j][i];
        }
    }
}
