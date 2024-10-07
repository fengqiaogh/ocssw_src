/**
 *  @brief determine endianess
 *  
 *  Determines if the host calling this routine is a big or little
 *  endian machine.
 *  
 *  @returns A 0 for BIG_ENDIAN machines, 1 for LITTLE_ENDIAN mahcines.
 */

int endianess(void) {
    int x = 1;
    return ( *(char *) &x == 1);
}
