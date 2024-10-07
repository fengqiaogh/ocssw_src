/*
 * Program to decode Channel Access Data Units following
 * CCSDS TM Sync and Channel Coding as proposed for use on Outernet.
 *
 * Mark McCrum 
 */

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     FutureTech     10/18/18 0.90  Continue processing after
//                                               decode error.
//  Liang Hong     SAIC           05/05/23 1.00  moved excess warning to verbose

#include "fec.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define I 			4		/* Interleave depth */
#define L_SYMBOL		223		/* RS Symbol size (bytes) */
#define L_CODEWORD              255		/* RS Codeword size (bytes) */
#define L_PARITY                32		/* Codeword parity bytes */
#define L_CODEBLOCK		(I * L_CODEWORD)
#define FRAME_SIZE 		(I * L_SYMBOL)  /* User data in a frame */
#define SYMBOLS_NO_CONV         (L_ASM + L_CODEBLOCK)
#define MAX_SYMBOLS		( 2 * (L_ASM + L_CODEBLOCK + 1))
						/* Max symbols in a CADU */
#define L_ASM			4		/* Length of sync vector */

/* Viterbi decoded config */
#define V_FRAMEBITS		(8 * SYMBOLS_NO_CONV) /* Bits in frame */

#define VERSION "1.00"

/**
 * Formatted hex output
 */
void formatHex(unsigned char* data, int length) {
  int col = 0;
  while(length) {
    printf("0x%02hhx", *(data++));
    length--;
    if (length) {
      printf(", ");
    }
    col += 4;
    if (col > 80) {
      col = 0;
      putchar('\n');
    }
  }
  putchar('\n');
}

/**
 * Unformatted hex output
 */
void outputHex(unsigned char* data, int length) {
  while(length) {
    printf("%02hhx", *(data++));
    length--;
  }
}

/*
 * LFSR for CCSDS Pseudo-random number generator
 *
 * lfsr Pointer to state of linear feedback shift register. Should initially
 *      be 0xFF
 *
 * o    Pointer to output byte, existing value will be xor'd with the next
 *      random number output by the generator.
 */
void prng(unsigned char* lfsr, unsigned char* o) {
  unsigned char rnd = 0U;
  unsigned char bit = 0U;
  unsigned int i;
  for (i = 0; i < 8; i++) {
    rnd = (rnd << 1) | (*lfsr & 0x1U);
    bit = (((*lfsr >> 3) ^ *lfsr) ^ (*lfsr >> 5)) ^ (*lfsr >>7);
    *lfsr = (*lfsr >> 1) | (bit << 7);
  }	
  *o = *o ^ rnd;
} 

void convDecode(unsigned char* data, unsigned char* coded) {
  void *vp = create_viterbi27(V_FRAMEBITS);
  int polys[] = {V27POLYB, -V27POLYA};	/* CCSDS Polynomial config */
  unsigned char symbols[2 * (V_FRAMEBITS + 6)]; /* Channel symbols inc.
                                                        tail bits */
  int i;

  /* build array of symbols (one element per bit) includes 6 tail bits*/
  for (i = 0; i < 2 * (V_FRAMEBITS + 6); i++) {
    symbols[i ] = ((coded[i/8] << (i % 8)) & 0x80) > 0 ? 255 : 0;
  }

  set_viterbi27_polynomial(polys);
  init_viterbi27(vp, 0);

  /* Decode block */
  update_viterbi27_blk(vp, symbols, V_FRAMEBITS + 6);
      
  /* Do Viterbi chainback */
  chainback_viterbi27(vp, data, V_FRAMEBITS,0);

  delete_viterbi27(vp);
}

void showHelp(char *progname) {
  printf("\nTest decoded for Outernet\n");
  printf("\nDecode Channel Access Data Units encoded via propsed use of");
  printf("\nCCSDS TM Sync and Channel Coding standard\n\n");
  printf("Use: %s [options]\n", progname);
  printf(" c Indicates that CADU has been convolutionally coded\n"); 
  printf("   A Viterbi decoder will be used\n");
  printf("Output:\n");
  printf(" Decoded frame\n");
  printf("Notes:\n");
  printf(" Errors from the Reed-Solomon decode are detected but not\n");
  printf(" corrected.\n"); 
  exit(0);
}


int main(int argc, char** argv) {
  unsigned char sync[] = {0x1A, 0xCF, 0xFC, 0x1D};  /* Attached sync marker */
  unsigned char symbols[MAX_SYMBOLS];		    /* CADU symbols         */
  unsigned char convdecoded[SYMBOLS_NO_CONV + 1];
						    /* frame after conv decoding */
  //unsigned char frame[ FRAME_SIZE ];		    /* Transfer frame       */
  unsigned char codeword[L_CODEWORD];               /* Single RS code word*/
  unsigned char lfsr = 0xFFU;	    /* Shift register for randomiser  */	
  //unsigned char csr = 0x00;         /* Shift register for convolutional coder */
  //unsigned char flush = 0x00;       /* Flush byte for convolutional coder */
  int rsErrCnt;			    /* RS error count */
  //int pad = 0;                      /* RS coding padding */
  int i, j;                         /* Loop counters */
  int d;                            /* Getopt option */
  int convolution = 0;              /* Do convolutional coding? */
  int verbose = 0;

  FILE* in;// = stdin;		    /* Input file */
  FILE* out;		            /* Output file */
  int caduSymbols;
  int errPos[L_PARITY];             /* Error positions */

  printf("outernetDecode %s (%s %s)\n", VERSION, __DATE__, __TIME__);

  if ( argc == 1) {
    printf("\nouternetDecode input_downlink_filename output_telemetry_filename\n");
    return 0;
  }

  /* Process command line options */
  while ((d = getopt(argc, argv, "cv")) != EOF) {
    switch(d) {
    case 'c':
      convolution = 1;
      break;
    case 'v':
      verbose = 1;
      break;
    case '?':
      showHelp(argv[0]);
    }
  }

  printf("MAX_SYMBOLS: %d\n", MAX_SYMBOLS);
  printf("FRAME_SIZE: %d\n", FRAME_SIZE);
  
  /* Number of bytes to read per frame depends on convolutional coding */
  caduSymbols = convolution ? MAX_SYMBOLS : SYMBOLS_NO_CONV;

  in = fopen(argv[optind], "r");
  out = fopen(argv[optind+1], "w");

  int ierr=0;
  while (fread(symbols, sizeof(char), caduSymbols, in) == caduSymbols) {

    lfsr = 0xFFU;	/* reset the randomiser shift register */
    memset(convdecoded, 0, SYMBOLS_NO_CONV);
    if (convolution) {
      /* Run a viterbi decoder on the data */
      convDecode(convdecoded, symbols);
    } else {
      memcpy(convdecoded, symbols, SYMBOLS_NO_CONV);
    }
    /* Check the attached sync marker */
    if(memcmp(sync, convdecoded, L_ASM) != 0) {
      fprintf(stderr, "ASM not found\n");
      //      break;
    }

    /* De-randomise the data */
    for (i = 0; i < ( L_CODEBLOCK); i++) {
      prng(&lfsr, &convdecoded[L_ASM + i]);
    }

    /* Check parity */
    for (i = 0; i < I; i++) {
      /* codeword interleaving */
      for (j = 0; j < L_SYMBOL; j++) {
        codeword[j] = convdecoded[L_ASM + i + j * I];
      }
		
      /* parity interleaving */
      for(j = 0; j < L_PARITY; j++) {
        codeword[j + L_SYMBOL] = 
          convdecoded[L_ASM + FRAME_SIZE + i + j * I];
      }
		
      /* run the RS decoder */
      rsErrCnt = decode_rs_ccsds(codeword, errPos, 0, 0);
      if (rsErrCnt == -1) {
        if (verbose > 0) fprintf(stderr, "Uncorrectable errors in frame\n");
        ierr = 1;
        //        break;
      } else if (rsErrCnt > 0 & verbose > 0) {
        /* TODO: Copy the corrections back into the frame */
        fprintf(stderr, "%d errors in block\n", rsErrCnt);
        while(rsErrCnt) {
          fprintf(stderr, "symbol: %i\n", 
                  errPos[--rsErrCnt]);			
        }
      }
    }
    if (ierr == 0) {
      /* Output the frame */
      //      outputHex(&convdecoded[L_ASM], FRAME_SIZE);
      //putchar('\n');
      fwrite(&convdecoded[L_ASM], sizeof(char), FRAME_SIZE, out);
    } else {
      ierr = 0;
    }
  }

  fclose(in);
  fclose(out);
  
  return 0;
}

