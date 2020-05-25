#ifndef FFT_DEFINE_H_
#define FFT_DEFINE_H_

#define DATA_TW 32         //Total binary bits used to represent data
#define DATA_FW 10          //Total binary bits used to represent fractional part of data
#define TF_TW 18           //Total binary bits used to represent twiddle factor
#define TF_FW 16           //Total binary bits used to represent fractional part of twiddle factor
#define MAX_FFT_LEN 8192   //Maximum FFT length
//#define DUMP_TO_FILE       //Dump data into text file
#define ROUNDING false     //Using truncate instead of rounding
//#define MICROSOFT_C

#endif
