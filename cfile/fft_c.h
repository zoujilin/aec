#ifndef FFT_C_H
#define FFT_C_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
/*
 * fft foward, input is Q10 int32_t, inplace output
 */
//void FFT_forward( int32_t* pRealData, int32_t* pImgData, int nLen );
void realFFT_forward( int32_t* poutData, int32_t* pinData, int nLen );

/*
 * inverse fft 
 */
//void FFT_Inverse( int32_t* pRealData, int32_t* pImgData, int nLen );
void realFFT_Inverse( int32_t* poutData, int32_t* pinData, int nLen );
}
#endif

#endif
