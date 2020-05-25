//#include "fft_c.h"
#include "fft.h"
#include <stdint.h>
#include "stdio.h"

extern "C" void realFFT_forward( int32_t* poutData, int32_t* pinData, int nLen )
{
	
	fft fft_c(nLen, DATA_TW, DATA_FW, TF_TW, TF_FW, false);
	for(int i=0; i< nLen; i++){
		fft_c.fnSetInput(i, (float)(pinData[i])/1024.0, 0.0 );
	}

	fft_c.fnGenerateTw();
	fft_c.fnFftCalculate();
	fft_c.fnReorder();

	for ( int i = 0; i < nLen/2+1; i++ )
	{
		double fReal = fft_c.fnGetOutput(i, 0 );
		double fImg = fft_c.fnGetOutput(i, 1 );

		poutData[i*2] = (int32_t)(fReal * 1024 );
		poutData[i*2+1] = (int32_t)( fImg * 1024 );
	}

	return;
}

extern "C" void realFFT_Inverse( int32_t* poutData, int32_t* pinData, int nLen )
{	

	fft fft_c(nLen,DATA_TW, DATA_FW, TF_TW, TF_FW,true);
	fft_c.fnSetInput(0, (float)(pinData[0])/1024.0, 0 );
	for(int i=1; i< nLen/2; i++){
		fft_c.fnSetInput(i, (float)(pinData[i*2])/1024.0, (float)(pinData[i*2+1])/1024.0 );
		fft_c.fnSetInput(nLen-i, (float)(pinData[i*2])/1024.0, -(float)(pinData[i*2+1])/1024.0 );
	}

	fft_c.fnGenerateTw();
	fft_c.fnFftCalculate();
	fft_c.fnReorder();

	for ( int i = 0; i < nLen; i++ )
	{
		double fReal = fft_c.fnGetOutput(i, 0 );
		//double fImg = fft_c.fnGetOutput(i, 1 );

		poutData[i] = (int32_t)(fReal * 1024 );
		//pImgData[i] = (int32_t)( fImg * 1024 );
	}

}