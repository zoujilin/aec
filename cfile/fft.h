#ifndef FFT_H_
#define FFT_H_

#include "dtype.h"
#include "fft_define.h"


class fft
{
private:
	int iFftLength;
	int iTotalNumOfStage;
	int iTotalNumOfNode;
	int iTotalWidth;
	int iFracWidth;
	int iTwTotalWidth;
	int iTwFracWidth;
	bool bInversedFFt;

	complex ** caStage;
	complex * caTw;

public:
	fft();
	fft(int len = 512, int dw = 16, int fw = 8, int twdw = 16, int twfw = 8, bool ifft = false);
	~fft();

	void fnSetInput(int idx, double i, double j);
	double fnGetOutput(int idx, int flag);
	void fnDisplaySettings();
	void fnGenerateTw();
	void fnFftCalculate();
	void fnReorder();
	void fnReorder(int *out);
};



#endif /* FFT_H_ */

