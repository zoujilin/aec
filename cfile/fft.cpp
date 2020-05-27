#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "fft.h"
#include "fft_define.h"

using namespace::std;

fft::fft()
{
	iFftLength = 512;
	iTotalWidth = 16;
	iFracWidth = 8;
	iTwTotalWidth = 16;
	iTwFracWidth = 8;
	bInversedFFt = false;

	iTotalNumOfNode = iFftLength / 2;
	iTotalNumOfStage = int(log(iFftLength) / log(2));

	caStage = new complex *[int(log(iFftLength) / log(2)) + 2];

	for (int i = 0; i<int(log(iFftLength) / log(2)) + 2; i++)
	{
		caStage[i] = new complex[iFftLength];
	}

	for (int i = 0; i<int(log(iFftLength) / log(2)) + 2; i++)
	{
		for (int j = 0; j < iFftLength; j++)
		{
			caStage[i][j].set(0.0, 0.0, iTotalWidth, iFracWidth, ROUNDING);
		}
	}

	caTw = new complex[MAX_FFT_LEN];
}

fft::fft(int len, int dw, int fw, int twdw, int twfw, bool ifft)
{
	iFftLength = len;
	iTotalWidth = dw;
	iFracWidth = fw;
	iTwTotalWidth = twdw;
	iTwFracWidth = twfw;
	bInversedFFt = ifft;

	iTotalNumOfNode = iFftLength / 2;
	iTotalNumOfStage = int(log(iFftLength) / log(2));

	caStage = new complex *[int(log(iFftLength) / log(2)) + 2];

	for (int i = 0; i<int(log(iFftLength) / log(2)) + 2; i++)
	{
		caStage[i] = new complex[iFftLength];
	}

	for (int i = 0; i<int(log(iFftLength) / log(2)) + 2; i++)
	{
		for (int j = 0; j < iFftLength; j++)
		{
			caStage[i][j].set(0.0, 0.0, iTotalWidth, iFracWidth, ROUNDING);
		}
	}

	caTw = new complex[MAX_FFT_LEN];
}

fft::~fft()
{
	//cout << "Cleaning..." << endl;

	delete[] caStage;
	delete[] caTw;
}

void fft::fnDisplaySettings()
{
	cout << "-------------------------------------------" << endl;
	cout << "FFT length:                      " << iFftLength << endl;
	cout << "Data total width:                " << iTotalWidth << endl;
	cout << "Data fractional width:           " << iFracWidth << endl;
	cout << "Twiddle factor total width:      " << iTwTotalWidth << endl;
	cout << "Twiddle factor fractional width: " << iTwFracWidth << endl;
	cout << "-------------------------------------------" << endl;
}

void fft::fnSetInput(int idx, double i, double j)
{
	caStage[0][idx].set(i, j, DATA_TW, DATA_FW, ROUNDING);
}

double fft::fnGetOutput(int idx, int flag)
{
	double fValue;

	if (flag == 0)
	{
		fValue = caStage[iTotalNumOfStage + 1][idx].getReal();
	}
	else
	{
		fValue = caStage[iTotalNumOfStage + 1][idx].getImag();
	}

	return fValue;
}
void fft::fnGenerateTw()
{
	double fSinValue;
	double fCosValue;

	//cout.setf(ios::left);
	//cout.fill(0);
	//cout.precision(16);

	//cout << "Generating twiddle factor ..." << endl;
	//cout << "Defining PI: " << M_PI << endl;

	for (int i = 0; i < MAX_FFT_LEN; i++)
	{
		fSinValue = -sin(2 * M_PI*i / MAX_FFT_LEN);
		fCosValue = cos(2 * M_PI*i / MAX_FFT_LEN);

		caTw[i].set(fCosValue, fSinValue, TF_TW, TF_FW, ROUNDING);
	}

	//std::cout.unsetf(std::ios_base::floatfield);

#ifdef DUMP_TO_FILE
	ofstream TW_FILE("tw.txt");
	for (int i = 0; i < MAX_FFT_LEN; i++)
	{
		TW_FILE << caTw[i] << endl;
	}

	TW_FILE.close();
#endif
}

void fft::fnFftCalculate()
{
	int iNodePartNum;
	int iNodeNum;
	int iData0Idx;
	int iData1Idx;
	int iTwIdx;
	int iMappedTxIdx;
	int iSign;
	double fSin;
	double fCos;

	complex cData0(0.0, 0.0, DATA_TW, DATA_FW, ROUNDING);
	complex cData1(0.0, 0.0, DATA_TW, DATA_FW, ROUNDING);
	complex cData0PlusData1(0.0, 0.0, DATA_TW + 1, DATA_FW, ROUNDING);
	complex cData0MinusData1(0.0, 0.0, DATA_TW + 1, DATA_FW, ROUNDING);
	complex cData0MinusData1MulTw(0.0, 0.0, DATA_TW + TF_TW + 2, DATA_FW + TF_TW, ROUNDING);
	complex cTw(0.0, 0.0, TF_TW, TF_FW, ROUNDING);

	iNodePartNum = 1;

	for (int i = 0; i < iTotalNumOfStage; i++)
	{
		//cout << "Perform calculation on stage... " << i << endl;

		iNodeNum = iTotalNumOfNode / iNodePartNum;
		for (int j = 0; j < iNodePartNum; j++)
		{
			for (int k = 0; k < iNodeNum; k++)
			{
				iData0Idx = j * (iNodeNum * 2) + k;
				iData1Idx = iData0Idx + iNodeNum;
				iTwIdx = k * iNodePartNum*MAX_FFT_LEN / iFftLength;

				cData0 = caStage[i][iData0Idx];
				cData1 = caStage[i][iData1Idx];

				if (i == 0 && bInversedFFt)
				{
					cData0.conjugate();
					cData1.conjugate();
				}

				cData0PlusData1 = cData0 + cData1;
				cData0MinusData1 = cData0 - cData1;

				if (iTwIdx <= 2048)
				{
					iMappedTxIdx = iTwIdx;
					iSign = 1;
				}
				else
				{
					iMappedTxIdx = 2048 - (iTwIdx - 2048);
					iSign = -1;
				}
				
				fCos = caTw[iMappedTxIdx].getReal();
				fSin = caTw[iMappedTxIdx].getImag();
				fCos = iSign * fCos;

				cTw.set(fCos, fSin, TF_TW, TF_FW, ROUNDING);

				cData0MinusData1MulTw = cData0MinusData1 * cTw;

				if (bInversedFFt)
				{
					cData0PlusData1 = cData0PlusData1 * 0.5;
					cData0MinusData1MulTw = cData0MinusData1MulTw * 0.5;
				}

				if (i == iTotalNumOfStage - 1 && bInversedFFt)
				{
					cData0PlusData1.conjugate();
					cData0MinusData1MulTw.conjugate();
				}

				caStage[i + 1][iData0Idx] = cData0PlusData1;
				caStage[i + 1][iData1Idx] = cData0MinusData1MulTw;
			}
		}
		iNodePartNum = iNodePartNum * 2;
	}
}

void fft::fnReorder()
{
	int iBw = int(log(iFftLength) / log(2));
	int iGain;
	int iFlag;
	int iReversed;

	for (int i = 0; i < iTotalNumOfNode * 2; i++)
	{
		caStage[iTotalNumOfStage + 1][i] = caStage[iTotalNumOfStage][i];
	}

	for (int i = 0; i < iTotalNumOfNode*2; i++)
	{
		iGain = int(pow(2, iBw - 1));
		iFlag = 0x00000001;
		iReversed = 0;
		for (int j = 0; j < iBw; j++)
		{
			iReversed = iReversed + iGain * ((i&iFlag) / iFlag);
			iGain = iGain / 2;
			iFlag = iFlag * 2;
		}

		if (i != iReversed)
		{
			caStage[iTotalNumOfStage + 1][iReversed] = caStage[iTotalNumOfStage][i];
		}
	}

#ifdef DUMP_TO_FILE
	char sFileName[512];

	for (int i = 0; i < iTotalNumOfStage + 2; i++)
	{
		if (i == 0)
		{
#ifdef MICROSOFT_C
			sprintf_s(sFileName, "input.txt");
#else
			sprintf(sFileName, "input.txt");
#endif
		}
		else if (i < iTotalNumOfStage + 1)
		{
#ifdef MICROSOFT_C
			sprintf_s(sFileName, "stage_%0d.txt", i);
#else
			sprintf(sFileName, "stage_%0d.txt", i);
#endif
		}
		else
		{
#ifdef MICROSOFT_C
			sprintf_s(sFileName, "output.txt");
#else
			sprintf(sFileName, "output.txt");
#endif
		}

		ofstream STAGE_FILE(sFileName);

		for (int j = 0; j < iFftLength; j++)
		{
			STAGE_FILE << caStage[i][j] << endl;
		}

		STAGE_FILE.close();
	}
#endif
}


void fft::fnReorder(int *poutdata)
{
	int iBw = int(log(iFftLength) / log(2));
	int iGain;
	int iFlag;
	int iReversed;

	for (int i = 0; i < iTotalNumOfNode * 2; i++)
	{
		caStage[iTotalNumOfStage + 1][i] = caStage[iTotalNumOfStage][i];
	}

	for (int i = 0; i < iTotalNumOfNode*2; i++)
	{
		iGain = int(pow(2, iBw - 1));
		iFlag = 0x00000001;
		iReversed = 0;
		for (int j = 0; j < iBw; j++)
		{
			iReversed = iReversed + iGain * ((i&iFlag) / iFlag);
			iGain = iGain / 2;
			iFlag = iFlag * 2;
		}

		if (i != iReversed)
		{
			caStage[iTotalNumOfStage + 1][iReversed] = caStage[iTotalNumOfStage][i];
		}
	}

#if 0
#ifdef DUMP_TO_FILE
	char sFileName[512];

	for (int i = 0; i < iTotalNumOfStage + 2; i++)
	{
		if (i == 0)
		{
#ifdef MICROSOFT_C
			sprintf_s(sFileName, "input.txt");
#else
			sprintf(sFileName, "input.txt");
#endif
		}
		else if (i < iTotalNumOfStage + 1)
		{
#ifdef MICROSOFT_C
			sprintf_s(sFileName, "stage_%0d.txt", i);
#else
			sprintf(sFileName, "stage_%0d.txt", i);
#endif
		}
		else
		{
#ifdef MICROSOFT_C
			sprintf_s(sFileName, "output.txt");
#else
			sprintf(sFileName, "output.txt");
#endif
		}

		ofstream STAGE_FILE(sFileName);

		for (int j = 0; j < iFftLength; j++)
		{
			STAGE_FILE << caStage[i][j] << endl;
		}

		STAGE_FILE.close();
	}
#endif
#endif
	for(int j=0; j<iFftLength; j++)
	{
		poutdata[2*j] = caStage[iTotalNumOfStage + 1][j].getReal_Int();
		poutdata[2*j +1] = caStage[iTotalNumOfStage + 1][j].getImag_Int();
	}
}
