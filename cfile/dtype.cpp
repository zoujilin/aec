#include <math.h>
#include "dtype.h"

complex::complex()
{
	this->fReal = 0.0;
	this->fImag = 0.0;
	this->iTotalWidth = 16;
	this->iFracWidth = 8;
	this->bRnd = false;

	this->quantize();
}

complex::complex(int i, int j,int tw, int fw,bool rnd)
{
	this->fReal = double(i);
	this->fImag = double(j);
	this->iTotalWidth = tw;
	this->iFracWidth = fw;
	this->bRnd = rnd;

	this->quantize();
}

complex::complex(float i, float j,int tw, int fw,bool rnd)
{
	this->fReal = double(i);
	this->fImag = double(j);
	this->iTotalWidth = tw;
	this->iFracWidth = fw;
	this->bRnd = rnd;

	this->quantize();
}

complex::complex(double i, double j,int tw, int fw,bool rnd)
{
	this->fReal = i;
	this->fImag = j;
	this->iTotalWidth = tw;
	this->iFracWidth = fw;
	this->bRnd = rnd;

	this->quantize();
}

complex::~complex()
{

}

void complex::set(int i, int j,int tw, int fw,bool rnd)
{
	this->fReal = double(i);
	this->fImag = double(j);
	this->iTotalWidth = tw;
	this->iFracWidth = fw;
	this->bRnd = rnd;

	this->quantize();
}

void complex::set(float i, float j,int tw, int fw,bool rnd)
{
	this->fReal = double(i);
	this->fImag = double(j);
	this->iTotalWidth = tw;
	this->iFracWidth = fw;
	this->bRnd = rnd;

	this->quantize();
}

void complex::set(double i, double j,int tw, int fw,bool rnd)
{
	this->fReal = i;
	this->fImag = j;
	this->iTotalWidth = tw;
	this->iFracWidth = fw;
	this->bRnd = rnd;

	this->quantize();
}

int complex::getTotalWidth()
{
	return iTotalWidth;
}

int complex::getFracWidth()
{
	return iFracWidth;
}

double complex::getReal()
{
	return fReal;
}

double complex::getImag()
{
	return fImag;
}

int complex::getImag_Int()
{
	return (int)(fImag * pow(2.0, iFracWidth));
}

int complex::getReal_Int()
{
	return (int)(fReal * pow(2.0, iFracWidth));
}

void complex::conjugate()
{
	fImag = -fImag;
	quantize();
}

void complex::quantize()
{
	double fRealGain;
	double fImagGain;
	double fUpLimit;
	double fLowLimit;

	fUpLimit = pow(2.0,this->iTotalWidth-1)-1;
	fLowLimit = pow(2.0,this->iTotalWidth-1)*(-1);

	//Quantize real part
	fRealGain = this->fReal * pow(2.0,this->iFracWidth);

	if(this->bRnd) {
		fRealGain = fRealGain + 0.5;
	}

   fRealGain = floor(fRealGain);

	if(fRealGain > fUpLimit)
	{
		fRealGain = fUpLimit;
	}
	else if(fRealGain < fLowLimit)
	{
		fRealGain = fLowLimit;
	}

	this->fReal = fRealGain/pow(2.0,this->iFracWidth);

	//Quantize imaginary part
	fImagGain = this->fImag * pow(2.0,this->iFracWidth);

	if(this->bRnd) {
		fImagGain = fImagGain + 0.5;
	}

	fImagGain = floor(fImagGain);

	if(fImagGain > fUpLimit)
	{
		fImagGain = fUpLimit;
	}
	else if(fImagGain < fLowLimit)
	{
		fImagGain = fLowLimit;
	}

	this->fImag = fImagGain/pow(2.0,this->iFracWidth);
}

ostream & operator<<(ostream & output, const complex & a)
{
	output.precision(16);

    output << a.fReal * pow(2.0, a.iFracWidth)  << " " << a.fImag * pow(2.0, a.iFracWidth);

	output.precision(6);

	return output;
}

complex operator+(const complex & a, const complex & b)
{
	double fR;
	double fI;
	int iTw;
	int iFw;

	fR = a.fReal + b.fReal;
	fI = a.fImag + b.fImag;

	if(a.iTotalWidth >= b.iTotalWidth)
	{
		iTw = a.iTotalWidth+1;
	}
	else
	{
		iTw = b.iTotalWidth+1;
	}

	if(a.iFracWidth >= b.iFracWidth)
	{
		iFw = a.iFracWidth;
	}
	else
	{
		iFw = b.iFracWidth;
	}

	complex cResult(fR,fI, iTw, iFw, a.bRnd);

	return cResult;
}

complex operator-(const complex & a, const complex & b)
{
	double fR;
	double fI;
	int iTw;
	int iFw;

	fR = a.fReal - b.fReal;
	fI = a.fImag - b.fImag;

	if (a.iTotalWidth >= b.iTotalWidth)
	{
		iTw = a.iTotalWidth + 1;
	}
	else
	{
		iTw = b.iTotalWidth + 1;
	}

	if (a.iFracWidth >= b.iFracWidth)
	{
		iFw = a.iFracWidth;
	}
	else
	{
		iFw = b.iFracWidth;
	}

	complex cResult(fR, fI, iTw, iFw, a.bRnd);

	return cResult;
}

complex operator*(const complex & a, const complex & b)
{
	double fK0 = a.fReal * b.fReal;
	double fK1 = a.fReal * b.fImag;
	double fK2 = a.fImag * b.fReal;
	double fK3 = a.fImag * b.fImag;

	double fR = fK0 - fK3;
	double fI = fK1 + fK2;

	int iTw = a.iTotalWidth + b.iTotalWidth + 1;
	int iFw = a.iFracWidth + b.iFracWidth;

	complex cResult(fR, fI, iTw, iFw, a.bRnd);

	return cResult;
}

complex operator*(const complex & a, double b)
{
	double fR = a.fReal * b;
	double fI = a.fImag * b;

	int iTw = a.iTotalWidth;
	int iFw = a.iFracWidth;

	complex cResult(fR, fI, iTw, iFw, a.bRnd);

	return cResult;
}

complex complex::operator=(const complex & a)
{
	this->fReal = a.fReal;
	this->fImag = a.fImag;
	this->quantize();

	return * this;
}

