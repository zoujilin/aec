#ifndef DTYPE_H_
#define DTYPE_H_

#include <iostream>
using namespace std;
class complex {

private:
	double fReal;
	double fImag;
	bool bRnd;
	int iTotalWidth;
	int iFracWidth;

public:
	//Constructor & De-constructor
	complex();
	complex(int i, int j, int tw, int fw, bool rnd=false);
	complex(float i, float j, int tw, int fw, bool rnd=false);
	complex(double i, double j, int tw, int fw, bool rnd=false);
	~complex();

	//Quantize
	void quantize();

	//Setter
	void set(int i, int j, int tw, int fw, bool rnd=false);
	void set(float i, float j, int tw, int fw, bool rnd=false);
	void set(double i, double j, int tw, int fw, bool rnd=false);

	//Getter
	int getTotalWidth();
	int getFracWidth();
	double getReal();
	double getImag();
	int getReal_Int();
	int getImag_Int();

	//Conjugation
	void conjugate();

	//Operator overload
	friend ostream & operator<<(ostream & output, const complex & a);
	friend complex operator+(const complex & a, const complex & b);
	friend complex operator-(const complex & a, const complex & b);
	friend complex operator*(const complex & a, const complex & b);
	friend complex operator*(const complex & a, double b);
	complex operator=(const complex & a);
};

#endif /* DTYPE_H_ */

