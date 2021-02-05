#ifndef __COMPLEX_H__
#define __COMPLEX_H__

#include <math.h>

//-----------------------------复数类--------------------------------
template<typename T>
class Complex{
public:
	Complex(T real=0,T imag=0){
		this->real=real;
		this->imag=imag;
	}
	Complex(const Complex& c){
		real=c.real;
		imag=c.imag;
	}
	~Complex(){}
public:
	inline Complex& operator=(const Complex& c){
		if(this==&c)return *this;
		real=c.real;
		imag=c.imag;
		return *this;
	}
	inline Complex operator+(const Complex& c){
		return Complex(real+c.real,imag+c.imag);
	}
	inline Complex operator-(const Complex& c){
		return Complex(real-c.real,imag-c.imag);
	}
	inline Complex operator*(const Complex& c){
		return Complex(real*c.real-imag*c.imag,real*c.imag+imag*c.real);
	}
	inline Complex operator/(const Complex& c){
		double d=c.real*c.real+c.imag*c.imag;
		return Complex((real*c.real+imag*c.imag)/d,(imag*c.real-real*c.imag)/d);
	}
	inline Complex& operator+=(const Complex &c){
		real+=c.real;
		imag+=c.imag;
		return *this;
	}
	inline Complex& operator-=(const Complex &c){
		real-=c.real;
		imag-=c.imag;
		return *this;
	}
	inline Complex& operator*=(const Complex &c){
		double r=real*c.real-imag*c.imag;
		imag=real*c.imag+imag*c.real;
		real=r;
		return *this;
	}
	inline Complex& operator/=(const Complex& c){
		double d=c.real*c.real+c.imag*c.imag;
		double r=(real*c.real+imag*c.imag)/d;
		imag=(imag*c.real-real*c.imag)/d;
		real=r;
		return *this;
	}
	inline bool operator==(const Complex& c){
		return (real==c.real&&imag==c.imag)
	}
	inline bool operator!=(const Complex& c){
		return (real!=c.real||imag!=c.imag);
	}
	inline friend Complex operator-(const Complex &c){
		return Complex(-c.real,-c.image);
	}
public:
	inline Complex& operator=(T real){
		this->real=real;imag=0;
		return *this;
	}
	inline Complex operator+(T real){
		return Complex(this->real+real,imag);
	}
	inline Complex operator-(T real){
		return Complex(this->real-real,imag);
	}
	inline Complex operator*(T d){
		return Complex(real*d,imag*d);
	}
	inline Complex operator/(T d){
		return Complex(real/d,imag/d);
	}
	inline Complex& operator+=(T real){
		this->real+=real;
		return *this;
	}
	inline Complex& operator-=(T real){
		this->real-=real;
		return *this;
	}
	inline Complex& operator*=(T d){
		real*=d;
		imag*=d;
		return *this;
	}
	inline Complex& operator/=(T d){
		real/=d;
		imag/=d;
		return *this;
	}
public:
	inline double mould()const{
		return sqrt(real*real+imag*imag);
	}
	inline double phase()const{
		return atan(imag/real);
	}

#ifdef COMPLEX_DEBUG
	void print()const{
		char cr=real>=0?' ':'-';
		char ci=imag>=0?'+':'-';
		printf("(%c%11.4e %c%11.4ei)",cr,abs(real),ci,abs(imag));
	}
	friend std::ostream& operator<<(std::ostream &userStr,const Complex &c){
		c.print();
		return userStr;
	}
#endif

public:
	T real,imag;
};



#endif
