#ifndef __FFT_H__
#define __FFT_H__

#include "Complex.h"
#include "MArray.h"

/******************************************
*			 傅里叶变换函数对象
*	可实现一维和二维傅里叶变换和反变换
*	长度（二维的宽和高）为2的整数次幂
//*****************************************/

//----------------------------------
typedef unsigned int UINT;
#define PI 3.141592653589793
//----------------------------------

//----------------傅里叶变换(FFT_fun)-----------------------
template<typename Complex>//---需要一个复数类参数
class FFT_fun{
public:
	template<typename T>//---T是可以赋值给Complex的类型
	void operator()(MArray<Complex>& X,const MArray<T>& x){
		//检查N
		UINT t,i,j,k,index;
		UINT wl=sizeof(UINT)*8-1;
		const UINT N=x.len;
		for(index=0,t=1;index<wl&&t!=N;++index,t<<=1);
		if(t!=N)return;
		//倒序
		int *a=new int[N];
		a[0]=0;
		for(i=1;i<=N/2;i*=2)
			for(int j=i;j<2*i;++j)
				a[j]=a[j-i]+N/2/i;
		X.reset(N);
		for(i=0;i<N;++i)
			X[i]=x[a[i]];
		//构造W
		Complex *W=new Complex[N/2];
		const double C=2*PI/N;
		for(i=0;i<N/2;++i){
			W[i].real=cos(C*i);
			W[i].imag=-sin(C*i);
		}
		//计算离散傅里叶变换
		Complex temp;
		for(i=1;i<N;i*=2)
			for(j=0;j<N;j+=2*i)
				for(k=j;k<j+i;++k){
					X[k+i]*=W[(k-j)*N/2/i];
					temp=X[k];
					X[k]=X[k]+X[k+i];
					X[k+i]=temp-X[k+i];
				}
		//释放内存
		delete[] a;
		delete[] W;
	}

	//-------------------二维变换----------------------
	template<typename T>//---T是可以赋值给Complex的类型
	void operator()(MArray2<Complex> &X,const MArray2<T> &x){
		//检查N，M
		UINT t,i,j,index_n,index_m;
		const UINT N=x.height(),M=x.width();
		const UINT wl=sizeof(UINT)*8-1;
		for(index_n=0,t=1;index_n<wl&&t!=N;++index_n,t<<=1);
		if(t!=N)return;	
		for(index_m=0,t=1;index_m<wl&&t!=M;++index_m,t<<=1);
		if(t!=M)return;
		X.reset(N,M);
		//倒序
		int *an=new int[N];
		int *am=(N==M)?an:new int[M];
		an[0]=0;
		for(i=1;i<=N/2;i*=2)
			for(j=i;j<2*i;++j)
				an[j]=an[j-i]+N/2/i;
		if(N!=M){
			am[0]=0;
			for(i=1;i<=M/2;i*=2)
				for(j=i;j<2*i;++j)
					am[j]=am[j-i]+M/2/i;
		}
		for(i=0;i<N;++i)
			for(j=0;j<M;++j)
				X[i][j]=x[an[i]][am[j]];
		//构造W
		const UINT half_N=N/2,half_M=M/2;
		Complex *Wn=new Complex[N/2];
		Complex *Wm=(N==M)?Wn:new Complex[M/2];
		Complex *Wnm=new Complex[half_N*half_M];
		const double Cn=2*PI/N;
		for(i=0;i<half_N;++i){
			Wn[i].real=cos(Cn*i);
			Wn[i].imag=-sin(Cn*i);
		}
		if(N!=M){
			const double Cm=2*PI/M;
			for(i=0;i<half_M;++i){
				Wm[i].real=cos(Cm*i);
				Wm[i].imag=-sin(Cm*i);
			}
		}
		for(i=0;i<half_N;++i)
			for(j=0;j<half_M;++j)
				Wnm[i*half_M+j]=Wn[i]*Wm[j];
		//计算傅里叶变换
		UINT *in,*im;
		UINT min;
		if(N<M){
			min=index_n;
			UINT len=M/N;
			in=new UINT[min];
			im=new UINT[min];
			in[0]=1;im[0]=len;
			for(i=1;i<min;++i){
				in[i]=in[i-1]*2;
				im[i]=M*in[i]/N;
			}
			Complex *W=new Complex[len/2]; 
			Complex temp;
			const double C=2*PI/len;
			for(i=0;i<len/2;++i){
				W[i].real=cos(C*i);
				W[i].imag=-sin(C*i);
			}
			UINT n,m,k;
			for(n=0;n<N;++n)
				for(m=0;m<M;m+=len)
					for(i=1;i<len;i*=2)
						for(j=m;j<m+len;j+=2*i)
							for(k=j;k<j+i;++k){
								X[n][k+i]*=W[(k-j)*len/2/i];
								temp=X[n][k+i];
								(X[n][k+i]-=(temp*2))+=X[n][k];
								X[n][k]+=temp;
							}
			delete[] W;
		}else if(M<N){
			min=index_m;
			UINT len=N/M;
			in=new UINT[min];
			im=new UINT[min];
			im[0]=1;in[0]=len;
			for(i=1;i<min;++i){
				im[i]=im[i-1]*2;
				in[i]=N*im[i]/M;
			}
			Complex *W=new Complex[len/2]; 
			Complex temp;
			const double C=2*PI/len;
			for(i=0;i<len/2;++i){
				W[i].real=cos(C*i);
				W[i].imag=-sin(C*i);
			}
			UINT n,m,k;
			for(n=0;n<N;n+=len)
				for(m=0;m<M;++m)
					for(i=1;i<len;i*=2)
						for(j=n;j<n+len;j+=2*i)
							for(k=j;k<j+i;++k){
								X[k+i][m]*=W[(k-j)*len/2/i];
								temp=X[k+i][m];
								(X[k+i][m]-=(temp*2))+=X[k][m];
								X[k][m]+=temp;
							}
			delete[] W;	
		}else{
			min=index_n;
			in=new UINT[min];
			im=in;in[0]=1;
			for(i=1;i<min;++i)
				in[i]=in[i-1]*2;
		}
		UINT n,m,ni,mi;
		UINT jn,jm,kn,km;
		Complex temp1,temp2,temp3,temp4;
		for(i=0;i<min;++i){
			ni=in[i];mi=im[i];
			for(jn=0;jn<N;jn+=2*ni)
				for(jm=0;jm<M;jm+=2*mi)
					for(kn=jn;kn<jn+ni;++kn){
						n=(kn-jn)*N/2/ni;
						for(km=jm;km<jm+mi;++km){
							m=(km-jm)*M/2/mi;
							X[kn+ni][km]*=Wn[n];
							X[kn][km+mi]*=Wm[m];
							X[kn+ni][km+mi]*=Wnm[n*half_M+m];
							temp1=X[kn][km];
							temp2=X[kn][km+mi];
							temp3=X[kn+ni][km];
							temp4=X[kn+ni][km+mi];
							((X[kn][km]+=temp2)+=temp3)+=temp4;
							(((X[kn][km+mi]-=temp2*2)+=temp1)+=temp3)-=temp4;
							(((X[kn+ni][km]-=temp3*2)+=temp1)+=temp2)-=temp4;
							((X[kn+ni][km+mi]+=temp1)-=temp2)-=temp3;
						}
					}
		}
		//释放内存
		delete[] an;
		delete[] Wn;
		delete[] Wnm;
		delete[] in;
		if(N!=M){
			delete[] am;
			delete[] Wm;
			delete[] im;
		}
	}

public:
	//计算大于等于n的2的幂次
	UINT get2n(UINT n){
		UINT N,i,wl=sizeof(UINT)*8-1;
		for(i=0,N=1;i<wl&&N<n;++i,N<<=1);
		return N;
	}
	//补零zeros-padding
	template<typename T1,typename T2>//--T2是可以给T1赋值的类型
	void zp(MArray<T1> out,const MArray<T2> &in){
		UINT wl=sizeof(int)*8-1;
		UINT i,N,n=in.length();
		for(i=0,N=1;i<wl&&N<n;++i,N<<=1);
		out.reset(N);
		for(i=0;i<n;++i)
			out[i]=in[i];
		for(i=n;i<N;++i)
			out[i]=0;
	}
	template<typename T1,typename T2>//--T2是可以给T1赋值的类型
	void zp(MArray2<T1>& out,const MArray2<T2>& in){
		const UINT n=in.height();
		const UINT m=in.width();
		UINT N,M,i,j;
		const UINT wl=sizeof(UINT)*8-1;
		for(i=0,N=1;i<wl-1&&N<n;++i,N<<=1);
		for(j=0,M=1;j<wl-1&&M<m;++j,M<<=1);
		out.reset(N,M);
		for(i=0;i<n;++i)
			for(j=0;j<m;++j)
				out[i][j]=in[i][j];
		if(m!=M){
			for(i=0;i<n;++i)
				for(j=m;j<M;++j)
					out[i][j]=0;
		}
		if(n!=N){
			for(i=n;i<N;++i)
				for(j=0;j<M;++j)
					out[i][j]=0;
		}
	}
};
//-------------//end(FFT_fun)------------------------//*/

//--------------------反变换(IFFT_fun)------------------------
template<typename Complex>
class IFFT_fun{
public:
	template<typename T1,typename T2>//---x接收反傅里叶变换实部的值
	void operator()(MArray<T1> &x,const MArray<T2> &X){
		UINT i,N=X.length();
		MArray<Complex> xi(N);
		FFT_fun<Complex> fft;
		fft(xi,X,N);
		for(i=0;i<N;++i)
			x[i]=xi[i].real/N;
		T1 temp;
		for(i=1;i<N/2;i++){
			temp=x[i];
			x[i]=x[N-i];
			x[N-i]=temp;
		}
	}

	template<typename T>
	void operator()(MArray<Complex> &x,const MArray<T> &X){
		UINT i,N=X.length();
		FFT_fun<Complex> fft;
		fft(x,X,N);
		x[0]/=N;x[N/2]/=N;
		Complex temp;
		for(i=1;i<N/2;++i){
			temp=x[i]/N;
			x[i]=x[N-i]/N;
			x[N-i]=temp;
		}
	}

	//--------------------二维变换-----------------
	template<typename T1,typename T2>//---x接收反傅里叶变换实部的值
	void operator()(MArray2<T1> &x,const MArray2<T2> &X){
		const UINT N=X.height();
		const UINT M=X.width();
		const UINT L=N*M;
		UINT i,j;
		T1 temp;
		MArray2<Complex> xi(N,M);
		FFT_fun<Complex> fft;
		fft(xi,X);
		x.reset(N,M);
		for(i=0;i<N;++i)
			for(j=0;j<M;++j)
				x[i][j]=(xi[i][j].real/L);
		for(i=1;i<N;++i)
			for(j=1;j<M/2;++j){
				temp=x[i][j];
				x[i][j]=x[N-i][M-j];
				x[N-i][M-j]=temp;
			}
		for(j=1;j<M/2;++j){
			temp=x[0][j];
			x[0][j]=x[0][M-j];
			x[0][M-j]=temp;
		}
		for(j=0;j<=M/2;j+=M/2)
			for(i=1;i<N/2;++i){
				temp=x[i][j];
				x[i][j]=x[N-i][j];
				x[N-i][j]=temp;
			}
	}

	template<typename T>
	void operator()(MArray2<Complex> &x,const MArray2<T> &X){
		FFT_fun<Complex> fft;
		const UINT N=X.height();
		const UINT M=X.width();
		const UINT L=N*M;
		UINT i,j;
		Complex temp;
		fft(x,X);
		for(i=1;i<N;++i)
			for(j=1;j<M/2;++j){
				temp=x[i][j]/L;
				x[i][j]=x[N-i][M-j]/L;
				x[N-i][M-j]=temp;
			}
		for(j=1;j<M/2;++j){
			temp=x[0][j]/L;
			x[0][j]=x[0][M-j]/L;
			x[0][M-j]=temp;
		}
		for(j=0;j<=M/2;j+=M/2)
			for(i=1;i<N/2;++i){
				temp=x[i][j]/L;
				x[i][j]=x[N-i][j]/L;
				x[N-i][j]=temp;
			}
		for(i=0;i<N;i+=N/2)
			for(j=0;j<M;j+=M/2)
				x[i][j]/=L;
	}\
};
//------------------//end(IFFT_fun)------------------------//*/

#endif