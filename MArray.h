#ifndef _MARRAY2_H__
#define _MARRAY2_H__

//-----------------------------
typedef unsigned int UINT;
//-----------------------------
//-------------------数组类(MArray)---------------------------
template<typename T>
class MArray{
public:
	MArray(){
		len=0;
		pA=NULL;
	}
	MArray(UINT len){
		if(len){
			this->len=len;
			pA=new T[len];
		}else{
			this->len=0;
			pA=NULL;
		}
	}
	MArray(const MArray<T>& arr){
		len=arr.len;
		if(len){
			pA=new T[len];
			for(int i=0;i<len;++i)
				pA[i]=arr.pA[i];
		}
		else pA=NULL;
	}
	~MArray(){
		if(pA!=NULL)delete[] pA;
	}
public:
	MArray& operator=(const MArray<T>& arr){
		if(this==&arr)return *this;
		reset(arr.len);
		save(arr.pA);
		return *this;
	}
	MArray& operator+=(const MArray<T>& arr){
		for(UINT i=0;i<len;++i)
			pA[i]+=arr.pA[i];
		return *this;
	}
	MArray& operator-=(const MArray<T>& arr){
		for(UINT i=0;i<len;++i)
			pA[i]-=arr.pA[i];
		return *this;
	}
	MArray& operator*=(const MArray<T>& arr){
		for(UINT i=0;i<len;++i)
			pA[i]*=arr.pA[i];
		return *this;
	}
	inline T& operator[](int n)const{
		return *(pA+n);
	}
	inline bool operator==(MArray<T> &arr)const{return arr.len==len;}
	inline bool operator!=(MArray<T> &arr)const{return arr.len!=len;}
public:
	void set(UINT len){
		if(len==this->len)return;
		if(len){
			T* p=new T[len];
			int min=(len<this->len)?len:this->len;
			for(int i=0;i<min;++i)
				p[i]=pA[i];
			this->len=len;
			if(pA!=NULL)delete[] pA;
			pA=p;
		}else{
			this->len=0;
			if(pA!=NULL)delete[] pA;
			pA=NULL;
		}
	}
	void reset(UINT len=0){
		if(len==this->len)return;
		if(len){
			this->len=len;
			if(pA!=NULL)delete[] pA;
			pA=new T[len];
		}else{
			this->len=0;
			if(pA!=NULL)delete[] pA;
			pA=NULL;
		}
	}
	template<typename pT>void save(pT p){
		for(UINT i=0;i<len;++i,++p)
			pA[i]=*p;
	}
	template<typename pT>void save(pT p,int len){
		int min=(len<this->len)?len:this->len;
		for(UINT i=0;i<min;++i,++p)
			pA[i]=*p;
	}
	template<typename pT>void save(pT frist,pT last){
		for(UINT i=0;frist!=last&&i<len;++frist,++i)
			pA[i]=*frist;
	}
public:
	inline UINT length()const{return len;}
	inline T* getPA()const{return pA;}
public:
	template<typename Oper_fun>void arrTraverse(Oper_fun oper){
		for(UINT i=0;i<len;++i)
			oper(pA[i]);
	}
protected:
	UINT len;
	T *pA;
};
//-----------------//end(MArray)-----------------------------//*/

//------------------二维数组类(MArray2)----------------------//
template<typename T>
class MArray2:public MArray<T>{
public:
	MArray2(){
		N=0;M=0;
	}
	MArray2(UINT N,UINT M):MArray<T>(N*M){
		if(N<=0||M<=0){
			this->N=0;this->M=0;
		}else{
			this->N=N;this->M=M;
		}
	}
	MArray2(const MArray2<T> &arr):MArray<T>(arr.len){
		N=arr.N;M=arr.M;
		save(arr.pA);
	}
	~MArray2(){}
public:
	inline T* operator[](int i)const{
		return pA+i*M;
	}
	inline T& operator()(int i,int j)const{
		return *(pA+i*M+j);
	}
	inline T& operator()(int n)const{
		return *(pA+n);
	}
	MArray2<T>& operator=(const MArray2<T> &arr){
		if(this==&arr)return *this;
		reset(arr.N,arr.M);
		save(arr.pA);
		return *this;
	}
	MArray2<T>& operator+=(const MArray2<T>& ma){
		for(UINT i=0;i<len;++i)
			pA[i]+=ma.pA[i];
		return *this;
	}
	MArray2<T>& operator*=(const MArray2<T>& arr){
		for(UINT i=0;i<len;++i)
			pA[i]*=arr.pA[i];
		return *this;
	}
	inline bool operator==(const MArray2<T>& arr)const{
		return arr.N==N&&arr.M==M;
	}
	inline bool operator!=(const MArray2<T>& arr)const{
		return arr.N!=N||arr.M!=M;
	}
public:
	void set(UINT len){
		if(len){
			N=1;M=len;
		}else{
			N=0;M=0;
		}
		MArray<T>::set(len);
	}
	void set(UINT N,UINT M){
		if(N&&M){
			if(M==this->M){
				if(N==this->N)return;
				else{
					this->N=N;
					MArray<T>::set(N*M);	
				}
			}else{
				len=N*M;
				T *p=new T[len];
				UINT minN=(N<this->N)?N:this->N;
				UINT minM=(M<this->M)?M:this->M;
				UINT i,j;
				for(i=0;i<minN;++i)
					for(j=0;j<minM;++j)
						p[i*M+j]=pA[i*this->M+j];
				if(pA!=NULL)delete[] pA;
				pA=p;
				this->M=M;this->N=N;
			}
		}
		else{
			this->N=0;this->M=0;
			if(pA!=NULL)delete[] pA;
			len=0;
			pA=NULL;
		}
	}
	void reset(UINT N=0,UINT M=0){
		if(N&&M){
			if(N*M!=len){
				delete[] pA;
				len=N*M;
				pA=new T[len];
			}
			this->N=N;this->M=M;
		}
		else{
			this->N=0;this->M=0;
			if(pA!=NULL)delete[] pA;
			len=0;
			pA=NULL;
		}
	}
public:
	inline UINT width()const{return M;}
	inline UINT height()const{return N;}
protected:
	UINT N,M;
};
//-------------------//end(MArray2)-------------------//*/

#endif
