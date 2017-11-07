#include "defs.h"

#if defined(_WIN32) || defined(_WIN64)

char *basename(char *name){
	char *result = name;
	while(*name){
		if(*name == '/' || *name == '\\') result = name + 1;
		name++;
	}
	return result;
}
#endif

bool check_is_null(allData *aD, Ipp32f *v, const char *name, const char *comm){
	if(v != NULL)return false;
	sprintf(aD->message,"can not allocate memory (name \"%s\", address: %p, %s)",name, (void *)v, comm);
	printError(aD);
	return true;
}


bool find_power(int origw, int *next, int *extpower){
	int u, w;
	u = 0;
	w = 1;
	while(w < origw){
		u++;
		w*=2;
	}
	*next = w;
	*extpower = u;
	return true;
}


bool nom(int a, int b, int *c){
	int n, m;
	n = 1;
	m = 2;
	while(m< a && m < b){
		while(a%m == 0 && b%m == 0){
			n*=m;
			a/=m;
			b/=m;
		}
		m++;
	}
	*c = n*a*b;
	return true;
}


void int_ratio(unsigned int m, unsigned int d, unsigned int *n){
	unsigned int r, u;
	r = m%d;
	u = m/d;
	if(r >0) u++;
	*n = d*u;
}

int cudaRound(int v, int d){
	int r, q;
	r = v%d;
	q = v/d;
	if(r>0)q++;
	q*=d;
	return q;
}

bool toLower(char *value){
	int len;
	len = (int)strlen(value);
	ippsLowercaseLatin_8u_I((Ipp8u*)value, len);
	return true;
}

bool toUpper(char *value){
	int len;
	len = (int)strlen(value);
	ippsUppercaseLatin_8u_I((Ipp8u*)value, len);
	return true;
}

bool replaceSlash(char *value){
	int len;
	len = (int)strlen(value);
	for(int i=0;i<len;i++){
		if(value[i] == '\\')value[i] = '/';
	}
	return true;
}


void formFilterRL2(Ipp32f *vec,Ipp32f omega, Ipp32f ds, int nx){
	Ipp32f *vt = ippsMalloc_32f(nx+1);
	Ipp32f s, ss, pi2;
	pi2= float(pi*pi);
	for(int i=0;i<nx+1;i++){
		s=omega*ds*((float)(i));
		if(fabs(s)<1.e-5){
			ss= s*(0.5f-s*s/24.0f);
		}else{
			ss= (1.0f-cos(s))/s;
		}
		vt[i]=ss;
	}
	ippsCopy_32f(vt+1,vec,nx);
	ippsSub_32f_I(vt,vec,nx);
	Ipp32f val = omega/(4.0f*pi2*float(nx));
	ippsMulC_32f_I(val,vec,nx);
	ippsFlip_32f(vec,vec+nx/2,nx/2);
	ippsFree(vt); vt = NULL;
}


Ipp32f sinc(Ipp32f x){
	if(fabs(x)<1.e-5){
		Ipp32f y=x*x;
		return (float)(1.-y/6.0+y*y/120.0);
	}else{
		return (sin(x)/x);
	}
}
