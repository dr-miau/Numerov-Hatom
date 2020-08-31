#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double potential(double E, float Z, int l, double r){
	double Ar;
	Ar=2.0*E+2.0*Z/r-l*(l+1.0)/(r*r);
	return Ar;
}

void numerov(int N, double ener, double r0, double h,double *f, float Z, int l,float rcut){	
	double A[3],num,denom,pos;
	int k,i,Ai;
	for(i=1;i<N;i++){
		k=0;
		while(k<3){
			pos=r0+h*(i+k-1);
			Ai=i+k-1;
			if(pos<rcut){
				A[k]=potential(ener,Z,l,r0+h*(Ai));
//				printf("%d %d %d %d\n",Ai,i,k,N);
			}
			else{
				A[k]=potential(ener,0,l,r0+h*(Ai));
			}
			k++;
		}
//		A1=potential(ener,Z,l,r0+h*i);
//		A2=potential(ener,Z,l,r0+h*(i+1));
		if(i==1){
			num=(2-(h*h*5*A[1]/6))*f[i];
		}
		else{
			num=(2-(h*h*5*A[1]/6))*f[i]-(1+(h*h*A[0]/12))*f[i-1];
		}
		denom=1+(h*h*A[2]/12);
		f[i+1]=num/denom;
	}
//		printf("%f\n",ener);
}

void normalize(double *f, double h, int N){
	double Integ,Aint;
	for(int j=0;j<N;j++){
		Integ=Integ+(h*0.5)*(f[j]*f[j]+f[j+1]*f[j+1]);
	}
	Aint=sqrt(1/Integ);
	for(int j=0;j<N+1;j++){
		f[j]=f[j]*Aint;
	}
}

int main(){
	int N,nstate,l;
	double h,f0=0,r0=0,f1,ff=0,rf=20,enero,enern,epos,eneg,ebiz,enerh,enerf,*f,ffo,ffn,sh;
	FILE *fp1;
	double tol=0.0001;
	float rcut,Z;

	scanf("%d",&N);
	scanf("%d",&l);
	scanf("%f",&Z);
	scanf("%d",&nstate);
	scanf("%f",&rcut);
//	rcut=2.0;
//	scanf("%f",&ener);
	enero=-1.0;
	enern=0;
	enerh=0.1;	
	h=(rf-r0)/N;
	f1=h;
	
	f=(double*)malloc((N+1)*sizeof(double));
	f[0]=f0;
	f[1]=f1;

//ciclo sobre estados
	for(int i=0;i<nstate;i++){
		if(i>0){
			enero=enerf+tol;
		}
	numerov(N,enero,r0,h,f,Z,l,rcut);
	ffo=f[N];

//	printf("%f\n",fabs(f[N]));	
	if(fabs(f[N])>tol){
	if(ffo>0){
		enern=enern+enerh;
		numerov(N,enern,r0,h,f,Z,l,rcut);
		ffn=f[N];
		epos=enero;
		while(ffn/ffo>0){
			enern=enern+enerh;
			numerov(N,enern,r0,h,f,Z,l,rcut);
			ffn=f[N];	
		}
		eneg=enern;
		while(fabs(f[N])>tol){
			ebiz=(epos+eneg)/2;
			enern=ebiz;
			numerov(N,ebiz,r0,h,f,Z,l,rcut);
			ffn=f[N];
			if(ffn>0){
				epos=ebiz;
			}
			if(ffn<0){
				eneg=ebiz;
			}
		}
		enerf=enern;
	}

	if(ffo<0){
		enern=enern+enerh;
		numerov(N,enern,r0,h,f,Z,l,rcut);
		ffn=f[N];
		eneg=enero;
		while(ffn/ffo>0){
			enern=enern+enerh;
			numerov(N,enern,r0,h,f,Z,l,rcut);
			ffn=f[N];
		}
		eneg=enern;
		while(fabs(f[N])>tol){
			ebiz=(epos+eneg)/2;
			enern=ebiz;
			numerov(N,ebiz,r0,h,f,Z,l,rcut);
			ffn=f[N];
			if(ffn>0){
				eneg=ebiz;
			}
			if(ffn<0){
				epos=ebiz;
			}
		}
		enerf=enern;
	}
	}

	normalize(f,h,N);

	printf("%lf\n\n",enerf);
	for(int i=0;i<N+1;i++){
	printf("%f %f\n",(r0+i*h),f[i]);
	}
	printf("\n\n");
	}
/*
	fp1=fopen("edo1.out","w");

        for(int i=0;i<N+1;i++){
                fprintf(fp1,"%f\n",f[i]);
        }
        fclose(fp1);
*/
	free(f);
	return 0;
}
