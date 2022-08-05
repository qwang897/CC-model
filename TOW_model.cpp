#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>

#define Fig5					0			//Fig5, FigS2 and FigS3		
#define FigS4					1			//1,2,3,4 means FigS4A,B,C,D

#define snap 					50000
#define print_parameters		0


int main() {
	int 	bind1=0,bind2=0,num=0,i=0,j=0,T=0,K=0,D=0,N=0,k=0,K1=0,imax=50,jmax=50,times=5;
	int 	n1=0,n2=0,n3=0,n4=0,n5=0,n6=0;
	double x1[100000]={0},x2[100000]={0},x3[100000]={0},x4[100000]={0},x5[100000]={0},x6[100000]={0};
	double s1=0,s2=0,pi1=0,v10=0,e10=0,e1=0;
	double pi2=0,v20=0,e20=0,e2=0;
	double r1=0,l1=0,l2=0,l3=0,l4=0,ltot=0,T0=0,K0=0,D0=0;
	double smax1[100000]={0},smax2[100000]={0},stot1=0,stot2=0;
	double l0=8,s0=0,k1,k2,F1=0,F2=0,vF1=0,vF2=0,Fs1,Fs2,Fd1,Fd2;
	double s_Fd2[10000]={0},s_1[10000]={0},s_2[10000]={0},s_3[10000]={0};
	double p1=0,p2=0,p3=0,p4=0,p5=0,p6=0,x,y,xy,mid;
	double S1,S2,S3,S4,S5,S6;
	double r02,r03,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15;
	double v1back,v2back,change_num;
	FILE* k_d1_d2=fopen("./output.txt","w");
	
	if(Fig5==1)			kmax=10000;
	if(FigS4==1)		kmax=100;

	srand((int)time(NULL));
	for(k=1;k<=kmax;k++){
		if(Fig5==1){
			r02=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r03=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r4 =(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r5 =(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r6 =(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r7 =(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r8 =(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r9 =(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r10=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r11=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r12=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r13=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r14=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
			r15=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;

			if(r02>=0)		v10=1000*(r02+1);
			else 			v10=1000/(-r02+1);
			if(r03>=0)		v20=400*(r03+1);
			else 			v20=400/(-r03+1);
			if(r4>=0)		e10=1*(r4+1);
			else 			e10=1/(-r4+1);
			if(r5>=0)		e20=0.8*(r5+1);
			else 			e20=0.8/(-r5+1);
			if(r6>=0)		Fs1=6*(r6+1);
			else 			Fs1=6/(-r6+1);
			if(r7>=0)		Fs2=1.1*(r7+1);
			else 			Fs2=1.1/(-r7+1);
			if(r8>=0)		Fd1=3*(r8+1);
			else 			Fd1=3/(-r8+1);
			if(r9>=0)		Fd2=0.75*(r9+1);
			else 			Fd2=0.75/(-r9+1);
			if(r10>=0)		k1 =0.2*(r10+1);
			else 			k1 =0.2/(-r10+1);
			if(r11>=0)		k2 =0.2*(r11+1);
			else 			k2 =0.2/(-r11+1);
			if(r12>=0)		pi1=10*(r12+1);
			else 			pi1=10/(-r12+1);
			if(r13>=0)		pi2=1.6*(r13+1);
			else 			pi2=1.6/(-r13+1);
			if(r14>=0)		v1back=72*(r14+1);
			else 			v1back=72/(-r14+1);
			if(r15>=0)		v2back=67*(r15+1);
			else		 	v2back=67/(-r15+1);
		}
		else{
			Fs1 = 6; 		Fs2 = 1.1;
			Fd1 = 3;    	Fd2 = 0.75;
			pi1 = 5;     	pi2 = 1.6;
			v10 = 1000;    	v20 = 1;
			v1back = 6;   	v2back = 72;
			e10 = 1;		e20 = 0.27;
			k1=0.2;			k2=0.05;
		}

		switch(FigS4){
			case 1:
				{change_num=k*0.001;k1=change_num;}
			break;
			case 2:
				{change_num=k*20;v10=change_num;}
			break;
			case 3:
				{change_num=k*0.1;pi1=change_num;}
			break;
			case 4:
				{change_num=k*0.01;e20=change_num;}
			break;
		}

		for(j=1;j<=1;j++){
			ltot=0;bind1=1;bind2=1;s1=0;s2=0;num=0;n1=0;n2=0;n3=0;n4=0;n5=0;n6=0;
			for(i=0;i<=100000;i++)					{x1[i]=0;x2[i]=0;x3[i]=0;x4[i]=0;x5[i]=0;x6[i]=0;}
			T=0;K=0;D=0;N=0;T0=0;D0=0;K0=0;
			//if(j==0)			pi2=0;				else				pi2=1.6;
			for(i=0;i<snap+10;i++)					{smax1[i]=0;smax2[i]=0;}
			while(1){
				if(bind1==0 && bind2==0)							s0=0;
				else												s0=(k1*s1+k2*s2)/(bind1*k1+bind2*k2);
				smax1[num]=s0;
				smax2[num]=s1;
				r1=rand()/(double)(RAND_MAX);
				if(bind1==0 && bind2==0){
					T=0;K=0;D=0;N=1;
					l1=pi1;		l2=0;		l3=pi2;		l4=0;		ltot=l1+l2+l3+l4;
					if(r1>=0 && r1<l1/ltot)							{s1++;bind1=1;}
					else if(r1>=(l1+l2)/ltot && r1<=(l1+l2+l3)/ltot){s2--;bind2=1;}
					}
				else if(bind1==1 && bind2==0){
					if(K==0)										K0=s1;
					T=0;K=1;D=0;N=0;
					l1=v10/l0;	l2=e10;	l3=pi2;		l4=0;			ltot=l1+l2+l3+l4;
					if(r1>=0 && r1<l1/ltot)							{s1++;}
					else if(r1>=l1/ltot && r1<(l1+l2)/ltot)			{bind1=0;x3[n3]=s1-K0;n3++;s1=0;}
					else if(r1>=(l1+l2)/ltot && r1<=(l1+l2+l3)/ltot){bind2=1;x4[n4]=s1-K0;n4++;s2=s1;}

					}
				else if(bind1==0 && bind2==1){
					if(D==0)										D0=s0;
					T=0;K=0;D=1;N=0;
					l1=pi1;		l2=0;		l3=v20/l0;	l4=e20;		ltot=l1+l2+l3+l4;
					if(r1>=0 && r1<l1/ltot)							{bind1=1;x5[n5]=s2-D0;n5++;s1=s2;}
					else if(r1>=(l1+l2)/ltot && r1<(l1+l2+l3)/ltot)	{s2--;}
					else if(r1>=(l1+l2+l3)/ltot && r1<=1)			{bind2=0;x6[n6]=s2-D0;n6++;s2=0;}
					}
				else if(bind1==1 && bind2==1){
					if(T==0)										T0=s0;
					T=1;K=0;D=0;N=0;
					s0=(k1*s1+k2*s2)/(bind1*k1+bind2*k2);
					F1=k1*l0*(s1-s0);								F2=k2*l0*(s0-s2);
					vF1=v10*(1-F1/Fs1);								vF2=v20*(1-F2/Fs2);
					if(F1>Fs1)										vF1=0;
					if(F2>Fs2)										vF2=0;
					e1=e10*exp(fabs(F1)/Fd1);						e2=e20*exp(fabs(F2)/Fd2);
					l1=vF1/l0;	l2=e1;		l3=vF2/l0;		l4=e2;	ltot=l1+l2+l3+l4;
					if(r1>=0 && r1<l1/ltot)							{s1++;}
					else if(r1>=l1/ltot && r1<(l1+l2)/ltot)			{bind1=0;x2[n2]=s2-T0;n2++;s1=0;}
					else if(r1>=(l1+l2)/ltot && r1<(l1+l2+l3)/ltot)	{s2--;}
					else if(r1>=(l1+l2+l3)/ltot && r1<=1)			{bind2=0;x1[n1]=s1-T0;n1++;s2=0;}
					}
				if(bind1==0 && bind2==0)			s0=0;
				else								s0=(k1*s1+k2*s2)/(bind1*k1+bind2*k2);
				if(bind1==0 && bind2==0)			{num++;s1=0;s2=0;s0=0;bind1=1;bind2=1;}
				if(num>=snap)						break;
				if(n1>99990 || n2>99990 || n3>99990|| n4>99990|| n5>99990|| n6>99990) break;

				}

			p1=n1/(n1+n2+0.0);	p2=1-p1;
			p3=n3/(n3+n4+0.0);	p4=1-p3;
			p5=n5/(n5+n6+0.0);	p6=1-p5;
			for(i=0,stot1=0;i<=n1;i++)											stot1+=x1[i];S1=l0*stot1/n1;
			for(i=0,stot1=0;i<=n2;i++)											stot1+=x2[i];S2=l0*stot1/n2;
			for(i=0,stot1=0;i<=n3;i++)											stot1+=x3[i];S3=l0*stot1/n3;
			for(i=0,stot1=0;i<=n4;i++)											stot1+=x4[i];S4=l0*stot1/n4;
			for(i=0,stot1=0;i<=n5;i++)											stot1+=x5[i];S5=l0*stot1/n5;
			for(i=0,stot1=0;i<=n6;i++)											stot1+=x6[i];S6=l0*stot1/n6;

			stot1=0;stot2=0;
			for(i=0;i<num;i++){
				stot1+=smax1[i];stot2+=smax2[i];}
			if(print_parameters)
				printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",v10,v20,e10,e20,Fs1,Fs2,Fd1,Fd2,k1,k2,pi1,pi2);
			}

		mid=0;
		for(i=0;i<imax;i++)
			for(j=0;j<jmax;j++){
				x=1;y=1;xy=1;
				for(K1=1;K1<=i;K1++)					x=x*K1;
				for(K1=1;K1<=j;K1++)					y=y*K1;
				for(K1=1;K1<=i+j;K1++)				xy=xy*K1;
				if(i==0)											x=1;
				if(j==0)											y=1;
				if(i==0 && j==0)							xy=1;
				mid+=xy/x/y*pow(p1*p4,i)*pow(p2*p5,j)*((p1*p3+p2*p6)*(i*(S1+S4)+j*(S2+S5))+p1*p3*(S1+S3)+p2*p6*(S2+S6));
				}
		printf("%lf\t%lf\n",mid,l0*stot1/num);
		//fprintf(k_d1_d2,"%lf\t%lf\t%lf\n",change_num,(1+(e10+pi1)/e20)/(1+(p2/p1/e10*e20)*(e10+pi2)/(pi1+e20)),(1+(e10+pi1)/e20)/(1+(e10+pi2)/(e10+e20)*(S1/v10*(e10+e20))));

		fprintf(k_d1_d2,"%lf\t%lf\n",change_num,l0*stot1/num/v10*e10);
		}

}
