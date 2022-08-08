#include <string>
#include <cmath>
#include <stdlib.h>
#include <ctime>
#include <cstdio>
#include <iostream>


#define Average_runlength_velocity		0		//Fig2------Should set Restart_from_the_origin to 0
#define Cooperative						0		//Fig3.A and Fig3.B
#define Balance_probability				1		//Fig3.C
#define Motors_number					0		//Fig3.D----Should set Np_ini=Np_max=2
#define Max_force						0
#define Show_progress_bar				0

#define Max_time_of_gillespie			100000
#define Max_cycle_of_gillespie			100
#define Max_repeat_of_gillespie			1
#define Parameter_is_random				0
#define Restart_from_the_origin			1		//default 1
#define Motors_can_move_backwards		0		//default 0
#define different_rate					0

#define Np_ini							2
#define Nm_ini							0
#define Np_max							2
#define Nm_max	 						5

#define time1							5

int main(){
	void parameter(double*);
	void gillespie(int,double*);
	void print_average_runlength_velocity(FILE*,FILE*);
	double input[30]={0};
	FILE* average=fopen("./average.txt","w");
	FILE* Output=fopen("./parameter_output.txt","w");
	srand((int)time(NULL));
	for(int i=0;i<Max_repeat_of_gillespie;i++){
		parameter(input);
		gillespie(Max_cycle_of_gillespie,input);
		if(Average_runlength_velocity==1)
			print_average_runlength_velocity(average,Output);
	}
	fclose(average);
}

void parameter(double input[30]){
	void specific(double*);
	void random(double*);
	switch(Parameter_is_random){
		case 0:
			specific(input);
		break;
		case 1:
			random(input);
		break;
	}
}
	void specific(double input[30]){
		double pi1=0,v10=0,e10=0,Fs1,Fs2,Fd1,Fd2,v1max,v2max,v1min,v2min;
		double pi2=0,v20=0,e20=0,k1,k2,v1back=0,v2back=0;
		if(Average_runlength_velocity==1){
			v10=380;        v20=806;
			e10=0.32;       e20=1.09;
			Fs1=2.50;       Fs2=2.50;
			Fd1=2.00;       Fd2=1.74;
			k1 =0.33;       k2=0.12;
			pi1=0.49;       pi2=0.59;
		}
		else{
			Fs1 = 6; 		Fs2 = 1.1;
			Fd1 = 3;    	Fd2 = 0.75;
			pi1 = 5;     	pi2 = 1.6;
			v10 = 1000;    	v20 = 650;
			v1back = 6;   	v2back = 72;
			e10 = 1;		e20 = 0.27;
			k1 = 0.2;		k2 = 0.05;
		}
		input[1]=v10;input[2]=v20;input[3]=e10;input[4]=e20;input[5]=Fs1;input[6]=Fs2;input[7]=Fd1;
		input[8]=Fd2;input[9]=k1;input[10]=k2;input[11]=pi1;input[12]=pi2;input[13]=v1back;input[14]=v2back;
		input[15]=v1max;input[16]=v2max;input[17]=v1min;input[18]=v2min;input[19]=Np_ini;input[20]=Nm_ini;
		input[21]=Np_max;input[22]=Nm_max;
	}
	void random(double input[30]){
		double times=time1;
		double r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19;
		double pi1=0,v10=0,e10=0,Fs1,Fs2,Fd1,Fd2;
		double pi2=0,v20=0,e20=0,k1,k2,v1back,v2back;
		r2 =(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
		r3 =(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
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
		r16=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
		r17=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
		r18=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;
		r19=(rand()/(double)(RAND_MAX)-0.5)*(times-1)*2;

		v10=463;        v20=632;
		e10=0.45;       e20=0.73;
		Fs1=2.50;       Fs2=2.50;
		Fd1=2.00;       Fd2=1.74;
		k1 =0.26;       k2=0.12;
		pi1=0.35;       pi2=0.40;
		v1back=8;       v2back=13;


		if(r2>=0)		v10=v10*(r2+1);
		else 			v10=v10/(-r2+1);
		if(r3>=0)		v20=v20*(r3+1);
		else 			v20=v20/(-r3+1);
		if(r4>=0)		e10=e10*(r4+1);
		else 			e10=e10/(-r4+1);
		if(r5>=0)		e20=e20*(r5+1);
		else 			e20=e20/(-r5+1);
		if(r6>=0)		Fs1=Fs1*(r6+1);
		else 			Fs1=Fs1/(-r6+1);
		if(r7>=0)		Fs2=Fs2*(r7+1);
		else 			Fs2=Fs2/(-r7+1);
		if(r8>=0)		Fd1=Fd1*(r8+1);
		else 			Fd1=Fd1/(-r8+1);
		if(r9>=0)		Fd2=Fd2*(r9+1);
		else 			Fd2=Fd2/(-r9+1);
		if(r10>=0)		k1 =k1*(r10+1);
		else 			k1 =k1/(-r10+1);
		if(r11>=0)		k2 =k2*(r11+1);
		else 			k2 =k2/(-r11+1);
		if(r12>=0)		pi1=pi1*(r12+1);
		else 			pi1=pi1/(-r12+1);
		if(r13>=0)		pi2=pi2*(r13+1);
		else 			pi2=pi2/(-r13+1);
		if(r14>=0)		v1back=v1back*(r14+1);
		else 			v1back=v1back/(-r14+1);
		if(r15>=0)		v2back=v2back*(r15+1);
		else		 	v2back=v2back/(-r15+1);
		input[1]=v10;input[2]=v20;input[3]=e10;input[4]=e20;input[5]=Fs1;input[6]=Fs2;input[7]=Fd1;
		input[8]=Fd2;input[9]=k1;input[10]=k2;input[11]=pi1;input[12]=pi2;input[13]=v1back;input[14]=v2back;
		input[19]=Np_ini;input[20]=Nm_ini;input[21]=Np_max;input[22]=Nm_max;
	}

    
void gillespie(int kmax,double input[30]){
	double  calculate_s0(int,int,int*,double,int*,double*,double*);
	void 	calculate_abeF(int,int,int*,double,double*,double*,double*,double*,double*,double*);
	void 	print(FILE*,double,double);
	void 	calculate_cooperative(int,int,int,double*,double(*)[100],double(*)[100]);
	void 	print_cooperative(double(*)[100],double(*)[100],double*);
	void 	calculate_balance_probability(int,int,int,int,double,double,double*,double*,double*,double*);
	void 	print_balance_probability(double*);
	void 	calculate_max_force(int,int,int,double*,double*,double*);
	void 	print_max_force(int,double*,double*);
	void 	calculate_motors_number(int,int,double,double(*)[100],double*,double*);
	void 	print_motors_number(double(*)[100],double*,double*);
	void 	print_final_length(int,double*,double*);
	int 	print_progress_bar(int,double);
	int 	i,j=0,k=0;
	int 	Np,Nm,np=0,nm=0,motor[2]={0},unbind=0,bind[100];
	double 	r1,r2,tau,T=0,v_motor=0,n_tau_av_aF[10]={0};
	double 	pi0[100],F[100],aF[100],bF[100],eF[100];
	double 	s0=0,s1,s[100],s_fin[100000];
	double 	ltot,lnow,l[300];
	double 	s_average[100][100],s_e[100][100];
	double 	kinesin_force[kmax],dynein_force[kmax];
	double 	numkd[100][100]={0},numtau[2]={0};

	FILE* fp=fopen("./output.txt","w");
	for(int i=1;i<=22;i++)	fprintf(fp,"%lf\t",input[i]);fprintf(fp,"\n");
	for(Np=Np_ini;Np<=Np_max;Np++)
	for(Nm=Nm_ini;Nm<=Nm_max;Nm++){
		for(k=0;k<100000;k++)	s_fin[k]=0;
		k=0;
		for(i=1;i<=Np+Nm;i++){
			if(i<=Np){			pi0[i]=input[11];}
			else if(i>Np){	pi0[i]=input[12];}}
		srand((int)time(NULL));
		while(1){
				r1=rand()/(double)(RAND_MAX);
				r2=rand()/(double)(RAND_MAX);
				if(r1==0)    r1=r1+0.000001;
				s0=calculate_s0(Np,Nm,motor,s0,bind,input,s);
				calculate_abeF(Np,Nm,bind,s0,s,F,input,aF,bF,eF);
				ltot=0;lnow=0;s1=s0;unbind=0;np=motor[0];nm=motor[1];
				for(i=1;i<=Np+Nm;i++)
						if(bind[i]==0)				unbind++;
				if(unbind==Np+Nm)					k++;
				s_fin[k]=s0;
				print(fp,T,s0);
				for(i=1;i<=Np+Nm;i++){
					l[3*i-2]=aF[i];l[3*i-1]=bF[i];l[3*i-0]=eF[i];
					if(std::isinf(l[3*i-2]))		l[3*i-2]=pow(10,100);
					if(std::isinf(l[3*i-1]))		l[3*i-1]=pow(10,100);
					if(std::isinf(l[3*i-0]))		l[3*i-0]=pow(10,100);
					if(bind[i]==0)					{l[3*i-2]=pi0[i];l[3*i-1]=0;l[3*i]=0;}
					ltot+=l[3*i]+l[3*i-1]+l[3*i-2];
					}
			for(i=1;i<=Nm+Np;i++){
					if(r2>=lnow/ltot && r2<(lnow+l[3*i-2])/ltot){
						if(i<=Np){
							if(bind[i]==0)			s[i]=round(s0);
							else					s[i]++;}
						else{
							if(bind[i]==0)			s[i]=round(s0);
							else					s[i]--;}
						bind[i]=1;
						break;}
					else if(r2>=(lnow+l[3*i-2])/ltot && r2<=(lnow+l[3*i-2]+l[3*i-1])/ltot){
						if(i<=Np)					s[i]--;
						else						s[i]++;
						bind[i]=1;
						break;}
					else if(r2>=(lnow+l[3*i-2]+l[3*i-1])/ltot && r2<=(lnow+l[3*i-2]+l[3*i-1]+l[3*i])/ltot){
						s[i]=0;bind[i]=0;
						break;}
					lnow+=l[3*i-2]+l[3*i-1]+l[3*i];}

				tau=(1.0/ltot)*log(1.0/r1);T+=tau;
				s0=calculate_s0(Np,Nm,motor,s0,bind,input,s);
				calculate_abeF(Np,Nm,bind,s0,s,F,input,aF,bF,eF);
				v_motor=(s0-s1)*8/tau;
				if(Max_force==1)					calculate_max_force(Np,Nm,k,F,kinesin_force,dynein_force);
				if(Balance_probability==1)			calculate_balance_probability(Np,Nm,np,nm,tau,v_motor,n_tau_av_aF,F,aF,input);
				if(Motors_number==1)				calculate_motors_number(np,nm,tau,numkd,numtau,input);
				if(Show_progress_bar==1)			j=print_progress_bar(j,T);
				if(Average_runlength_velocity==1 && T>Max_time_of_gillespie)	break;
				if(k>kmax)  break;
			}
		if(Balance_probability==1)					print_balance_probability(n_tau_av_aF);
		if(Max_force==1)							print_max_force(kmax,kinesin_force,dynein_force);
		if(Cooperative==1)							calculate_cooperative(Np,Nm,kmax,s_fin,s_average,s_e);
	}
	if(Cooperative==1)								print_cooperative(s_average,s_e,input);
	if(Motors_number==1)							print_motors_number(numkd,numtau,input);
	fclose(fp);
}
	double calculate_s0(int Np,int Nm,int motor[2],double s0,int bind[100],double input[30],double s[100]){
		int i,np,nm;
		double Sp_mid1,Sm_mid1;
		np=0;nm=0;Sp_mid1=0;Sm_mid1=0;
		for(i=1;i<=Np;i++){
			if(bind[i]!=0){										np++;Sp_mid1+=s[i];}}
		for(i=Np+1;i<=Np+Nm;i++){
			if(bind[i]!=0){										nm++;Sm_mid1+=s[i];}}
		if(np==0 && nm==0){
			if(Restart_from_the_origin==1)		s0=0;}
		else																s0=(input[9]*Sp_mid1+input[10]*Sm_mid1)/(np*input[9]+nm*input[10]);
		motor[0]=np;motor[1]=nm;
		return s0;
	}
	void 	calculate_abeF(int Np,int Nm,int bind[100],double s0,double s[100],double F[100],double input[30],double aF[100],double bF[100],double eF[100]){
		int i;
		double mid1,mid2[100];
		double Fs[100],Fd[100],v0[100],e0[100],VB[100],VF[100],vF[100],qF[100],vmin[100],vmax[100],l0=8,q0=800;
		for(i=1;i<=Np+Nm;i++){
				if(i<=Np)							F[i]=input[9]*l0*(s[i]-s0);
				else								F[i]=input[10]*l0*(0-s[i]+s0);
				if(bind[i]==0)						F[i]=0;
				if(F[i]==Fs[i])						F[i]+=0.000001;}
		for(i=1;i<=Np+Nm;i++){
			if(i<=Np){		Fs[i]=input[5];Fd[i]=input[7];e0[i]=input[3];VB[i]=input[13];VF[i]=input[1];v0[i]=input[1];vmin[i]=input[17];vmax[i]=input[15];}
			else if(i>Np){	Fs[i]=input[6];Fd[i]=input[8];e0[i]=input[4];VB[i]=input[14];VF[i]=input[2];v0[i]=input[2];vmin[i]=input[18];vmax[i]=input[16];}}

		switch(1){
		case 1:
			for(i=1;i<=Np+Nm;i++){
				if(different_rate==1){
					if(i<=Np){
						if(F[i]<Fs[i]){
							vF[i]=VF[i]*(1-pow(F[i]/Fs[i],1));
							eF[i]=e0[i]*exp(F[i]/Fd[i]);}
						else{
							vF[i]=0;
							eF[i]=e0[i]*(1.535+0.186*F[i]);}}
					else{
						if(F[i]<Fs[i]){
							vF[i]=VF[i]*(1-pow(F[i]/Fs[i],1));
							eF[i]=e0[i]*exp(F[i]/Fd[i]);}
						else{
							vF[i]=0;
							eF[i]=e0[i]*(1/(0.254*(1-exp(-F[i]/1.96646))));}}}
				else{
					if(F[i]<=Fs[i])					vF[i]=VF[i]*(1-F[i]/Fs[i]);
					else							vF[i]=-VB[i]*(1-F[i]/Fs[i]);
					eF[i]=e0[i]*exp(fabs(F[i])/Fd[i]);}
				if(Motors_can_move_backwards==0){
					if(F[i]<Fs[i])					{bF[i]=0;aF[i]=vF[i]/l0;}
					if(F[i]>=Fs[i])					{bF[i]=0;aF[i]=0;}}
				else if(Motors_can_move_backwards==1){
					if(F[i]<Fs[i])					{bF[i]=0;aF[i]=vF[i]/l0;}
					if(F[i]>=Fs[i])					{bF[i]=vF[i]/l0;aF[i]=0;}}
				if(std::isinf(eF[i]))				eF[i]=pow(10,100);}
		break;
		case 2:
			for(i=1;i<=Np+Nm;i++){
				qF[i]=pow(q0,1-F[i]/Fs[i]);
				mid1=(vmin[i]-v0[i])/(v0[i]-vmax[i]);
				mid2[i]=pow((vmax[i]/vmin[i])*(-mid1),F[i]/Fs[i]);
				vF[i]=(vmax[i]*mid1+vmin[i]*mid2[i])/(mid1+mid2[i]);
				aF[i]=(vF[i]/l0)*(qF[i]/(qF[i]-1));
				bF[i]=(vF[i]/l0)*(1/(qF[i]-1));
				if(Motors_can_move_backwards==0)	bF[i]=0;
				eF[i]=e0[i]*exp(fabs(F[i])/Fd[i]);
				if(std::isinf(eF[i]))				eF[i]=pow(10,100);
			}
		break;
		}
	}
	void 	print(FILE* fp,double T,double s0){
	switch(1){
		case 1:	fprintf(fp,"%lf\t%lf\n",T,s0+0.0);	break;
		case 2:	fprintf(fp,"%lf\n",T);				break;
	}
}
	void calculate_cooperative(int Np,int Nm,int kmax,double s_fin[],double s_average[100][100],double s_e[100][100]){
	int k=0;
	double s_sum=0;
	s_sum=0;
	for(k=1;k<kmax;k++)		s_sum+=s_fin[k];
	s_average[Np][Nm]=s_sum/kmax;
	s_sum=0;
	for(k=1;k<kmax;k++)		s_sum+=pow(s_fin[k]-s_average[Np][Nm],2);
	s_e[Np][Nm]=sqrt(s_sum/(kmax*(kmax-1)));
}
	void print_cooperative(double s_average[100][100],double s_e[100][100],double input[30]){
	int Np=0,Nm=0;
	for(Np=Np_ini;Np<=Np_max;Np++){
		for(Nm=Nm_ini;Nm<=Nm_max;Nm++){
			printf("%.2f±%.2f  \n",s_average[Np][Nm]*8,s_e[Np][Nm]*8);}
			printf("\n");}
}
	void calculate_balance_probability(int Np,int Nm,int np,int nm,double tau,double v_motor,double n_tau_av_aF[10],double F[100],double aF[100],double input[30]){
	int i;
	double l0=8,v[100]={0},average_vp=0,average_Fp=0,cop_v_sum=0;
	double average_vm=0,average_Fm=0,average_v_sum=0;
	double average_vp_sum=0,average_vm_sum=0,average_Fm_sum=0,average_Fp_sum=0;
	if(nm!=Nm || np!=Np)							return;
	for(i=1;i<=Np;i++){
		if(F[i]<input[5])							v[i]=aF[i]*l0;
		if(F[i]>=input[5])							v[i]=-aF[i]*l0;
		average_vp+=v[i];average_Fp+=F[i];
		if(fabs((v[i]-v_motor)/v_motor)<0.2)		cop_v_sum+=1.0/(np+nm)*tau;
		}
	for(i=Np+1;i<=Nm+Np;i++){
		if(F[i]<input[6])							v[i]=-aF[i]*l0;
		if(F[i]>=input[6])							v[i]=aF[i]*l0;
		average_vm+=v[i];average_Fm+=F[i];
		if(fabs((v[i]-v_motor)/v_motor)<0.2)		cop_v_sum+=1/(np+nm)*tau;
		}
	average_vp_sum+=average_vp/np*tau;average_vm_sum+=average_vm/nm*tau;
	average_Fp_sum+=average_Fp/np*tau;average_Fm_sum+=average_Fm/nm*tau;
	average_v_sum+=average_vp/np-average_vm/nm;
	n_tau_av_aF[1]++;n_tau_av_aF[2]+=tau;
	n_tau_av_aF[3]+=cop_v_sum;
	n_tau_av_aF[4]+=average_Fp/np*tau;
}
	void print_balance_probability(double n_tau_av_aF[10]){
		printf("Velocity balance probability\t%.3f\n",n_tau_av_aF[3]/n_tau_av_aF[2]);
}
	void calculate_max_force(int Np,int Nm,int k,double F[100],double kinesin_force[],double dynein_force[]){
		int i;
    double kinesin_force_tot=0,dynein_force_tot=0;
    for(i=1;i<=Np;i++)								kinesin_force_tot+=F[i];
    for(i=Np+1;i<=Np+Nm;i++)						dynein_force_tot +=F[i];
    if(kinesin_force[k]<=kinesin_force_tot)			kinesin_force[k]=kinesin_force_tot;
    if(dynein_force[k] <=dynein_force_tot )			dynein_force[k] =dynein_force_tot ;
}
	void print_max_force(int kmax,double kinesin_force[],double dynein_force[]){
	int k;
	double sum1=0,sum2=0;
	for(k=1;k<=kmax;k++){		sum1+=kinesin_force[k];sum2+=dynein_force[k];}
	sum1=sum1/kmax;sum2=sum2/kmax;
	printf("Max force\t%.3f\t%.3f\n",sum1,sum2);
}
	void calculate_motors_number(int np,int nm,double tau,double numkd[100][100],double numtau[2],double input[30]){
	int i,i1;
	for(i=0;i<=Np_max;i++)
			for(i1=0;i1<=Nm_max;i1++)
				if(i==np && i1==nm){
					numkd[i][i1]+=tau;
					numtau[1]+=tau;
			}
}
	void print_motors_number(double numkd[100][100],double numtau[2],double input[30]){
	int i,i1;
	for(i=0;i<=Np_max;i++){
			for(i1=0;i1<=Nm_max;i1++)
					printf("%.3f\t",numkd[i][i1]/numtau[1]);
			printf("\n");
			}
}


void print_average_runlength_velocity(FILE* average,FILE* Output){
	double 	pause_x=30,pause_t=0.23,reverse_t=0.16,v_min=50;
	int i,j,k,num=0,pause_start=0,n=0,plus=0,minus=0;
	double t[30000]={0},x[30000]={0},x0=0,t0=0;
	double txt[50000],run_plus[50000]={0},run_minus[50000]={0};
	double run_plus_sum=0,run_minus_sum=0,run_plus_e=0,run_minus_e=0;
	double run[300000]={0},tau[300000]={0};
	FILE *fp;
	fp=fopen("./output.txt","r");
	i=1;
	while(fscanf(fp, "%lf", &txt[i]) != EOF)
		{i++;
		if(i>22)		break;}
	double l1=txt[1],l2=txt[2],l3=txt[3],l4=txt[4],l5=txt[5],l6=txt[6],l7=txt[7],l8=txt[8],l9=txt[9];
	double l10=txt[10],l11=txt[11],l12=txt[12],l13=txt[13],l14=txt[14],l15=txt[19],l16=txt[20];
	for(k=0;k<20;k++){
		i=0;
		while(fscanf(fp, "%lf", &txt[i]) != EOF)
			{i++;
			if(i>=20000)		break;}
		for(i=0;i<20000;i++){
			t[i]=txt[2*i];x[i]=txt[2*i+1]*8;}
		for(i=1;i<10000;i++){
			if((x[i]-x[i-1])*(x[i+1]-x[i]) >= 0){
				run[num]+=x[i]-x[i-1];tau[num]+=t[i]-t[i-1];
				for(n=i;fabs(x0)<pause_x;n++){
					x0+=x[n]-x[n-1];t0+=t[n]-t[n-1];
					if(fabs(x0+x[n+1]-x[n])>pause_x){
						if(fabs(x0/t0)<=v_min && t0>pause_t){
							if(fabs(run[num]+x0)<pause_x && fabs((run[num]+x0)/(t[num]+t0))<v_min){
								run[num]+=x0;tau[num]+=t0;num=num+1;}
							else{
								run[num+1]=x0;tau[num+1]=t0;num=num+2;}
							i=n;
							break;
						}
						else if(t0>pause_t){
							for(j=n;t0>pause_t;j--){
								x0-=x[j]-x[j-1];t0-=t[j]-t[j-1];
								if(fabs(x0/t0)<=v_min){
									if(fabs(run[num]+x0)<pause_x && fabs((run[num]+x0)/(t[num]+t0))<v_min){
										run[num]+=x0;tau[num]+=t0;num=num+1;}
									else{
										run[num+1]=x0;tau[num+1]=t0;num=num+2;}
									i=j;
									break;
								}
							}
							if(i==j)	break;}}}
				x0=0;t0=0;
			}
			else{
				if(fabs(run[num])<=pause_x && tau[num]<reverse_t){
					run[num-1]+=run[num];tau[num-1]+=tau[num];
					run[num]=0;tau[num]=0;
				}
				else
					num++;
				x0=0;t0=0;
			}
			if(num>299990)	break;
		}
		if(num>299990)	break;
	}
	for(i=0;i<=num;i++)
		if(fabs(run[i])>30 && fabs(run[i+1])>30 && (run[i]*run[i+1])>=0){
			run[i+1]=run[i]+run[i+1];tau[i+1]=tau[i]+tau[i+1];run[i]=0;tau[i]=0;}
	for(i=0;i<=num;i++){
		if(fabs(run[i])<30 && !(tau[i]>0.23 && fabs(run[i]/tau[i])<50))	tau[i]=0;
		if(fabs(run[i])>30 && !(tau[i]>0.16 && fabs(run[i]/tau[i])>50))	tau[i]=0;
	}
	for(i=0,j=0;i<=num;i++)
		if(tau[i]>0)		{tau[j]=tau[i];run[j]=run[i];j++;}num=j;

	int 	pause_num=0;
	double 	pause_duration[50000]={0},duration_run=0,time_between_pause[50000]={0};
	double 	sum_pause_duration=0,sum_time_between_pause=0,error_pause_duration=0;
	double	quick_reversal=0,slow_reversal=0;

	for(i=0,pause_start=0;i<=num;i++){
		duration_run+=tau[i];
		if(fabs(run[i])<30)			{pause_duration[pause_num]=tau[i];pause_start=1;pause_num++;}
		else						pause_start=0;
		if(pause_start==0){
			time_between_pause[pause_num]+=tau[i];
			for(j=i-1;fabs(run[j])<30;j--)
				time_between_pause[pause_num]+=tau[j];
			}
		if((run[i]*run[i+1])<0){
			if(fabs(run[i])<30)		slow_reversal++;
			else					quick_reversal++;}}
	for(i=0;i<pause_num;i++){
		sum_pause_duration+=pause_duration[i];}
	for(i=0;i<pause_num;i++)
		error_pause_duration+=pow(pause_duration[i]-sum_pause_duration/pause_num,2);
	error_pause_duration=sqrt(error_pause_duration/(pause_num*(pause_num-1)));

	for(i=0;i<=num;i++){
		if(run[i]>pause_x	){		run_plus[plus]=run[i];plus++;}
		if(run[i]<-pause_x	){		run_minus[minus]=run[i];minus++;}
	}
	for(i=0;i<plus;i++)				run_plus_sum+=run_plus[i];
	for(i=0;i<minus;i++)			run_minus_sum+=run_minus[i];
	for(i=0;i<plus;i++)
		run_plus_e+=pow(run_plus[i]-run_plus_sum/plus,2);
	run_plus_e=sqrt(run_plus_e/(plus*(plus-1)));
	for(i=0;i<minus;i++)
		run_minus_e+=pow(run_minus[i]-run_minus_sum/minus,2);
	run_minus_e=sqrt(run_minus_e/(minus*(minus-1)));

	printf("Np=%.0f \tNm=%.0f\nv10=%.0f \tv20=%.0f\ne10=%.2f\te20=%.2f\n",l15,l16,l1,l2,l3,l4);
	printf("Fs1=%.2f\tFs2=%.2f\nFd1=%.2f\tFd2=%.2f\nk1 =%.2f\t",l5,l6,l7,l8,l9);
	printf("k2=%.2f\npi1=%.2f\tpi2=%.2f\nv1back=%.0f\tv2back=%.0f\n\n",l10,l11,l12,l13,l14);

	printf("/*\n%% duration paused  = %.0f%%\n",sum_pause_duration*100/duration_run);
	printf("Pause duration     = %.3f±%.3f\n",sum_pause_duration/pause_num,error_pause_duration);
	printf("runlength  plus=%.3f±%.3f\nrunlength minus=%.3f±%.3f\n*/\n",run_plus_sum/plus,run_plus_e,run_minus_sum/minus,run_minus_e);

	fprintf(average,"v10=%.0f \tv20=%.0f\ne10=%.2f\te20=%.2f\n",l1,l2,l3,l4);
	fprintf(average,"Fs1=%.2f\tFs2=%.2f\nFd1=%.2f\tFd2=%.2f\nk1 =%.2f\t",l5,l6,l7,l8,l9);
	fprintf(average,"k2=%.2f\npi1=%.2f\tpi2=%.2f\nv1back=%.0f\tv2back=%.0f\n\n",l10,l11,l12,l13,l14);
	fprintf(average,"runlength  plus=%.3f±%.3f\nrunlength minus=%.3f±%.3f\n\n\n",run_plus_sum/plus,run_plus_e,run_minus_sum/minus,run_minus_e);

}

int 	 print_progress_bar(int i,double T){
	if((T+1)/3000>i)	{printf("%d%%\n",10*i);i++;}
	if(i==30)	{printf("100%%\n");i++;}
	return i;
}


























