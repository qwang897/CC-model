#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

#define Max_step_of_gillespie 			200000
#define Max_cycle_of_gillespie 			10000

#define Np_ini							6
#define Nm_ini							6
#define Np_max							6
#define Nm_max	 						6

#define choose_bidirection 				0
	#define bid_sum						0

#define printf_parameters				1
#define printf_delta					1
	#define printf_error				0
#define times 							5

int main(){
    double Fsp,Fsm,Fdp,Fdm;
    int i,j,k,I,m=0,Np,Nm,np,nm,num=0,s_C;
    int k_p=0,k_m=0;
    double Fc,Vc,gammam,gammap;
    double x = 0.0;
    double r1,r2,t = 0.0,T=0,tau,lam,l1,l2,l3,l4,l0;
    double gamma0p,gamma0m,pip,pim,Vfp,Vfm,Vbm,Vbp,s_max_e[100][100]={0};
    double V[201]={0},s_max[10000]={0},s_max_sum=0,s_max_average[100][100]={0};
    double run[1001]={0},runlength=0,direction=0,Direction[100]={0},runtau=0,Arv=0;
    double pause_t=0,pause_x=0,s_max_p[10000]={0},s_max_m[10000]={0};
    double s_max_average_m[100][100]={0},s_max_average_p[100][100]={0};
    double s_max_e_m[100][100]={0},s_max_e_p[100][100]={0};
    double x_p[10000]={0},x_m[10000]={0},x01=0,xp[100][100]={0},xm[100][100]={0};
    double sum1=0,sum2=0;
    double x0 = 0.0,T0,arv[10000],atau[10000]={0},runl[10000]={0};
    double x_p_sum[10000]={0},x_m_sum[10000]={0};
    double t_sum[10000]={0},T_sum[100][100]={0};
    double r21,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,mid1,mid2,mid3,mid4;
    k=0;
    srand((int)time(NULL));
    for(i= 1;i< 5;i++ )     Direction[i] = 0;
    for(i= 0;i< 101;i++ )   run[i] = 0;
    for(i = 0;i < 10000;i++ ){
    	arv[i] = 0;runl[i] = 0;}

	Fsp = 6;				Fsm = 1.1;
	Fdp = 3;    			Fdm = 0.75;
	pip = 5;    			pim = 1.6;
	Vfp = 1000;    			Vfm = 650;
	Vbp = 67;   			Vbm = 72;
	gamma0p = 1;			gamma0m = 0.27;

	if(printf_parameters==1)
		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",Vfp,Vfm,gamma0p,gamma0m,Fsp,Fsm,Fdp,Fdm,Vbp,Vbm,pip,pim);



	for(Np=Np_ini;Np<=Np_max;Np++)
		for(Nm=Nm_ini;Nm<=Nm_max;Nm++){
		for(i=1;i<10000;i++)					
			{s_max[i]=0;s_max_p[i]=0;s_max_m[i]=0;t_sum[i]=0;
			x_p[i]=0;x_m[i]=0;x_p_sum[i]=0;x_m_sum[i]=0;}


		nm = 0;np = 0;t=0;num=0;x=0;i=0;k_p=0;k_m=0;
		for(j = 0;j < Max_step_of_gillespie;j++){
			s_max[num]=x;
			r1 = rand()/(double)(RAND_MAX);
			r2 = rand()/(double)(RAND_MAX);
			if(r1 == 0)    r1 = r1 + 0.000001;
			// printf("%lf\t%lf\n",r1,r2);

			if(nm == 0)                    Fc = 0;
			else{
				if((np*Fsp) >= (nm*Fsm))
					{lam = 1/(1 + (np*Fsp*Vbm)/(nm*Fsm*Vfp));
					Fc = lam*np*Fsp + (1-lam)*nm*Fsm;}
				else
					{lam = 1/(1 + (np*Fsp*Vfm)/(nm*Fsm*Vbp));
					Fc = lam*np*Fsp + (1-lam)*nm*Fsm;}

			}
			if(np != 0)                  gammap = gamma0p*exp(Fc/(np*Fdp));
			else                         gammap = 0;
			if(nm != 0)                  gammam = gamma0m*exp(Fc/(nm*Fdm));
			else                         gammam = 0;

			l1 = gammap*np;                               //plus motor unbind
			l2 = pip*(Np-np);                             //plus motor bind
			l3 = gammam*nm;                               //minus motor unbind
			l4 = pim*(Nm-nm);                             //minus motor bind
			l0 = l1 + l2 + l3 + l4;

			if((r2 >= 0.0) && (r2 < l1/l0))
				np--;
			if((r2 >= l1/l0) && (r2 < (l1+l2)/l0))
				np++;
			if((r2 >= (l1+l2)/l0) && (r2 < (l1+l2+l3)/l0))
				nm--;
			if((r2 >= (l1+l2+l3)/l0) && (r2 < (l1+l2+l3+l4)/l0))
				nm++;

			if(nm == 0&&np == 0)            Vc = 0;
			else{
				if((np*Fsp) >= (nm*Fsm))    Vc = (np*Fsp-nm*Fsm)/(np*Fsp/Vfp+ nm*Fsm/Vbm);
				else                      	Vc = (np*Fsp-nm*Fsm)/(np*Fsp/Vbm+ nm*Fsm/Vfm);}
	//again
			if(nm == 0)                    	Fc = 0;
			else{
				if((np*Fsp) >= (nm*Fsm))
					{lam = 1/(1 + (np*Fsp*Vbm)/(nm*Fsm*Vfp));
					Fc = lam*np*Fsp + (1-lam)*nm*Fsm;}
				else
					{lam = 1/(1+ (np*Fsp*Vfm)/(nm*Fsm*Vbp));
					Fc = lam*np*Fsp+ (1-lam)*nm*Fsm;}}
			if(np!= 0)                    	gammap = gamma0p*exp(Fc/(np*Fdp));
			else                        	gammap = 0;
			if(nm!= 0)                  	gammam = gamma0m*exp(Fc/(nm*Fdm));
			else                        	gammam = 0;
			l1 = gammap*np;l2 = pip*(Np-np);l3 = gammam*nm;l4 = pim*(Nm-nm);
			l0 = l1 + l2 + l3 + l4;

			tau = (1.0/l0)*log(1.0/r1);
			x = x + tau*Vc;
			t  = t + tau;
			t_sum[num]+=tau;


		if(np==0 && nm==0){
			if(bid_sum==1){
					for(i=0;i<10000;i++)	{x_p_sum[num]+=x_p[i];x_m_sum[num]+=x_m[i];}
					for(i=0;i<10000;i++)	{x_p[i]=0;x_m[i]=0;}
					direction=0;k_p=0;k_m=0;}
					num++;
					if(choose_bidirection!=1 || bid_sum==1)	x=0;
				}



				if(choose_bidirection==1){
					if(nm == 0&&np == 0)          		{Vc = 0;
						if(bid_sum==1)					x01=0;}
					else{
						if((np*Fsp) >= (nm*Fsm))		Vc = (np*Fsp-nm*Fsm)/(np*Fsp/Vfp+ nm*Fsm/Vbm);
						else                      		Vc = (np*Fsp-nm*Fsm)/(np*Fsp/Vbm+ nm*Fsm/Vfm);}
						if(Vc>0){
							if(direction==-1){			k_p++;x01=x;}direction=1;}
						if(Vc<0){
							if(direction==1){			k_m++;x01=x;}direction=-1;}
						if(direction==1)				x_p[k_p]=x-x01;
						if(direction==-1)				x_m[k_m]=x-x01;}

			if(num>Max_cycle_of_gillespie)				break;
		}

	if(choose_bidirection==1 && bid_sum==1)
		for(k=1;k<=Max_cycle_of_gillespie;k++)		{x_p[k]=x_p_sum[k];x_m[k]=x_m_sum[k];}


	s_max_sum=0;
	for(k=1;k<Max_cycle_of_gillespie;k++)		s_max_sum+=s_max_p[k];
	s_max_average_p[Np][Nm]=s_max_sum/Max_cycle_of_gillespie;
	s_max_sum=0;
	for(k=1;k<Max_cycle_of_gillespie;k++)		s_max_sum+=pow(s_max_p[k]-s_max_average_p[Np][Nm],2);
	s_max_e_p[Np][Nm]=sqrt(s_max_sum/(Max_cycle_of_gillespie*(Max_cycle_of_gillespie-1)));
	s_max_sum=0;
	for(k=1;k<Max_cycle_of_gillespie;k++)		s_max_sum+=s_max_m[k];
	s_max_average_m[Np][Nm]=s_max_sum/Max_cycle_of_gillespie;
	s_max_sum=0;
	for(k=1;k<Max_cycle_of_gillespie;k++)		s_max_sum+=pow(s_max_m[k]-s_max_average_m[Np][Nm],2);
	s_max_e_m[Np][Nm]=sqrt(s_max_sum/(Max_cycle_of_gillespie*(Max_cycle_of_gillespie-1)));

	s_max_sum=0;mid1=0;
	for(k=1;k<Max_cycle_of_gillespie;k++)
		if(s_max[k]>=0)			{s_max_sum+=s_max[k];mid1++;}
	s_max_average_p[Np][Nm]=s_max_sum/mid1;
	s_max_sum=0;mid1=0;
	for(k=1;k<Max_cycle_of_gillespie;k++)
		if(s_max[k]<=0)			{s_max_sum+=s_max[k];mid1++;}
	s_max_average_m[Np][Nm]=s_max_sum/mid1;

	sum1=0;sum2=0;
	for(k=1;k<=Max_cycle_of_gillespie;k++){		sum1+=x_p[k];sum2+=x_m[k];}
	xp[Np][Nm]=sum1/Max_cycle_of_gillespie;xm[Np][Nm]=sum2/Max_cycle_of_gillespie;
	sum1=0;
	for(k=1;k<=Max_cycle_of_gillespie;k++)		sum1+=t_sum[k];
	T_sum[Np][Nm]=sum1/Max_cycle_of_gillespie;


}			//Np,Nm


	s_max_average[0][0]=0;s_max_e[0][0]=0;


	if(printf_delta==1){
		if(printf_error==1){
		for(Np=Np_ini;Np<=Np_max;Np++){
				for(Nm=Nm_ini;Nm<=Nm_max;Nm++){
						printf("%.2f_+_%.2f  \t",s_max_average_p[Np][Nm],s_max_e_p[Np][Nm]);}
						printf("\n");}
						printf("\n");
			for(Np=Np_ini;Np<=Np_max;Np++){
				for(Nm=Nm_ini;Nm<=Nm_max;Nm++){
						printf("%.2f_+_%.2f  \t",s_max_average_m[Np][Nm],s_max_e_m[Np][Nm]);}
						printf("\n");}}
						printf("\ndistance max +++\n");
		for(Np=Np_ini;Np<=Np_max;Np++){
			for(Nm=Nm_ini;Nm<=Nm_max;Nm++){
					printf("%.2f  \t",s_max_average_p[Np][Nm]);}
					printf("\n");}
					printf("\ndistance max ---\n");
		for(Np=Np_ini;Np<=Np_max;Np++){
			for(Nm=Nm_ini;Nm<=Nm_max;Nm++){
					printf("%.2f  \t",s_max_average_m[Np][Nm]);}
					printf("\n");}}
	if(choose_bidirection	==1){
					printf("\ndistance +++\n");
		for(Np=Np_ini;Np<=Np_max;Np++){
			for(Nm=Nm_ini;Nm<=Nm_max;Nm++){
					printf("%.2f  \t",xp[Np][Nm]);}
					printf("\n");}
					printf("\ndistance ---\n");
		for(Np=Np_ini;Np<=Np_max;Np++){
			for(Nm=Nm_ini;Nm<=Nm_max;Nm++){
					printf("%.2f  \t",xm[Np][Nm]);}
					printf("\n");}printf("\n");}



}



