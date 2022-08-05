function can_not_back
syms x;
e10=1;          e20=0.27;
v10=1000;     v20=650;
Fs1=6;          Fs2=1.1;
Fd1=3;          Fd2=0.75;
pi1=5;          pi2=1.6;
l=8;
k1=0.2;         k2=0.05;
k=k1*k2/(k1+k2);
choose=1;
number=50;
curve=1;
if curve==1
	vmax=747.4;vmin=-88.8;
end
if curve==2
	vmax=1480;vmin=-7.4;
end

[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17]= ...,
    textread('./output.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');


N=length(x13);
sum_stot(1)=0;
stot0(1)=0;
for n=1:N
        v10=x1(n);
        v20=x2(n);
        e10=x3(n);
        e20=x4(n);
        Fs1=x5(n);
        Fs2=x6(n);
        Fd1=x7(n);
        Fd2=x8(n);
        k1=x9(n);
        k2=x10(n);
        pi1=x11(n);
        pi2=x12(n);  
        if choose==2
            v20=-x2(n);
            Fs2=-x6(n);
        end
        k=k1*k2/(k1+k2);
        
        if curve==0
            vF1=v10*(1-x/Fs1)*((abs(Fs1-x))/(Fs1-x)+1)/2;
            vF2=v20*(1-x/Fs2)*((abs(Fs2-x))/(Fs2-x)+1)/2;
            e1=e10.*exp(abs(x)./Fd1);
            e2=e20.*exp(abs(x)./Fd2);
            f=k*(vF1+vF2)/(e1+e2)-x;
            F=vpasolve(f);
            a1=v10./l.*(1-F/Fs1)*((abs(Fs1-F))/(Fs1-F)+1)/2;
            a2=v20./l.*(1-F/Fs2)*((abs(Fs2-F))/(Fs2-F)+1)/2;
            e1=e10.*exp(abs(F)./Fd1);
            e2=e20.*exp(abs(F)./Fd2);
            fs1=a1/v10;
            fs2=a2/v20;
            fd1=exp(abs(F)./Fd1);
            fd2=exp(abs(F)./Fd2);
        end
        
        if curve~=0
            e1=e10.*exp(abs(x)./Fd1);
            e2=e20.*exp(abs(x)./Fd2);
            mid11=(vmin-v10)/(v10-vmax);
            mid12=((vmax/vmin)*(-mid11))^(x/Fs1);
            vF1=(vmax*mid11+vmin*mid12)/(mid11+mid12);
            mid21=(vmin-v20)/(v20-vmax);
            mid22=((vmax/vmin)*(-mid21))^(x/Fs2);
            vF2=(vmax*mid21+vmin*mid22)/(mid21+mid22);
            f=k*(vF1+vF2)/(e1+e2)-x;
            F=vpasolve(f);
            e1=e10.*exp(abs(F)./Fd1);
            e2=e20.*exp(abs(F)./Fd2);
            mid11=(vmin-v10)/(v10-vmax);
            mid12=((vmax/vmin)*(-mid11))^(F/Fs1);
            a1=(vmax*mid11+vmin*mid12)/(mid11+mid12);
            mid21=(vmin-v20)/(v20-vmax);
            mid22=((vmax/vmin)*(-mid21))^(F/Fs2);
            a2=(vmax*mid21+vmin*mid22)/(mid21+mid22);
            fs1=a1/v10;
            fs2=a2/v20;
            fd1=exp(abs(F)./Fd1);
            fd2=exp(abs(F)./Fd2);
        end

        

        
        mid1=fd2.*(e20+pi1);
        mid2=fs1.*(e10+pi2)+e1+e2;
        mid3=fd1.*(e10+pi2);
        mid4=fs2.*(e20+pi1)+e1+e2;
        mid5=(e1+e2).*(e1+e2+fd2.*pi1+fd1.*pi2);
        
        s1=mid1.*mid2./mid5.*(v10./e10);
        s2=mid3.*mid4./mid5.*(v20./e20);
        s=s1-s2;
        sum_stot(n)=s;
        stot0(n)=x13(n);
        mid00=[stot0;sum_stot];
        mid00=mid00';
        save('matlab_output.txt','mid00','-ascii')
        


end
mid00
%     n=100:100:5000
%     plot(n,sum_s1);
%     hold on;
%     plot(n,sum_s2);hold on;
%     plot(n,sum_F);
%     sum_p1'


     x=-10000:10000:50000;
     y=-10000:10000:50000;
%     sum_stot
%     x16=vpa(x13',6)
%     x17=vpa(sum_stot,6)
%     mid00=[x13';sum_stot]
%     mid00=mid00'
%     save('output.txt','mid00','-ascii')
     plot(x13,sum_stot,'o',x,y)
     
     
     
