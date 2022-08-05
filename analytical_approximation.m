function accurate_k
syms x10;
global Fs1
global Fs2
global Fd1
global Fd2
global ktot
global e10
global e20
global v10
global v20
global l
global n
global sum
global fs1
global fs2
global f
[x1,x2,x3]=textread('./output.txt','%f %f %f');

e10=1;          e20=0.27;
v10=1000;     v20=650;
Fs1=6;          Fs2=1.1;
Fd1=3;          Fd2=0.75;
pi1=5;          pi2=1.6;
l=8;
k1=0.2;         k2=0.05;


lam(1)=0;
NN=10;
small_p=0;
for N=1:NN
    ktot=k1*k2/(k1+k2);
    l=8;k=2;
    sum=linspace(0,0,20);
    j=linspace(0,0,20);
ktot=100*N;
    j(1)=-1;
    Pk_sum=0;Pd_sum=0;Sk=0;Sd=0;
    for x=0:14
        for y=0:10
            for i=0:x+y
                Fi=i*ktot*l;
                fs1(i+1)=1-Fi/Fs1;
                fs2(i+1)=1-Fi/Fs2;
                if Fi>Fs1
                    fs1(i+1)=10^(-8);
                end
                if Fi>Fs2
                    fs2(i+1)=10^(-8);
                end
                f(i+1)=fs2(i+1)./fs1(i+1);
            end
            mid1=((v10/l)^x)*((v20/l)^y);
            mid2=tot1(x,y);
            n=0;
            for i=1:20
                sum(i)=0;j(i+1)=0;
            end
            sum(y+2)=1;mid3=0;
            if y==0
                mid3=1;
            end
                if y==1
                    mid3=lengthsum1(x,y);
                end
                if y==2
                    mid3=lengthsum2(x,y);
                end
                if y==3
                    mid3=lengthsum3(x,y);
                end
                if y==4
                    mid3=lengthsum4(x,y);
                end
                if y==5
                    mid3=lengthsum5(x,y);
                end
                if y==6
                    mid3=lengthsum6(x,y);
                end
                if y==7
                    mid3=lengthsum7(x,y);
                end
                if y==8
                    mid3=lengthsum8(x,y);
                end            
                if y==9
                    mid3=lengthsum9(x,y);
                end            
                if y==10
                    mid3=lengthsum10(x,y);
                end            
                if y==11
                    mid3=lengthsum11(x,y);
                end
                if y==12
                    mid3=lengthsum12(x,y);
                end
                if y==13
                    mid3=lengthsum13(x,y);
                end
                if y==14
                    mid3=lengthsum14(x,y);
                end
                if y==15
                    mid3=lengthsum15(x,y);
                end            
            res=mid1*mid2*mid3;
            Pkxy=res;
            F=(x+y)*ktot*l;
            Pdxy=res*e10/e20*exp((abs(F))*(1/Fd1-1/Fd2));
            if isnan(Pdxy)
                Pdxy=0;
            end
            Pk_sum=Pk_sum+Pkxy;
            Pd_sum=Pd_sum+Pdxy;
            Sk=Sk+Pkxy*x*l;
            Sd=Sd-Pdxy*y*l;
        end
    end

    Sk=Sk/Pk_sum;
    Sd=Sd/Pd_sum;
    p1=Pk_sum
    p2=Pd_sum
    s1=Sk;
    s2=Sd;
    p1=vpa(Pk_sum/(Pk_sum+Pd_sum),5)
    p2=vpa(1-p1,5)
    
    
    p10=e20/(e10+e20);
    p20=1-p10;
    s10=v10/(e10+e20);
    s20=-v20/(e10+e20);
    p3=e10/(e10+pi2);
    p4=1-p3;
    p5=pi1/(e20+pi1);
    p6=1-p5;
    s3=v10/(e10+pi2);
    s4=s3;
    s5=-v20/(e20+pi1);
    s6=s5;

    mid1=p1*p3*(1-p2*p5)*(s1+s3)+p1*p3*p2*p5*(s2+s5);
    mid2=p2*p6*(1-p1*p4)*(s2+s6)+p1*p2*p4*p6*(s1+s4);
    mid3=(1-p1*p4-p2*p5)^2;
    s=(mid1+mid2)/mid3;
    
    lam(N)=(1+(e10+pi2)/e20)/(1+p2/p1/p20*p10*(e10+pi2)/pi1+e20)


end
lam=lam'
N=1:NN;
plot(N,lam);

end 

function lengthsum1=lengthsum1(x,y)
global f
lengthsum1=0;
for j1=0:x+y-1

                    mid5=f(j1+1);
                    lengthsum1=lengthsum1+mid5;
           
end
end
function lengthsum2=lengthsum2(x,y)
global f
lengthsum2=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
       mid5=f(j1+1)*f(j2+1);
       lengthsum2=lengthsum2+mid5;
  end
end
end
function lengthsum3=lengthsum3(x,y)
global f
lengthsum3=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1

                    mid5=f(j1+1)*f(j2+1)*f(j3+1);
                    lengthsum3=lengthsum3+mid5;
    end
  end
end
end
function lengthsum4=lengthsum4(x,y)
global f
lengthsum4=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1

                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1);
                    lengthsum4=lengthsum4+mid5;
              end
    end
  end
end
end
function lengthsum5=lengthsum5(x,y)
global f
lengthsum5=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1
        for j5=j4+1:x+y-1
                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1)*f(j5+1);
                    lengthsum5=lengthsum5+mid5;
              end

      end
    end
  end
end
end
function lengthsum6=lengthsum6(x,y)
global f
lengthsum6=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1
        for j5=j4+1:x+y-1
            for j6=j5+1:x+y-1

                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1)*f(j5+1)*f(j6+1);
                    lengthsum6=lengthsum6+mid5;
              end
        end
      end
    end
  end
end
end
function lengthsum7=lengthsum7(x,y)
global f
lengthsum7=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1
        for j5=j4+1:x+y-1
            for j6=j5+1:x+y-1
                for j7=j6+1:x+y-1

                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1)*f(j5+1)*f(j6+1)*f(j7+1);
                    lengthsum7=lengthsum7+mid5;
                end
            end
        end
      end
    end
  end
end
end
function lengthsum8=lengthsum8(x,y)
global f
lengthsum8=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1
        for j5=j4+1:x+y-1
            for j6=j5+1:x+y-1
                for j7=j6+1:x+y-1
                    for j8=j7+1:x+y-1

                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1)*f(j5+1)*f(j6+1)*f(j7+1)*f(j8+1);
                    lengthsum8=lengthsum8+mid5;

                    end
                end
            end
        end
      end
    end
  end
end
end
function lengthsum9=lengthsum9(x,y)
global f
lengthsum9=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1
        for j5=j4+1:x+y-1
            for j6=j5+1:x+y-1
                for j7=j6+1:x+y-1
                    for j8=j7+1:x+y-1
                        for j9=j8+1:x+y-1
                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1)*f(j5+1)*f(j6+1)*f(j7+1)*f(j8+1) ...,
                        *f(j9+1);
                    lengthsum9=lengthsum9+mid5;
                            end
                    end
                end
            end
        end
      end
    end
  end
end
end
function lengthsum10=lengthsum10(x,y)
global f
lengthsum10=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1
        for j5=j4+1:x+y-1
            for j6=j5+1:x+y-1
                for j7=j6+1:x+y-1
                    for j8=j7+1:x+y-1
                        for j9=j8+1:x+y-1
                            for j10=j9+1:x+y-1
                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1)*f(j5+1)*f(j6+1)*f(j7+1)*f(j8+1) ...,
                        *f(j9+1)*f(j10+1);
                    lengthsum10=lengthsum10+mid5;
                            end
                        end
                    end
                end
            end
        end
      end
    end
  end
end
end
function lengthsum11=lengthsum11(x,y)
global f
lengthsum11=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1
        for j5=j4+1:x+y-1
            for j6=j5+1:x+y-1
                for j7=j6+1:x+y-1
                    for j8=j7+1:x+y-1
                        for j9=j8+1:x+y-1
                            for j10=j9+1:x+y-1
                                for j11=0:x+y-1
                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1)*f(j5+1)*f(j6+1)*f(j7+1)*f(j8+1) ...,
                        *f(j9+1)*f(j10+1);
                    lengthsum11=lengthsum11+mid5;
                                end
                            end
                        end
                    end
                end
            end
        end
      end
    end
  end
end
end
function lengthsum12=lengthsum12(x,y)
global f
lengthsum12=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1
        for j5=j4+1:x+y-1
            for j6=j5+1:x+y-1
                for j7=j6+1:x+y-1
                    for j8=j7+1:x+y-1
                        for j9=j8+1:x+y-1
                            for j10=j9+1:x+y-1
                                for j11=0:x+y-1
                                  for j12=j1+1:x+y-1
                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1)*f(j5+1)*f(j6+1)*f(j7+1)*f(j8+1) ...,
                        *f(j9+1)*f(j10+1)*f(j11+1)*f(j12+1);
                    lengthsum12=lengthsum12+mid5;
                                  end
                                end
                            end
                        end
                    end
                end
            end
        end
      end
    end
  end
end
end
function lengthsum13=lengthsum13(x,y)
global f
lengthsum13=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1
        for j5=j4+1:x+y-1
            for j6=j5+1:x+y-1
                for j7=j6+1:x+y-1
                    for j8=j7+1:x+y-1
                        for j9=j8+1:x+y-1
                            for j10=j9+1:x+y-1
                                for j11=0:x+y-1
                                  for j12=j1+1:x+y-1
                                    for j13=j2+1:x+y-1
                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1)*f(j5+1)*f(j6+1)*f(j7+1)*f(j8+1) ...,
                        *f(j9+1)*f(j10+1)*f(j11+1)*f(j12+1)*f(j13+1);
                    lengthsum13=lengthsum13+mid5;
                                        end
                                  end
                                end
                            end
                        end
                    end
                end
            end
        end
      end
    end
  end
end
end
function lengthsum14=lengthsum14(x,y)
global f
lengthsum14=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1
        for j5=j4+1:x+y-1
            for j6=j5+1:x+y-1
                for j7=j6+1:x+y-1
                    for j8=j7+1:x+y-1
                        for j9=j8+1:x+y-1
                            for j10=j9+1:x+y-1
                                for j11=0:x+y-1
                                  for j12=j1+1:x+y-1
                                    for j13=j2+1:x+y-1
                                      for j14=j3+1:x+y-1
                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1)*f(j5+1)*f(j6+1)*f(j7+1)*f(j8+1) ...,
                        *f(j9+1)*f(j10+1)*f(j11+1)*f(j12+1)*f(j13+1)*f(j14+1);
                    lengthsum14=lengthsum14+mid5;
                                        end
                                      end
                                    end
                                  end
                                end
                            end
                        end
                    end
                end
            end
        end
      end
    end
  end
end
function lengthsum15=lengthsum15(x,y)
global f
lengthsum15=0;
for j1=0:x+y-1
  for j2=j1+1:x+y-1
    for j3=j2+1:x+y-1
      for j4=j3+1:x+y-1
        for j5=j4+1:x+y-1
            for j6=j5+1:x+y-1
                for j7=j6+1:x+y-1
                    for j8=j7+1:x+y-1
                        for j9=j8+1:x+y-1
                            for j10=j9+1:x+y-1
                                for j11=0:x+y-1
                                  for j12=j1+1:x+y-1
                                    for j13=j2+1:x+y-1
                                      for j14=j3+1:x+y-1
                                        for j15=j4+1:x+y-1
                    mid5=f(j1+1)*f(j2+1)*f(j3+1)*f(j4+1)*f(j5+1)*f(j6+1)*f(j7+1)*f(j8+1) ...,
                        *f(j9+1)*f(j10+1)*f(j11+1)*f(j12+1)*f(j13+1)*f(j14+1)*f(j15+1);
                    lengthsum15=lengthsum15+mid5;
                                        end
                                      end
                                    end
                                  end
                                end
                            end
                        end
                    end
                end
            end
        end
      end
    end
  end
end
end

function tot1=tot1(x,y)
global Fs1
global Fs2
global Fd1
global Fd2
global ktot
global e10
global e20
global v10
global v20
global l
tot1=1;
for i=0:x+y
    F=i*ktot*l;
    a1=v10/l*(1-F/Fs1);
    a2=v20/l*(1-F/Fs2);
    if F>Fs1
        a1=10^(-8);
    end
    if F>Fs2
        a2=10^(-8);
    end
    e1=e10*exp(abs(F)/Fd1);
    e2=e20*exp(abs(F)/Fd2);
    tot1=tot1*(a1*l/v10)/(a1+a2+e1+e2);
end
tot1=tot1*e2/(a1*l/v10);
if(isinf(e1) || isinf(e2))
    tot1=0;
end
end




