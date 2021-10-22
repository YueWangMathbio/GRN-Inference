% This code is for the simulation in Section 10 of the paper 
% "Inference on the Structure of Gene Regulatory Networks"
% arXiv:2107.13099

% This code outputs one file, which contains the success rate
% for the method to determine autocatalysis.


clear all
close all

A=zeros(101,100);
k=5;% base growth
c=2;% degrade
for i=1:101
    b=i/100-0.01;%extra growth
    A(i,1)=1;
    for j=2:100
        A(i,j)=A(i,j-1)*(k+b*(j-2))/c/(j-1);
    end
    su=sum(A(i,:));
    A(i,:)=A(i,:)/su;% stationary distribution of X
end

C=zeros(101,2); % help to calculate VMR
for i=1:101
    q=A(i,:);
    sl=0:99;
    m1=sum(sl.*q);
    m2=sum(sl.^2.*q);
    sig=m2-m1^2;
    cen=sl-m1;
    cm4=sum(cen.^4.*q);
    vv=cm4-sig^2;
    C(i,1)=log(sig/m1);
    C(i,2)=sqrt(vv/sig+sig/m1);
end

n1=15;
n2=501;
y=zeros(n1,n2);
for z1=1:n1 %11    
    for z2=1:n2 %51        
        xx=zeros(1,101); % success rate
        n=32*2^z1; % sample size
        T=(z2-1)/(n2-1)*0.50; % threshold
        for i=1:1
            bo=sqrt(n)*(log(1+T)-C(i,1))/C(i,2);
            xx(i)=normcdf(bo);
        end
        for i=2:101
            bo=sqrt(n)*(log(1+T)-C(i,1))/C(i,2);
            xx(i)=1-normcdf(bo);
        end
        y(z1,z2)=xx(1)/2+sum(xx(2:101))/200; % overall success rate
    end
end


save acf.dat y -ascii














