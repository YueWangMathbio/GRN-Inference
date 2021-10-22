% This code is for the simulation in Subsection 8.1 of the paper 
% "Inference on the Structure of Gene Regulatory Networks"
% arXiv:2107.13099

% This code outputs one file, which contains the success rate
% for the method of path blocking with gene expression data.






%%%%%%%%%%% path blocking with gene expression data

clear all
close all
n=1000000;


%to=0.1;


for k=1:51
    k
    for j=1:101
        c1=0;%true po
        c2=0;%false po
        c3=0;%false ne
        c4=0;%true ne
        err=k/100-0.01;
        to=j/100-0.01;
        for i=1:n
            
            b=rand(2,1);
            A=-ones(2,2);
            p=rand;
            q=rand;
            if q<p
                A(1,2)=rand*2-1;
            else
                A(1,2)=0;
            end
            r=rand;
            if r<p
                A(2,1)=rand*2-1;
            else
                A(2,1)=0;
            end
            s=-(A\b);
            x1=s(2);
            x2=b(2);
            
            te=rand*2-1;
            y1=x1*(1+te*err);
            te2=rand*2-1;
            y2=x2*(1+te2*err);
            m1=min(y1,y2);
            m2=max(y1,y2);
            if m1/m2>=1-to
                if abs(x2-x1)<1e-10
                    c4=c4+1;
                else
                    c3=c3+1;
                end
            else
                if abs(x2-x1)>1e-10
                    c1=c1+1;
                else
                    c2=c2+1;
                end
            end
        end
        
        trate(k,j)=(c1+c4)/n;
    end
end
%hold on
% contourf(0:0.01:1,0:0.01:0.5,trate)
% xlabel('tolerate threshold')
% ylabel('error rate')
save trate.dat trate -ascii
