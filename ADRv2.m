% This code is for the simulation in Subsection 8.3 of the paper 
% "Inference on the Structure of Gene Regulatory Networks"
% arXiv:2107.13099

% This code outputs four files, which contain the four performance measures
% for the method of ancestor-descendant relations.


close all
clear all
n=10;
%p=0.04;


rep=10000;
pir=zeros(51,51);
pio=zeros(51,51);
nir=zeros(51,51);
nio=zeros(51,51);

for z1=1:51
    z1
    for z2=1:51
        C0=zeros(3,3);
        %z2
        per=z1/100-0.01;
        ner=z2/100-0.01;

        counter1=0;
        counter2=0;
        counter3=0;
        counter4=0;
        for z3=1:rep
            
%generate connectivity matrix A            
            A=eye(n);
            cc=randi([1*n 3*n]);
            for i=1:cc
                
                An=A^n;
                r=rand;
                
                c=0;
                while c<100
                    c=c+1;
                    r1=randi(n);
                    r2=randi(n);
                    if An(r1,r2)==0 && A(r2,r1)==0
                        A(r2,r1)=1;
                        break;
                    end
                end
            end
            
%generate ancestor-descendant matrix B            
            B=A^n;
            for i=1:n
                for j=1:n
                    if B(i,j)>0
                        B(i,j)=1;
                    end
                end
            end
            
%reconstruct connectivity matrix R            
            R=-ones(n,n);
            for i=1:n
                for j=1:n
                    if i==j
                        R(i,j)=1;
                        continue
                    end
                    if B(i,j)==0
                        R(i,j)=0;
                    else
                        te=B(i,:)*B(:,j);
                        if te>2
                            R(i,j)=0.5;
                        else
                            R(i,j)=1;
                        end
                    end
                end
            end
            
%pollute B            
            B1=B;
            for i=1:n
                for j=1:n
                    if i==j
                        continue
                    end
                    if B(i,j)==0
                        te=rand;
                        if te<ner
                            B1(i,j)=1;
                        end
                    else
                        te=rand;
                        if te<per
                            B1(i,j)=0;
                        end
                    end
                end
            end
            
%reconstruct R1 with B1            
            R1=-ones(n,n);
            for i=1:n
                for j=1:n
                    if i==j
                        R1(i,j)=1;
                        continue
                    end
                    if B1(i,j)==0
                        R1(i,j)=0;
                    else
                        te=B1(i,:)*B1(:,j);
                        if te>2
                            R1(i,j)=0.5;
                        else
                            R1(i,j)=1;
                        end
                    end
                end
            end
            
            
            pc=0;%R=R1=1
            nc=0;%R=R1=0
            rp=0;%R=1
            rn=0;%R=0
            op=0;%R1=1
            on=0;%R1=0
            for i=1:n
                for j=1:n
                    if i==j
                        continue;
                    end
                    if R(i,j)==0
                        rn=rn+1;
                    else
                        if R(i,j)==1
                            rp=rp+1;
                        end
                    end
                    if R1(i,j)==0
                        on=on+1;
                    else
                        if R1(i,j)==1
                            op=op+1;
                        end
                    end
                    if R(i,j)==1 && R1(i,j)==1
                        pc=pc+1;
                    end
                    if R(i,j)==0 && R1(i,j)==0
                        nc=nc+1;
                    end
                end
            end
            if rp==0
                rp=1;
            end
            if rn==0
                rn=1;
            end
            if op==0
                op=1;
            end
            if on==0
                on=1;
            end
            
            counter1=counter1+pc/rp; 
            counter2=counter2+pc/op; 
            counter3=counter3+nc/rn; 
            counter4=counter4+nc/on; 

%             C=zeros(3,3);
%             for i=1:n
%                 for j=1:n
%                     if i==j
%                         continue;
%                     end
%                     if R(i,j)==1
%                         if R1(i,j)==1
%                             C(1,1)=C(1,1)+1;
%                         end
%                         if R1(i,j)==0.5
%                             C(1,2)=C(1,2)+1;
%                         end
%                         if R1(i,j)==0
%                             C(1,3)=C(1,3)+1;
%                         end
%                     end
%                     if R(i,j)==0.5
%                         if R1(i,j)==1
%                             C(2,1)=C(2,1)+1;
%                         end
%                         if R1(i,j)==0.5
%                             C(2,2)=C(2,2)+1;
%                         end
%                         if R1(i,j)==0
%                             C(2,3)=C(2,3)+1;
%                         end
%                     end
%                     if R(i,j)==0
%                         if R1(i,j)==1
%                             C(3,1)=C(3,1)+1;
%                         end
%                         if R1(i,j)==0.5
%                             C(3,2)=C(3,2)+1;
%                         end
%                         if R1(i,j)==0
%                             C(3,3)=C(3,3)+1;
%                         end
%                     end
%                 end
%             end
%             C0=C0+C;
        end
        pir(z1,z2)=counter1/rep;
        pio(z1,z2)=counter2/rep;
        nir(z1,z2)=counter3/rep;
        nio(z1,z2)=counter4/rep;
        %C0/rep
    end
end
% subplot(2,2,1)
% contourf(0:0.01:0.5,0:0.01:0.5,pir)
% xlabel('ner') 
% ylabel('per') 
% title('SEN')
% subplot(2,2,2)
% contourf(0:0.01:0.5,0:0.01:0.5,pio)
% xlabel('ner') 
% ylabel('per') 
% title('PPV')
% subplot(2,2,3)
% contourf(0:0.01:0.5,0:0.01:0.5,nir)
% xlabel('ner') 
% ylabel('per')
% title('SPE')
% subplot(2,2,4)
% contourf(0:0.01:0.5,0:0.01:0.5,nio)
% xlabel('ner') 
% ylabel('per') 
% title('NPV')


save adrpir.dat pir -ascii 
save adrpio.dat pio -ascii 
save adrnir.dat nir -ascii 
save adrnio.dat nio -ascii 








%save sradr100.dat sr -ascii


