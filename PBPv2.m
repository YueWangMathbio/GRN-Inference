% This code is for the simulation in Subsection 8.2 of the paper 
% "Inference on the Structure of Gene Regulatory Networks"
% arXiv:2107.13099

% This code outputs four files, which contain the four performance measures
% for the method of path blocking with phenotype data.



clear all
close all

rep=100000;
sr=zeros(51,51);
for z1=1:51
    z1
    for z2=1:51
        C0=zeros(3,3);
        per=z1/100-0.01;
        ner=z2/100-0.01;
        counter=0;
        for z3=1:rep
            
            p=rand;
            A=zeros(4,4);
            %A0=A^3;
            %while A0(1,4)*A0(2,4)*A0(3,4)<1
            A=zeros(4,4);
            for i=1:3
                for j=1:4
                    te=rand;
                    if te<p
                        A(i,j)=1;
                    end
                end
            end
            for i=1:4
                A(i,i)=1;
            end
            
            %A
            
            A0=A^3;
            %end
            
            % if A0(1,4)*A0(2,4)*A0(3,4)<1
            % break
            %A0=A^3;
            A1=A;
            A1(1,:)=0;
            A1(:,1)=0;
            
            A2=A;
            A2(2,:)=0;
            A2(:,2)=0;
            
            A3=A;
            A3(3,:)=0;
            A3(:,3)=0;
            
            A12=A2;
            A12(1,:)=0;
            A12(:,1)=0;
            
            A13=A3;
            A13(1,:)=0;
            A13(:,1)=0;
            
            A23=A2;
            A23(3,:)=0;
            A23(:,3)=0;
            
            
            A1=A1^3;
            A2=A2^3;
            A3=A3^3;
            A12=A12^3;
            A13=A13^3;
            A23=A23^3;
            
            B1=-ones(1,4); % 2, 3, 2+3, emptyset % 1 means blocking, 0 means not blocking
            if A0(1,4)*A2(1,4)>0
                B1(1)=0;
            else
                B1(1)=1;
            end
            if A0(1,4)*A3(1,4)>0
                B1(2)=0;
            else
                B1(2)=1;
            end
            if A0(1,4)*A23(1,4)>0
                B1(3)=0;
            else
                B1(3)=1;
            end
            if A0(1,4)==0
                B1(4)=1;
            else
                B1(4)=0;
            end
            
            B2=-ones(1,3); % 1, 3, 1+3
            if A0(2,4)*A1(2,4)>0
                B2(1)=0;
            else
                B2(1)=1;
            end
            if A0(2,4)*A3(2,4)>0
                B2(2)=0;
            else
                B2(2)=1;
            end
            if A0(2,4)*A13(2,4)>0
                B2(3)=0;
            else
                B2(3)=1;
            end
            if A0(2,4)==0
                B2(4)=1;
            else
                B2(4)=0;
            end
            
            B3=-ones(1,3); % 1, 2, 1+2
            if A0(3,4)*A1(3,4)>0
                B3(1)=0;
            else
                B3(1)=1;
            end
            if A0(3,4)*A2(3,4)>0
                B3(2)=0;
            else
                B3(2)=1;
            end
            if A0(3,4)*A12(3,4)>0
                B3(3)=0;
            else
                B3(3)=1;
            end
            if A0(3,4)==0
                B3(4)=1;
            else
                B3(4)=0;
            end
            
            BO1=B1;
            BO2=B2;
            BO3=B3;
            
            %[B1;B2;B3]
            
            R=-ones(4,4);
            for i=1:4
                R(4,i)=0;
                R(i,i)=1;
            end
            
            %reconstruct R for V1
            if sum(B1)==0
                R(1,4)=1;
                R(1,2)=0.5;
                R(1,3)=0.5;
            else
                R(1,4)=0;
                if B1(2)==1 && B2(2)==0
                    R(1,2)=0;
                else
                    if B1(1)==1
                        R(1,2)=1;
                    end
                    if B1(2)==0 && B1(3)==1
                        R(1,2)=1;
                    end
                end
                if R(1,2)<0
                    R(1,2)=0.5;
                end
                
                if B1(1)==1 && B3(2)==0
                    R(1,3)=0;
                else
                    if B1(2)==1
                        R(1,3)=1;
                    end
                    if B1(1)==0 && B1(3)==1
                        R(1,3)=1;
                    end
                end
                if R(1,3)<0
                    R(1,3)=0.5;
                end
            end
            if B1(4)==1
                if B2(4)==1
                    R(1,2)=0.5;
                else
                    R(1,2)=0;
                end
                if B3(4)==1
                    R(1,3)=0.5;
                else
                    R(1,3)=0;
                end
            end
            
            
            if sum(B2)==0
                R(2,4)=1;
                R(2,1)=0.5;
                R(2,3)=0.5;
            else
                R(2,4)=0;
                if B2(2)==1 && B1(2)==0
                    R(2,1)=0;
                else
                    if B2(1)==1
                        R(2,1)=1;
                    end
                    if B2(2)==0 && B2(3)==1
                        R(2,1)=1;
                    end
                end
                if R(2,1)<0
                    R(2,1)=0.5;
                end
                
                if B2(1)==1 && B3(1)==0
                    R(2,3)=0;
                else
                    if B2(2)==1
                        R(2,3)=1;
                    end
                    if B2(1)==0 && B2(3)==1
                        R(2,3)=1;
                    end
                end
                if R(2,3)<0
                    R(2,3)=0.5;
                end
            end
            if B2(4)==1
                if B1(4)==1
                    R(2,1)=0.5;
                else
                    R(2,1)=0;
                end
                if B3(4)==1
                    R(2,3)=0.5;
                else
                    R(2,3)=0;
                end
            end
            
            if sum(B3)==0
                R(3,4)=1;
                R(3,1)=0.5;
                R(3,2)=0.5;
            else
                R(3,4)=0;
                if B3(2)==1 && B1(1)==0
                    R(3,1)=0;
                else
                    if B3(1)==1
                        R(3,1)=1;
                    end
                    if B3(2)==0 && B3(3)==1
                        R(3,1)=1;
                    end
                end
                if R(3,1)<0
                    R(3,1)=0.5;
                end
                
                if B3(1)==1 && B2(1)==0
                    R(3,2)=0;
                else
                    if B3(2)==1
                        R(3,2)=1;
                    end
                    if B3(1)==0 && B3(3)==1
                        R(3,2)=1;
                    end
                end
                if R(3,2)<0
                    R(3,2)=0.5;
                end
            end
            if B3(4)==1
                if B1(4)==1
                    R(3,1)=0.5;
                else
                    R(3,1)=0;
                end
                if B2(4)==1
                    R(3,2)=0.5;
                else
                    R(3,2)=0;
                end
            end
            R0=R;
            
            % generate polluted blocking information
            for i=1:4
                bt=B1(i);
                if bt==0
                    te=rand;
                    if te<ner
                        B1(i)=1;
                    end
                else
                    te=rand;
                    if te<per
                        B1(i)=0;
                    end
                end
            end
            
            for i=1:4
                bt=B2(i);
                if bt==0
                    te=rand;
                    if te<ner
                        B2(i)=1;
                    end
                else
                    te=rand;
                    if te<per
                        B2(i)=0;
                    end
                end
            end
            
            for i=1:4
                bt=B3(i);
                if bt==0
                    te=rand;
                    if te<ner
                        B3(i)=1;
                    end
                else
                    te=rand;
                    if te<per
                        B3(i)=0;
                    end
                end
            end
            
            %[B1;B2;B3]
            
            R=-ones(4);
            for i=1:4
                R(4,i)=0;
                R(i,i)=1;
            end
            
            %reconstruct R for V1
            if sum(B1)==0
                R(1,4)=1;
                R(1,2)=0.5;
                R(1,3)=0.5;
            else
                R(1,4)=0;
                if B1(2)==1 && B2(2)==0
                    R(1,2)=0;
                else
                    if B1(1)==1
                        R(1,2)=1;
                    end
                    if B1(2)==0 && B1(3)==1
                        R(1,2)=1;
                    end
                end
                if R(1,2)<0
                    R(1,2)=0.5;
                end
                
                if B1(1)==1 && B3(2)==0
                    R(1,3)=0;
                else
                    if B1(2)==1
                        R(1,3)=1;
                    end
                    if B1(1)==0 && B1(3)==1
                        R(1,3)=1;
                    end
                end
                if R(1,3)<0
                    R(1,3)=0.5;
                end
            end
            if B1(4)==1
                if B2(4)==1
                    R(1,2)=0.5;
                else
                    R(1,2)=0;
                end
                if B3(4)==1
                    R(1,3)=0.5;
                else
                    R(1,3)=0;
                end
            end
            
            
            if sum(B2)==0
                R(2,4)=1;
                R(2,1)=0.5;
                R(2,3)=0.5;
            else
                R(2,4)=0;
                if B2(2)==1 && B1(2)==0
                    R(2,1)=0;
                else
                    if B2(1)==1
                        R(2,1)=1;
                    end
                    if B2(2)==0 && B2(3)==1
                        R(2,1)=1;
                    end
                end
                if R(2,1)<0
                    R(2,1)=0.5;
                end
                
                if B2(1)==1 && B3(1)==0
                    R(2,3)=0;
                else
                    if B2(2)==1
                        R(2,3)=1;
                    end
                    if B2(1)==0 && B2(3)==1
                        R(2,3)=1;
                    end
                end
                if R(2,3)<0
                    R(2,3)=0.5;
                end
            end
            if B2(4)==1
                if B1(4)==1
                    R(2,1)=0.5;
                else
                    R(2,1)=0;
                end
                if B3(4)==1
                    R(2,3)=0.5;
                else
                    R(2,3)=0;
                end
            end
            
            if sum(B3)==0
                R(3,4)=1;
                R(3,1)=0.5;
                R(3,2)=0.5;
            else
                R(3,4)=0;
                if B3(2)==1 && B1(1)==0
                    R(3,1)=0;
                else
                    if B3(1)==1
                        R(3,1)=1;
                    end
                    if B3(2)==0 && B3(3)==1
                        R(3,1)=1;
                    end
                end
                if R(3,1)<0
                    R(3,1)=0.5;
                end
                
                if B3(1)==1 && B2(1)==0
                    R(3,2)=0;
                else
                    if B3(2)==1
                        R(3,2)=1;
                    end
                    if B3(1)==0 && B3(3)==1
                        R(3,2)=1;
                    end
                end
                if R(3,2)<0
                    R(3,2)=0.5;
                end
            end
            if B3(4)==1
                if B1(4)==1
                    R(3,1)=0.5;
                else
                    R(3,1)=0;
                end
                if B2(4)==1
                    R(3,2)=0.5;
                else
                    R(3,2)=0;
                end
            end
            
            R1=R;
            %R
            R=R0;
            C=zeros(3,3);
            for i=1:3
                for j=1:4
                    if i==j
                        continue;
                    end
                    if R(i,j)==1
                        if R1(i,j)==1
                            C(1,1)=C(1,1)+1;
                        end
                        if R1(i,j)==0.5
                            C(1,2)=C(1,2)+1;
                        end
                        if R1(i,j)==0
                            C(1,3)=C(1,3)+1;
                        end
                    end
                    if R(i,j)==0.5
                        if R1(i,j)==1
                            C(2,1)=C(2,1)+1;
                        end
                        if R1(i,j)==0.5
                            C(2,2)=C(2,2)+1;
                        end
                        if R1(i,j)==0
                            C(2,3)=C(2,3)+1;
                        end
                    end
                    if R(i,j)==0
                        if R1(i,j)==1
                            C(3,1)=C(3,1)+1;
                        end
                        if R1(i,j)==0.5
                            C(3,2)=C(3,2)+1;
                        end
                        if R1(i,j)==0
                            C(3,3)=C(3,3)+1;
                        end
                    end
                end
            end
            C0=C0+C;
        end
        C0=C0/rep;
        pir(z1,z2)=C0(1,1)/sum(C0(1,:));
        pio(z1,z2)=C0(1,1)/sum(C0(:,1));
        nir(z1,z2)=C0(3,3)/sum(C0(3,:));
        nio(z1,z2)=C0(3,3)/sum(C0(:,3));
    end
end
%contourf(sr)




%save sr.dat sr -ascii

subplot(2,2,1)
contourf(0:0.01:0.5,0:0.01:0.5,pir)
xlabel('ner') 
ylabel('per') 
title('SEN')
subplot(2,2,2)
contourf(0:0.01:0.5,0:0.01:0.5,pio)
xlabel('ner') 
ylabel('per') 
title('PPV')
subplot(2,2,3)
contourf(0:0.01:0.5,0:0.01:0.5,nir)
xlabel('ner') 
ylabel('per')
title('SPE')
subplot(2,2,4)
contourf(0:0.01:0.5,0:0.01:0.5,nio)
xlabel('ner') 
ylabel('per') 
title('NPV')

save pbppir.dat pir -ascii
save pbppio.dat pio -ascii
save pbpnir.dat nir -ascii
save pbpnio.dat nio -ascii
% n=1e5




