% This code is for the simulation in Subsection 8.4 of the paper 
% "Inference on the Structure of Gene Regulatory Networks"
% arXiv:2107.13099

% This code outputs four files, which contain the four performance measures
% for the method of combining conditional independence and 
% ancestor-descendant relations.

close all
clear all

rep=100000;
pir=zeros(51,51);
pio=zeros(51,51);
nir=zeros(51,51);
nio=zeros(51,51);
for z1=1:51
    z1
    for z2=1:51
        p=z1/100-0.01;
        q=z2/100-0.01;
        counter1=0;
        counter2=0;
        counter3=0;
        counter4=0;
        for i=1:rep
            
            % initialize DAG A, conditional ind C, ADR matrix R
            % C=[12 12|3
            %    23 23|1
            %    31 31|2]
            % R=[12 21
            %    23 32
            %    31 13]
            zb=randi([1 4]);
            if zb==1
                A=[-1 1 1
                    0 -1 1
                    0 0 -1];
                C=[0 0
                    0 0
                    0 0];
                R=[1 0
                    1 0
                    0 1];
            end
            if zb==2
                A=[-1 1 0
                    0 -1 1
                    0 0 -1];
                C=[0 0
                    0 0
                    0 1];
                R=[1 0
                    1 0
                    0 1];
            end
            if zb==3
                A=[-1 0 1
                    0 -1 1
                    0 0 -1];
                C=[1 0
                    0 0
                    0 0];
                R=[0 0
                    1 0
                    0 1];
            end
            if zb==4
                A=[-1 1 1
                    0 -1 0
                    0 0 -1];
                C=[0 0
                    0 1
                    0 0];
                R=[1 0
                    0 0
                    0 1];
            end
            
            % pollute C, R
            for w1=1:3
                for w2=1:2
                    if C(w1,w2)==1
                        te=rand;
                        if te<p
                            C(w1,w2)=0;
                        end
                    else
                        te=rand;
                        if te<q
                            C(w1,w2)=1;
                        end
                    end
                    if R(w1,w2)==1
                        te=rand;
                        if te<q
                            R(w1,w2)=0;
                        end
                    else
                        te=rand;
                        if te<p
                            R(w1,w2)=1;
                        end
                    end
                end
            end
            
            
            
            % infer A
            A2=-ones(3,3);
            if sum(C(1,:))>0
                A2(1,2)=0;
                A2(2,1)=0;
            else
                if R(1,1)==R(1,2)
                    te=rand;
                    if te<0.5
                        A2(1,2)=1;
                        A2(2,1)=0; 
                    else
                        A2(1,2)=0;
                        A2(2,1)=1; 
                    end
                else
                    if R(1,1)==1
                        A2(1,2)=1;
                        A2(2,1)=0;
                    else
                        A2(1,2)=0;
                        A2(2,1)=1;
                    end
                end
            end
            if sum(C(2,:))>0
                A2(2,3)=0;
                A2(3,2)=0;
            else
                if R(2,1)==R(2,2)
                    te=rand;
                    if te<0.5
                        A2(2,3)=1;
                        A2(3,2)=0; 
                    else
                        A2(2,3)=0;
                        A2(3,2)=1; 
                    end
                else
                    if R(2,1)==1
                        A2(2,3)=1;
                        A2(3,2)=0;
                    else
                        A2(2,3)=0;
                        A2(3,2)=1;
                    end
                end
            end
            if sum(C(3,:))>0
                A2(3,1)=0;
                A2(1,3)=0;
            else
                if R(3,1)==R(3,2)
                    te=rand;
                    if te<0.5
                        A2(3,1)=1;
                        A2(1,3)=0; 
                    else
                        A2(3,1)=0;
                        A2(1,3)=1; 
                    end
                else
                    if R(3,1)==1
                        A2(3,1)=1;
                        A2(1,3)=0;
                    else
                        A2(3,1)=0;
                        A2(1,3)=1;
                    end
                end
            end
            
            
            % calculate measures
            pp=0;
            pn=0;
            np=0;
            nn=0;
            for s1=1:3
                for s2=1:3
                    if s1==s2
                        continue;
                    end
                    if A(s1,s2)==1
                        if A2(s1,s2)==1
                            pp=pp+1;
                        else
                            pn=pn+1;
                        end
                    else
                        if A2(s1,s2)==1
                            np=np+1;
                        else
                            nn=nn+1;
                        end
                    end
                end
            end
            counter1=counter1+pp;
            counter2=counter2+pn;
            counter3=counter3+np;
            counter4=counter4+nn;
        end
        pir(z1,z2)=counter1/(counter1+counter2);
        pio(z1,z2)=counter1/(counter1+counter3);
        nir(z1,z2)=counter4/(counter3+counter4);
        nio(z1,z2)=counter4/(counter2+counter4);
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

save cipir.dat pir -ascii
save cipio.dat pio -ascii
save cinir.dat nir -ascii
save cinio.dat nio -ascii














