% This code is for the simulation in Subsection 8.5 of the paper 
% "Inference on the Structure of Gene Regulatory Networks"
% arXiv:2107.13099

% This code just shows that the partially observed ODE method does not
% work in reality with a small noise. Required parameters cannot be solved.



% solve for the original function
fun = @r1;
x0 = [2,1.5,-0.5,2.5,0,3.5];
%x0=[1.5 1.001 -1 2.001 -0.5 3.001];
x1 = fsolve(fun,x0)

fun2 = @r2;
y0 = [1.5,1.5,-0.5,2.5,0.5,3.5];
%y0=[1 1.001 -1 2.001 0 3.001];
y1 = fsolve(fun2,y0)



% solve for the function after Laplace transform
fun = @r3;
z0 = [2,1.5,-0.5,2.5,0,3.5];
%x0=[1.5 1.001 -1 2.001 -0.5 3.001];
x2 = fsolve(fun,z0)

fun2 = @r4;
w0 = [1.5,1.5,-0.5,2.5,0.5,3.5];
%y0=[1 1.001 -1 2.001 0 3.001];
y2 = fsolve(fun2,w0)





function F = r1(x)
err=0.01;
for z=1:6
    s1(z)=3/2*exp(-z)-exp(-2*z)-1/2*exp(-3*z)+7/3;
end
for i=1:6
    F(i) = x(1)*exp(-x(2)*i)+x(3)*exp(-x(4)*i)+x(5)*exp(-x(6)*i)+7/3-s1(i)*(1+err*(2*rand-1));
end
end

function F = r2(x)
err=0.01;
for z=1:6
    s2(z)=exp(-z)-exp(-2*z)+7/3;
end
for i=1:6
    F(i) = x(1)*exp(-x(2)*i)+x(3)*exp(-x(4)*i)+x(5)*exp(-x(6)*i)+7/3-s2(i)*(1+err*(2*rand-1));
end
end




function F = r3(x)
err=0.01;
for z=1:6
    s1(z)=3/2/(z+1)-1/(z+2)-1/2/(z+3)+7/3/z;
end
for i=1:6
    F(i) = x(1)/(i+x(2))+x(3)/(i+x(4))+x(5)/(i+x(6))+7/3/i-s1(i)*(1+err*(2*rand-1));
end
end

function F = r4(x)
err=0.01;
for z=1:6
    s2(z)=1/(z+1)-1/(z+2)+7/3/z;
end
for i=1:6
    F(i) = x(1)/(i+x(2))+x(3)/(i+x(4))+x(5)/(i+x(6))+7/3/i-s2(i)*(1+err*(2*rand-1));
end
end









