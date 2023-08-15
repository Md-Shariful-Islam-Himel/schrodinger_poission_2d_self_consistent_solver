
function [N, R] = poission(x, del_x, Na, Vi, Vl, nz)
%DEFINING CONSTANTS AND MATRIX
e=1.602*10^-19;
epsilon=8.854e-12;
m=round(x/del_x);
A=eye(m-2, m-2);



%N MATRIX FORMATION
A(boolean(eye(m-2)))=-2;
did1 = 2:m-1:(numel(A));
did2 = m-1:m-1:(numel(A));
for l = [1:m-3]
    A(did1(l)) = 1;
    A(did2(l)) = 1;
end

A= A/(del_x)^2;




%Q MATRIX FORMATION
Q=zeros((m-2), 1);
R=zeros((m-2), 1);
oxide=round(0.25*(m-2));


for k = [1:oxide]
    er=3.9;
    Q(k)=(e.*(-nz(k)-Na(k)))/-(epsilon*er);
    R(k)= e.*(nz(k)+Na(k));
end

for k = [(oxide+1):((m-2)-oxide)]
    er=11.7;
    Q(k)=(e.*(-nz(k)-Na(k)))/-(epsilon*er);
    R(k)= e.*(nz(k)+Na(k));
end

for k = [(((m-2)-oxide)+1):(m-2)]
    er=3.9;
    Q(k)=(e.*(-nz(k)-Na(k)))/-(epsilon*er);
    R(k)= e.*(nz(k)+Na(k));
end



Q(1)=Q(1)+(Vi/((del_x)^2));
Q(m-3)= Q(m-3)+(Vl/((del_x)^2));

P= linsolve(A, Q);
N= [Vi; P; Vl];

end