close all
n=6;
v=zeros(n-1, 1)
for i=1:n-1
v(i)=i;
end
v1=zeros(n, 1)
for i=1:n
v1(i)=i;
end
[Charge] = Self_Consistent_Solver_Function_for_CV(n);
close all
figure
plot (0.1*v1, real(Charge));
xlabel ("Gate Voltage");
ylabel ("Charge");

c=diff(real(Charge));
figure
plot (0.1*v, c);
xlabel ("Gate Voltage");
ylabel ("Capacitance");
