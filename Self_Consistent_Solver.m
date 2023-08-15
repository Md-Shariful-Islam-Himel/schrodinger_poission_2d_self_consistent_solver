clc
clear
close all


e=1.602*10^-19;
c=0;
x=1*10^-8;
del_x=1*10^-11;
m=round(x/del_x);
a=1e-10;
mass=9.1e-31;
red_w= a*[1:m];
w=a*(1:round(1.2*(x/del_x)));
nz=zeros(((round(x/del_x))-2), 1);
Vprev=0;
Vnew=100;
K=1.3807e-23;
hbar = 1.05*10^-34;
T=300;


Na=ones((round(x/del_x)), 1);
Na=1*10^22 * Na;
Vi=2;
Vg=e*Vi;
Vl=2;
Na((round(x/del_x)), :)=[];
Na(1, :)=[];
[v, R]=poission(x, del_x, Na, Vi, Vl, nz);

reduced_red_w=red_w;
reduced_red_w(m)=[];
reduced_red_w(1)=[];
subplot(1,3,1)
plot(reduced_red_w, R);
title ("Charge density at 1st run");
xlabel ("Z");
ylabel ("Charge density");

O=Vi+zeros(0.1*m, 1);
Q=Vl+zeros(0.1*m, 1);
Potential= [O; v; Q];
subplot(1,3,2)
plot(w, Potential);
title ("Potential profile at 1st run");
xlabel ("Z");
ylabel ("Potential");


Energy=-e.*Potential;
for(i=(round(0.22*length(w))):(round(0.783*length(w))))
    Energy(i)=Energy(i)-3.1e-19;
end
subplot(1,3,3)
plot (w, Energy);
title ("Energy band diagram at 1st run");
xlabel ("Z");
ylabel ("Energy");

v=-e.*v;
[shi, E]=schrodinger(x, del_x, v);
subplot(1,3,1)
plot (red_w, shi(:, 1));
hold on
plot (red_w, shi(:, 2));
hold on
plot (red_w, shi(:, 3));
hold on
title ("Wavefunction(Shi) at 1st run");
xlabel ("Z");
ylabel ("Shi");

subplot(1,3,3)
plot (red_w, E(:, 1));
hold on
plot (red_w, E(:, 2));
hold on
plot (red_w, E(:, 3));
title ("E at 1st run");
xlabel ("Z");
ylabel ("E");

shisq=shi^2;
norm_shi=shisq/trapz(shisq);
O=zeros(0.1*m, 1);
Q=zeros(0.1*m, 1);
nshi= [O; norm_shi; Q];
subplot(1,3,2)
plot (w, nshi);
title ("Normalized wavefunction at 1st run");
xlabel ("Z");
ylabel ("Shi(normalized)");


while (((Vnew-Vprev).^2)>1e-200)
    c=c+1;
    if (c==1)
    v=v/-e;
    end
    Vprev=v;
    set(groot,'defaultFigureVisible','off');
    Vi=v(1);
    Vl=v(round(x/del_x));
    v(round((x/del_x)), :)=[];
    v(1, :)=[];
reduced_norm_shi=norm_shi;
reduced_norm_shi(round((x/del_x)), :)=[];
reduced_norm_shi(1, :)=[];

reduced_E=E;
reduced_E(round((x/del_x)), :)=[];
reduced_E(1, :)=[];

for j=1:6
    nij=((K*T*mass)/(3.1416*hbar^2)).*log(1+ exp((0+Vg-(reduced_E(:, j)))/(K*T)));
    nz=nz+(nij.*reduced_norm_shi);
end
[v, Q]=poission(x, del_x, Na, Vi, Vl, nz);
v=-e.* v;
[shi, E]=schrodinger(x, del_x, v);
shisq=shi^2;
norm_shi=shisq/trapz(shisq);
Vnew=v;
end


set(groot,'defaultFigureVisible','on');

figure
%subplot (1, 3, 1)
plot(reduced_red_w,Q);
title ("Charge density after " +c+ " iterations");
xlabel ("Z");
ylabel ("Charge density");


O=Vi+zeros(0.1*m, 1);
Y=Vl+zeros(0.1*m, 1);
Energy= [O; v; Y];
Volt=Energy;
for i=(round(0.22*length(w))):(round(0.783*length(w)))
    Energy(i)=Energy(i)-3.1e-23;
end
for i=1:(round(0.18*length(w)))
    Energy(i)=0;
end
for i=(round(0.82*length(w))): length(w)
    Energy(i)=0;
end
%subplot(1, 3, 3)
figure
plot(w, Energy);
title ("Energy band diagram after " +c+ " iterations");
xlabel ("Z");
ylabel ("Energy");

v=Volt/-e;
%subplot(1, 3, 2)
figure
plot(w, v);
title ("Potential after " +c+ " iterations");
ylim([-3e-4 0.4e-4]);
xlabel ("Z");
ylabel ("Potential");

figure
%subplot(1,3,1)
plot (red_w, shi(:, 1));
hold on
plot (red_w, shi(:, 2));
hold on
plot (red_w, shi(:, 3));
hold on
title ("Wavefunction(Shi) after " +c+ " iterations");
xlabel ("Z");
ylabel ("Shi");

subplot(1,3,3)
plot (red_w, E(:, 1));
hold on
plot (red_w, E(:, 2));
hold on
plot (red_w, E(:, 3));
title ("E at after " +c+ " iterations");
xlabel ("Z");
ylabel ("E");


O=zeros(0.1*m, 1);
U=zeros(0.1*m, 1);
nshi= [O; norm_shi; U];
%subplot(1,3,2)
figure
plot (w, nshi);
title ("Normalized wavefunction after " +c+ " iterations");
xlabel ("Z");
ylabel ("Shi(Normalized)");


figure
plot (reduced_red_w, nz);
title ("inversion charge after " +c+ " iterations");
xlabel ("Z");
ylabel ("Inversion Charge Density");
