function [shi, E]= schrodinger(x, del_x, v)
hbar = 1.05*10^-34;
mass = 9.1*10^-31;
number_of_data_points = round(x/del_x);

kinetic = zeros(number_of_data_points, 1);

for i = 1:number_of_data_points
    kinetic(i,i) = 2;
    if i > 1
        kinetic(i, i-1) = -1;
        kinetic(i-1, i) = -1;
    end
end


t = hbar^2/(2*mass*(del_x)^2);
kinetic = kinetic*t;
k=1;

 for (i=1: (0.2*number_of_data_points))
 if(length(v))>1000
     v(1, :)=[];
     v((length(v)),: )=[];
 end
 end
h = kinetic + v;
figure
[shi, E] = eig(h);
end