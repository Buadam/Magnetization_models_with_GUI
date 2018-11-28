function [Mx_T,Mz_T,Mav_T]=Calc_MT_2D(Ei,T,B,Theta,B0)
k_B=1.38E-23;

Z=squeeze(sum(exp(-Ei/T),1));     %Z=sum(-Ei/kB*T), Z(B,Theta)

[dZr,dZt] = gradient(log(Z)',B,Theta);
dZt = dZt./(repmat(B(:)',length(Theta),1));

dlnZ_dBz=dZr'.*cos(Theta)-dZt'.*sin(Theta);
dlnZ_dBx=dZr'.*sin(Theta)+dZt'.*cos(Theta);

%%
%Magnetization in the polar reference frame in mu_B units 
Mx=k_B*T*dlnZ_dBx/9.27e-24; %Mx(B,Theta)
Mz=k_B*T*dlnZ_dBz/9.27e-24; %Mz(B,Theta)

%%
M_z=Mz(:,1); %Hard axis magnetization Mz(B,theta=0)
M_x=Mx(:,end); %Easy-plane magnetization Mx(B,theta=90°)
M_av=2/pi*(trapz(Theta,(Mx.*sin(Theta)+Mz.*cos(Theta)),2)); %averaging to the polar angle

Mx_T=interp1(B,M_x,B0);
Mz_T=interp1(B,M_z,B0);
Mav_T=interp1(B,M_av,B0);


end