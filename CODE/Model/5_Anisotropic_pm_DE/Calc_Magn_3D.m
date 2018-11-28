function [M_x,M_y,M_z,M_av]=Calc_Magn_3D(Ei,T,B,Theta,Phi)
k_B=1.38E-23;
[Bg,Thg,Phg]=meshgrid(B,Theta,Phi);


Z=squeeze(sum(exp(-Ei/T),1));     %Z=sum(-Ei/kB*T), Z(B,Theta,Phi)

[dZr,dZt,dZp] = gradient(permute(log(Z),[2,1,3]),B,Theta,Phi); %[dZ/dB,dZ/dTheta,dZ/dPhi](B,Theta,Phi)
%dZt = dZt./(repmat(B(:)',length(Theta),1,length(Phi)));
%dZp = dZp./permute(repmat(sin(Theta(:)'),length(B),1,length(Phi)),[2,1,3]);
dZt=dZt./Bg;
dZp=dZp./(Bg.*sin(Thg));
%dZp(isnan(dZp))=0;

dlnZ_dBx=dZr.*sin(Thg).*cos(Phg)+dZt.*cos(Thg).*cos(Phg)-dZp.*sin(Thg);
dlnZ_dBy=dZr.*sin(Thg).*sin(Phg)+dZt.*cos(Thg).*sin(Phg)+dZp.*cos(Thg);
dlnZ_dBz=dZr.*cos(Thg)-dZt.*sin(Thg);


%%
%Magnetization in the polar reference frame in mu_B units 
Mx=permute(dlnZ_dBx,[1,3,2])*k_B*T/9.27e-24; %Mx(Theta,Phi,B)
My=permute(dlnZ_dBy,[1,3,2])*k_B*T/9.27e-24; %Mx(Theta,Phi,B)
Mz=permute(dlnZ_dBz,[1,3,2])*k_B*T/9.27e-24; %Mz(Theta,Phi,B)

%%

M_z=squeeze(Mz(1,1,:)); %Hard axis magnetization Mz(theta=0,phi=0°,B)
M_x=squeeze(Mx(end,1,:)); %Easy-plane magnetization Mx(theta=90°,phi=0°,B)
M_y=squeeze(My(end,end,:)); %Easy-plane magnetization My(theta=90°,phi=90°,B)

Mz_B=Mx.*(sin(Theta')).*cos(Phi)+My.*(sin(Theta')).*sin(Phi)+Mz.*cos(Theta');
M_av=4/pi^2*squeeze(trapz(Theta,(trapz(Phi,Mz_B,2)),1));

end