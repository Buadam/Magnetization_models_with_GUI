function OUTPUT = anis_model_Spherical_3D(T,param_fit,param_fixed,fixed,B,xdata,plotdata,mode)

if nargin<6
    mode='plot';
end
if nargin<5
    plotdata='off';
end
if nargin<4
    Bmax=10;
    NB=100;
    B=[0 logspace(-3,log10(Bmax),NB)];
else
    NB=length(B);
    Bmax=max(B);
end
if nargin<3
    if strcmp(mode,'plot') 
        xdata=[];
    else
        xdata=0:0.1:10;
    end
end
%% -Get parameters
if nargin<2
    D=2; E=0; gx=2; gy=2; gz=2;
else

    %Sort fixed and fitted parameters
    i=1;
    j=1;
    for l=1:length(fixed)
        if fixed(l)==1 %is the l-th parameter fixed?
            p(l)=param_fixed(i); 
            i=i+1;
        else
            p(l)=param_fit(j);
            j=j+1;
        end
    end
    D=p(1);
    E=p(2);
    gx=p(3);
    gy=p(4);
    gz=p(5);
end
if nargin<1
    T=2;
end

%}
%CONSTANTS
k_B=1.38E-23;
mu_B=0.6720345504; %In Kelvin units

%Polar Grid
Nth=31;
Nphi=31;

Theta=linspace(0.001,pi/2,Nth);
Phi=linspace(0,pi/2,Nphi);

[Bg,Thg,Phg]=meshgrid(B,Theta,Phi);

%%
%Define spin-Hamiltonian
%Energy in Kelvin units
Ei=zeros(6,size(B,2),size(Theta,2),size(Phi,2));
for k=1:size(B,2)
for l=1:size(Theta,2)
for m=1:size(Phi,2)
    
    H=[25/4*D+mu_B*gz*5/2*B(k)*cos(Theta(l)) mu_B*sqrt(5)/2*(B(k)*sin(Theta(l))*cos(Phi(m))*gx-1i*B(k)*sin(Theta(l))*sin(Phi(m))*gy) sqrt(10)*E 0 0 0;
             mu_B*sqrt(5)/2*(B(k)*sin(Theta(l))*cos(Phi(m))*gx+1i*B(k)*sin(Theta(l))*sin(Phi(m))*gy) 9/4*D+mu_B*gz*3/2*B(k)*cos(Theta(l)) mu_B*sqrt(2)*(B(k)*sin(Theta(l))*cos(Phi(m))*gx-1i*B(k)*sin(Theta(l))*sin(Phi(m))*gy) 3*sqrt(2)*E 0 0;
             sqrt(10)*E mu_B*sqrt(2)*(B(k)*sin(Theta(l))*cos(Phi(m))*gx+1i*B(k)*sin(Theta(l))*sin(Phi(m))*gy) 1/4*D+mu_B*gz*1/2*B(k)*cos(Theta(l)) mu_B*3/2*(B(k)*sin(Theta(l))*cos(Phi(m))*gx-1i*B(k)*sin(Theta(l))*sin(Phi(m))*gy) 3*sqrt(2)*E 0;
             0 3*sqrt(2)*E mu_B*3/2*(B(k)*sin(Theta(l))*cos(Phi(m))*gx+1i*B(k)*sin(Theta(l))*sin(Phi(m))*gy) 1/4*D-mu_B*gz*1/2*B(k)*cos(Theta(l)) mu_B*sqrt(2)*(B(k)*sin(Theta(l))*cos(Phi(m))*gx-1i*B(k)*sin(Theta(l))*sin(Phi(m))*gy) sqrt(10)*E;
             0 0 3*sqrt(2)*E mu_B*sqrt(2)*(B(k)*sin(Theta(l))*cos(Phi(m))*gx+1i*B(k)*sin(Theta(l))*sin(Phi(m))*gy) 9/4*D-mu_B*gz*3/2*B(k)*cos(Theta(l)) mu_B*sqrt(5)/2*(B(k)*sin(Theta(l))*cos(Phi(m))*gx-1i*B(k)*sin(Theta(l))*sin(Phi(m))*gy);
             0 0 0 sqrt(10)*E mu_B*sqrt(5)/2*(B(k)*sin(Theta(l))*cos(Phi(m))*gx+1i*B(k)*sin(Theta(l))*sin(Phi(m))*gy) 25/4*D-mu_B*gz*5/2*B(k)*cos(Theta(l))];
    

%Diagonalize H to find energy eigenvalues (i=1..6)
Ei(:,k,l,m)=eig(H);      %Ei(i=1..6,B,Theta,Phi)
    
end
end
end
%{
Z=squeeze(sum(exp(-Ei/T),1));     %Z=sum(-Ei/kB*T), Z(B,Theta,Phi)

[dZr,dZt,dZp] = gradient(permute(log(Z),[2,1,3]),B,Theta,Phi); %[dZ/dB,dZ/dTheta,dZ/dPhi](B,Theta,Phi)
dZt=dZt./Bg;
dZp=dZp./(Bg.*sin(Thg));

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
M_x=squeeze(Mx(Nth,1,:)); %Easy-plane magnetization Mx(theta=90°,phi=0°,B)
M_y=squeeze(My(Nth,Nphi,:)); %Easy-plane magnetization Mx(theta=90°,phi=90°,B)

Mz_B=Mx.*(sin(Theta')).*cos(Phi)+My.*(sin(Theta')).*sin(Phi)+Mz.*cos(Theta');
M_av=4/pi^2*squeeze(trapz(Theta,(trapz(Phi,Mz_B,2)),1));

Bsym=sort([-B(2:end) 0 B(2:end)]);
M_av_sym=[-flipud(M_av(2:end));0;M_av(2:end)];
M_av_x=interp1(Bsym,M_av_sym,xdata);
%}

%%T=2K
T=2;
[M_x,M_y,M_z,M_av]=Calc_Magn_3D(Ei,T,B,Theta,Phi);
%Interpolate to original measurements
Bsym=sort([-B(2:end) 0 B(2:end)]);
M_av_sym=[-flipud(M_av(2:end));0;M_av(2:end)];
M_av_x=interp1(Bsym,M_av_sym,xdata);

%%T=5K
T=5;

[M_x_5K,M_y_5K,M_z_5K,M_av_5K]=Calc_Magn_3D(Ei,T,B,Theta,Phi);
%M_av_5K_int=interp1(B,M_av_5K,xdata(sizes(1)+sizes(2)+sizes(3)+1:sizes(1)+sizes(2)+sizes(3)+sizes(4)));

%%
%{
%%T-dependence
Tv=1:1:300;
Mx_T=zeros(size(Tv));
My_T=zeros(size(Tv));
Mz_T=zeros(size(Tv));
Mav_T=zeros(size(Tv));
B0=0.5;
for t=1:300
    [Mx_T(t),My_T(t),Mz_T(t),Mav_T(t)]=Calc_MT_3D(Ei,Tv(t),B,Theta,Phi,B0);
end

Mx_T_int=interp1(Tv,Mx_T,xdata(sizes(1)+sizes(2)+sizes(3)+sizes(4)+1:sizes(1)+sizes(2)+sizes(3)+sizes(4)+sizes(5)));
Mav_T_int=interp1(Tv,Mav_T,xdata(sizes(1)+sizes(2)+sizes(3)+sizes(4)+sizes(5)+1:end));
%}

switch mode
case 'plot'
    OUTPUT=[M_x,M_y,M_z,M_av,M_x_5K,M_y_5K,M_z_5K,M_av_5K];
case 'fit'
    OUTPUT=M_av_x;
end
%%

switch plotdata
case 'on'
%PLOTTING RESULTS


[BXg,BYg,BZg]=sph2cart(Phg,pi/2-Thg,Bg);
BXg=permute(BXg,[2,1,3]);
BYg=permute(BYg,[2,1,3]);
BZg=permute(BZg,[2,1,3]);

dlnZ_dBx=permute(dlnZ_dBx,[2,1,3]);
dlnZ_dBy=permute(dlnZ_dBy,[2,1,3]);
dlnZ_dBz=permute(dlnZ_dBz,[2,1,3]);

Mx=permute(Mx,[3,1,2]);
My=permute(My,[3,1,2]);
Mz=permute(Mz,[3,1,2]);

if strcmp(plotdata,'on')

  figure(1)
    subplot(6,1,1)
    cla
    box on
    hold on
    contour(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(Ei(1,:,:,1)))
    text(0.9*Bmax,0.7*Bmax,'E_1(K)')

    subplot(6,1,2)
    cla
    box on
    hold on
    contour(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(Ei(2,:,:,1)))
    text(0.9*Bmax,0.7*Bmax,'E_2(K)')

    subplot(6,1,3)
    cla
    box on
    hold on
    contour(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(Ei(3,:,:,1)))
    text(0.9*Bmax,0.7*Bmax,'E_3(K)')
    ylabel('B_z(T)')

    subplot(6,1,4)
    cla
    box on
    hold on
    contour(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(Ei(4,:,:,1)))
    text(0.9*Bmax,0.7*Bmax,'E_4(K)')

    subplot(6,1,5)
    cla
    box on
    hold on
    contour(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(Ei(5,:,:,1)))
    text(0.9*Bmax,0.7*Bmax,'E_5(K)')

    subplot(6,1,6)
    cla
    box on
    hold on
    contour(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(Ei(6,:,:,1)))
    text(0.9*Bmax,0.7*Bmax,'E_6(K)')
    xlabel('B_x(T)')
    
    
    figure(11)
    subplot(6,1,1)
    cla
    box on
    hold on
    contour(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(Ei(1,:,:,Nphi)))
    text(0.9*Bmax,0.7*Bmax,'E_1(K)')

    subplot(6,1,2)
    cla
    box on
    hold on
    contour(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(Ei(2,:,:,Nphi)))
    text(0.9*Bmax,0.7*Bmax,'E_2(K)')

    subplot(6,1,3)
    cla
    box on
    hold on
    contour(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(Ei(3,:,:,Nphi)))
    text(0.9*Bmax,0.7*Bmax,'E_3(K)')
    ylabel('B_z(T)')

    subplot(6,1,4)
    cla
    box on
    hold on
    contour(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(Ei(4,:,:,Nphi)))
    text(0.9*Bmax,0.7*Bmax,'E_4(K)')

    subplot(6,1,5)
    cla
    box on
    hold on
    contour(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(Ei(5,:,:,Nphi)))
    text(0.9*Bmax,0.7*Bmax,'E_5(K)')

    subplot(6,1,6)
    cla
    box on
    hold on
    contour(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(Ei(6,:,:,Nphi)))
    text(0.9*Bmax,0.7*Bmax,'E_6(K)')
    xlabel('B_y(T)')

    figure(2)
    subplot(1,3,1)
    cla
    box on
    hold all
    plot(B,squeeze(Ei(1,:,1,1)),'k',...
    B,squeeze(Ei(2,:,1,1)),'r',...
    B,squeeze(Ei(3,:,1,1)),'g',...
    B,squeeze(Ei(4,:,1,1)),'b',...
    B,squeeze(Ei(5,:,1)),'y',...
    B,squeeze(Ei(6,:,1,1)),'m')
    legend('E_1(Bz)(B_x=0,B_y=0)','E_2(Bz)(B_x=0,B_y=0)','E_3(Bz)(B_x=0,B_y=0)','E_4(Bz)(B_x=0,B_y=0)','E_5(Bz)(B_x=0,B_y=0)','E_6(Bz)(B_x=0,B_y=0)')
    xlabel('B_z(T)')
    ylabel('E(K)')

    subplot(1,3,2)
    cla
    box on
    hold all
    plot(B,squeeze(Ei(1,:,Nth,1)),'k',...
    B,squeeze(Ei(2,:,Nth,1)),'r',...
    B,squeeze(Ei(3,:,Nth,1)),'g',...
    B,squeeze(Ei(4,:,Nth,1)),'b',...
    B,squeeze(Ei(5,:,Nth,1)),'y',...
    B,squeeze(Ei(6,:,Nth,1)),'m')

    legend('E_1(Bx)(B_z=0,B_y=0)','E_2(Bx)(B_z=0,B_y=0)','E_3(Bx)(B_z=0,B_y=0)','E_4(Bx)(B_z=0,B_y=0)','E_5(Bx)(B_z=0,B_y=0)','E_6(Bx)(B_z=0,B_y=0)')
    xlabel('B_x(T)')
    ylabel('E(K)')
    
    subplot(1,3,3)
    cla
    box on
    hold all
    plot(B,squeeze(Ei(1,:,Nth,Nphi)),'k',...
    B,squeeze(Ei(2,:,Nth,Nphi)),'r',...
    B,squeeze(Ei(3,:,Nth,Nphi)),'g',...
    B,squeeze(Ei(4,:,Nth,Nphi)),'b',...
    B,squeeze(Ei(5,:,Nth,Nphi)),'y',...
    B,squeeze(Ei(6,:,Nth,Nphi)),'m')

    legend('E_1(By)(B_x=0,B_z=0)','E_2(By)(B_x=0,B_z=0)','E_3(By)(B_x=0,B_z=0)','E_4(By)(B_x=0,B_z=0)','E_5(By)(B_x=0,B_z=0)','E_6(By)(B_x=0,B_z=0)')
    xlabel('B_y(T)')
    ylabel('E(K)')

    figure(3)
    subplot(2,2,[1,2])
    surf(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(log(Z(:,:,1))))
    hold on
    contour(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(log(Z(:,:,1))))
    xlabel('Bx')
    ylabel('Bz')
    zlabel('Partition function, Z| B_y=0')

    subplot(2,2,3)
    surf(squeeze(squeeze(BXg(:,:,1))),squeeze(BZg(:,:,1)),squeeze(dlnZ_dBx(:,:,1)))
    xlabel('Bx')
    ylabel('Bz')
    zlabel('\partial(lnZ)/\partial(Bx) | B_y=0')
    subplot(2,2,4)
    surf(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(dlnZ_dBz(:,:,1)))
    xlabel('By')
    ylabel('Bz')
    zlabel('\partial(lnZ)/\partial(Bz) | B_y=0')
    
    figure(33)
    subplot(2,2,[1,2])
    surf(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(log(Z(:,:,Nphi))))
    hold on
    contour(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(log(Z(:,:,Nphi))))
    xlabel('By')
    ylabel('Bz')
    zlabel('Partition function, Z| B_x=0')

    subplot(2,2,3)
    surf(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(dlnZ_dBy(:,:,Nphi)))
    xlabel('By')
    ylabel('Bz')
    zlabel('\partial(lnZ)/\partial(By) | B_x=0')
    subplot(2,2,4)
    surf(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(dlnZ_dBz(:,:,Nphi)))
    xlabel('Bx')
    ylabel('Bz')
    zlabel('\partial(lnZ)/\partial(Bz) | B_x=0')

    figure(4678)
    subplot(1,2,1)
    cla
    surf(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(Mx(:,:,1)))
    hold on
    contour(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(Mx(:,:,1)))
    xlabel('B_x(T)')
    ylabel('B_z(T)')
    zlabel('M_x(\mu_B) | B_y=0')

    subplot(1,2,2)
    cla
    surf(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(Mz(:,:,1)))
    hold on
    contour(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),squeeze(Mz(:,:,1)))
    xlabel('B_x(T)')
    ylabel('B_z(T)')
    zlabel('M_z(\mu_B)| B_y=0')
    
    figure(44)
    subplot(1,2,1)
    cla
    surf(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(My(:,:,Nphi)))
    hold on
    contour(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(My(:,:,Nphi)))
    xlabel('B_y(T)')
    ylabel('B_z(T)')
    zlabel('M_y(\mu_B) | B_x=0')

    subplot(1,2,2)
    cla
    surf(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(Mz(:,:,Nphi)))
    hold on
    contour(squeeze(BYg(:,:,Nphi)),squeeze(BZg(:,:,Nphi)),squeeze(Mz(:,:,Nphi)))
    xlabel('B_y(T)')
    ylabel('B_z(T)')
    zlabel('M_z(\mu_B)| B_x=0')
    
 figure(6)
    subplot(1,2,1)
    cla
    surf(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),Mx(:,:,1))
    hold on
    contour(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),Mx(:,:,1))
    xlabel('B_x(T)')
    ylabel('B_z(T)')
    zlabel('M_x(\mu_B)')

    subplot(1,2,2)
    cla
    surf(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),Mz(:,:,1))
    hold on
    contour(squeeze(BXg(:,:,1)),squeeze(BZg(:,:,1)),Mz(:,:,1))
    xlabel('B_x(T)')
    ylabel('B_z(T)')
    zlabel('M_x(\mu_B)')   
    
 figure(7)
    cla
    h_ha=plot(B,M_ha,'k');
    hold on
    h_ep_x=plot(B,M_ep_x,'r');
    h_ep_y=plot(B,M_ep_y,'b');
    h_av=plot(B,M_av,'g-');

    xlabel('B(T)')
    ylabel('M(\mu_B)')
    legend([h_ha,h_ep_x,h_ep_y,h_av],'Mz(B)','Mx(B)','My(B)','M_av(B)')
    
end
end