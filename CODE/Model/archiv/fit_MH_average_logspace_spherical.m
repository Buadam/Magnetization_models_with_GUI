function OUTPUT = fit_MH_average_logspace_spherical(T,param,xdata,sizes,plotdata,B)

if nargin<6
    Bmax=15;
    NB=100;
    B=[0 logspace(-3,log10(Bmax),NB)];
end
if nargin<5
    plotdata='off';
end
if nargin<3
    x=0:0.1:10;
    xdata=[x';x';x'];
    sizes=size(x,2)*[1,1,1];
end
if nargin<2
    D=5.8*1.439;
    E=0;
else
    D=param(1);
    E=param(2);
end
if nargin<1
    T=2;
end


%CONSTANTS
k_B=1.38E-23;
gx=2;gz=2; %no g-tensor anisotropy assumed
mu_B=0.6720345504; %In Kelvin units

%Polar Grid

Nth=101;
Nphi=61;
Theta=linspace(0,pi/2,Nth);
Phi=linspace(0,2*pi,Nphi);


%%
%Define spin-Hamiltonian
%Energy in Kelvin units
H=zeros(6,6,size(B,2),size(Theta,2));

for k=1:size(B,2)
for l=1:size(Theta,2)
%for m=1:size(Phi,2)
    H(:,:,k,l)=[25/4*D+mu_B*gz*5/2*B(k)*cos(Theta(l)) mu_B*gx*sqrt(5)/2.*B(k)*sin(Theta(l)) sqrt(10)*E 0 0 0;
             mu_B*gx*sqrt(5)/2.*B(k)*sin(Theta(l)) 9/4*D+mu_B*gz*3/2*B(k)*cos(Theta(l)) mu_B*gx*sqrt(2)*B(k)*sin(Theta(l)) 3*sqrt(2)*E 0 0;
             sqrt(10)*E mu_B*gx*sqrt(2)*B(k)*sin(Theta(l)) 1/4*D+mu_B*gz*1/2*B(k)*cos(Theta(l)) mu_B*gx*3/2*B(k)*sin(Theta(l)) 3*sqrt(2)*E 0;
             0 3*sqrt(2)*E mu_B*gx*3/2.*B(k)*sin(Theta(l)) 1/4*D-mu_B*gz*1/2*B(k)*cos(Theta(l)) mu_B*gx*sqrt(2)*B(k)*sin(Theta(l)) sqrt(10)*E;
             0 0 3*sqrt(2)*E mu_B*gx*sqrt(2)*B(k)*sin(Theta(l)) 9/4*D-mu_B*gz*3/2*B(k)*cos(Theta(l)) mu_B*gx*sqrt(5)/2*B(k)*sin(Theta(l));
             0 0 0 sqrt(10)*E mu_B*gx*sqrt(5)/2*B(k)*sin(Theta(l)) 25/4*D-mu_B*gz*5/2*B(k)*cos(Theta(l))];   
%end
end
end
%Diagonalize H to find energy eigenvalues (i=1..6)
for k=1:size(B,2)
for l=1:size(Theta,2)
    
    Ei(:,k,l)=eig(H(:,:,k,l));      %Ei(i=1..6,Bx,Bz)
    
end
end

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

M_ha=Mz(:,1); %Hard axis magnetization Mz(B,theta=0)
M_ep=Mx(:,Nth); %Easy-plane magnetization Mx(B,theta=90°)
M_av=2/pi*(trapz(Theta,(Mx.*sin(Theta)+Mz.*cos(Theta)),2)); %averaging to the polar angle

M_ep_x=interp1(B,M_ep,xdata(1:sizes(1)));
M_av_x=interp1(B,M_av,xdata(sizes(1)+1:sizes(1)+sizes(2)));
M_av_x2=interp1(B,M_av,xdata(sizes(1)+sizes(2)+1:end));


OUTPUT=[M_ep_x;M_av_x;M_av_x2];


%%

%PLOTTING RESULTS

[Bg,Thg]=meshgrid(B,Theta);
[BZg,BXg]=pol2cart(Thg,Bg);
BXg=permute(BXg,[2,1]);
BZg=permute(BZg,[2,1]);

if strcmp(plotdata,'on')
    
    figure(1)
    subplot(6,1,1)
    cla
    box on
    hold on
    contour(BXg,BZg,squeeze(Ei(1,:,:)))
    text(0.9*Bmax,0.7*Bmax,'E_1(K)')

    subplot(6,1,2)
    cla
    box on
    %surf(BXg,BZg,squeeze(Ei(2,:,:)))
    hold on
    contour(BXg,BZg,squeeze(Ei(2,:,:)))
    text(0.9*Bmax,0.7*Bmax,'E_2(K)')

    subplot(6,1,3)
    cla
    box on
    %surf(BXg,BZg,squeeze(Ei(3,:,:)))
    hold on
    contour(BXg,BZg,squeeze(Ei(3,:,:)))
    text(0.9*Bmax,0.7*Bmax,'E_3(K)')
    ylabel('B_z(T)')

    subplot(6,1,4)
    cla
    box on
    %surf(BXg,BZg,squeeze(Ei(4,:,:)))
    hold on
    contour(BXg,BZg,squeeze(Ei(4,:,:)))
    text(0.9*Bmax,0.7*Bmax,'E_4(K)')

    subplot(6,1,5)
    cla
    box on
    %surf(BXg,BZg,squeeze(Ei(5,:,:)))
    hold on
    contour(BXg,BZg,squeeze(Ei(5,:,:)))
    text(0.9*Bmax,0.7*Bmax,'E_5(K)')

    subplot(6,1,6)
    cla
    box on
    %surf(BXg,BZg,squeeze(Ei(6,:,:)))
    hold on
    contour(BXg,BZg,squeeze(Ei(6,:,:)))
    text(0.9*Bmax,0.7*Bmax,'E_6(K)')
    xlabel('B_x(T)')

    figure(2)
    subplot(1,2,1)
    cla
    box on
    hold all
    plot(B,squeeze(Ei(1,:,1)),'k',...
    B,squeeze(Ei(2,:,1)),'r',...
    B,squeeze(Ei(3,:,1)),'g',...
    B,squeeze(Ei(4,:,1)),'b',...
    B,squeeze(Ei(5,:,1)),'y',...
    B,squeeze(Ei(6,:,1)),'m')
    legend('E_1(B_x=0,B_z)','E_2(B_x=0,B_z)','E_3(B_x=0,B_z)','E_4(B_x=0,B_z)','E_5(B_x=0,B_z)','E_6(B_x=0,B_z)')
    xlabel('B_z(T)')
    ylabel('E(K)')

    subplot(1,2,2)
    cla
    box on
    hold all
    plot(B,squeeze(Ei(1,:,round(Nth/2))),'k',...
    B,squeeze(Ei(2,:,round(Nth/2))),'r',...
    B,squeeze(Ei(3,:,round(Nth/2))),'g',...
    B,squeeze(Ei(4,:,round(Nth/2))),'b',...
    B,squeeze(Ei(5,:,round(Nth/2))),'y',...
    B,squeeze(Ei(6,:,round(Nth/2))),'m')

    legend('E_1(B_x,B_z=0)','E_2(B_x,B_z=0)','E_3(B_x,B_z=0)','E_4(B_x,B_z=0)','E_5(B_x,B_z=0)','E_6(B_x,B_z=0)')
    xlabel('B_x(T)')
    ylabel('E(K)')

    figure(3)
    subplot(2,2,[1,2])
    surf(BXg,BZg,log(Z))
    hold on
    contour(BXg,BZg,log(Z))
    xlabel('Bx')
    ylabel('Bz')
    zlabel('Partition function, Z')

    subplot(2,2,3)
    surf(BXg,BZg,dlnZ_dBx)
    xlabel('Bx')
    ylabel('Bz')
    zlabel('\partial(lnZ)/\partial(Bx)')
    subplot(2,2,4)
    surf(BXg,BZg,dlnZ_dBz)
    xlabel('Bx')
    ylabel('Bz')
    zlabel('\partial(lnZ)/\partial(Bz)')

    figure(4)
    subplot(1,2,1)
    cla
    surf(BXg,BZg,Mx)
    hold on
    contour(BXg,BZg,Mx)
    xlabel('B_x(T)')
    ylabel('B_z(T)')
    zlabel('M_x(\mu_B)')

    subplot(1,2,2)
    cla
    surf(BXg,BZg,Mz)
    hold on
    contour(BXg,BZg,Mz)
    xlabel('B_x(T)')
    ylabel('B_z(T)')
    zlabel('M_z(\mu_B)')

    figure(7)
    cla
    h_ha=plot(B,M_ha,'k');
    hold on
    h_ep=plot(B,M_ep,'r');
    h_av=plot(B,M_av,'b');

    xlabel('B(T)')
    ylabel('M(\mu_B)')
    legend([h_ep,h_ha,h_av],'Model, easy-plane', 'Model, hard-axis', 'Model, average')

end
