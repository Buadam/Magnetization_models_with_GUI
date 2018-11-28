function OUTPUT = Anis_2D_single(T,param_fit,param_fixed,fixed,Bx,Bz)

if nargin<6
    Bx=5;
end
if nargin<5
    Bz=5;
end
if nargin<4
    fixed=[0,0,0];
end
if nargin<3
    param_fixed=[];
end
if nargin<2
    D=10; gx=2; gz=2;
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
    gx=p(2);
    gz=p(3);
end
if nargin<1
    T=2;
end

dBx=0.001;
dBz=0.001;

%CONSTANTS
k_B=1.38E-23;
mu_B=0.6720345504; %In Kelvin/T units

%%
%Define spin-Hamiltonian
%Energy in Kelvin units

    Hx2=[25/4*D+mu_B*gz*5/2*Bz mu_B*gx*sqrt(5)/2*(Bx+dBx) 0 0 0 0;
             mu_B*gx*sqrt(5)/2.*(Bx+dBx) 9/4*D+mu_B*gz*3/2*Bz mu_B*gx*sqrt(2)*(Bx+dBx) 0 0 0;
             0 mu_B*gx*sqrt(2)*(Bx+dBx) 1/4*D+mu_B*gz*1/2*Bz mu_B*gx*3/2*(Bx+dBx) 0 0;
             0 0 mu_B*gx*3/2.*(Bx+dBx) 1/4*D-mu_B*gz*1/2*Bz mu_B*gx*sqrt(2)*(Bx+dBx) 0;
             0 0 0 mu_B*gx*sqrt(2)*(Bx+dBx) 9/4*D-mu_B*gz*3/2*Bz mu_B*gx*sqrt(5)/2*(Bx+dBx);
             0 0 0 0 mu_B*gx*sqrt(5)/2*(Bx+dBx) 25/4*D-mu_B*gz*5/2*Bz]; 
         
    Hx1=[25/4*D+mu_B*gz*5/2*Bz mu_B*gx*sqrt(5)/2*(Bx-dBx) 0 0 0 0;
             mu_B*gx*sqrt(5)/2.*(Bx-dBx) 9/4*D+mu_B*gz*3/2*Bz mu_B*gx*sqrt(2)*(Bx-dBx) 0 0 0;
             0 mu_B*gx*sqrt(2)*(Bx-dBx) 1/4*D+mu_B*gz*1/2*Bz mu_B*gx*3/2*(Bx-dBx) 0 0;
             0 0 mu_B*gx*3/2*(Bx-dBx) 1/4*D-mu_B*gz*1/2*Bz mu_B*gx*sqrt(2)*(Bx-dBx) 0;
             0 0 0 mu_B*gx*sqrt(2)*(Bx-dBx) 9/4*D-mu_B*gz*3/2*Bz mu_B*gx*sqrt(5)/2*(Bx-dBx);
             0 0 0 0 mu_B*gx*sqrt(5)/2*(Bx-dBx) 25/4*D-mu_B*gz*5/2*Bz]; 
         
    Hz2=[25/4*D+mu_B*gz*5/2*(Bz+dBz) mu_B*gx*sqrt(5)/2*Bx 0 0 0 0;
             mu_B*gx*sqrt(5)/2.*Bx 9/4*D+mu_B*gz*3/2*(Bz+dBz) mu_B*gx*sqrt(2)*Bx 0 0 0;
             0 mu_B*gx*sqrt(2)*Bx 1/4*D+mu_B*gz*1/2*(Bz+dBz) mu_B*gx*3/2*Bx 0 0;
             0 0 mu_B*gx*3/2*Bx 1/4*D-mu_B*gz*1/2*(Bz+dBz) mu_B*gx*sqrt(2)*Bx 0;
             0 0 0 mu_B*gx*sqrt(2)*Bx 9/4*D-mu_B*gz*3/2*(Bz+dBz) mu_B*gx*sqrt(5)/2*Bx;
             0 0 0 0 mu_B*gx*sqrt(5)/2*Bx 25/4*D-mu_B*gz*5/2*(Bz+dBz)]; 
         
    Hz1=[25/4*D+mu_B*gz*5/2*(Bz-dBz) mu_B*gx*sqrt(5)/2*Bx 0 0 0 0;
             mu_B*gx*sqrt(5)/2.*Bx 9/4*D+mu_B*gz*3/2*(Bz-dBz) mu_B*gx*sqrt(2)*Bx 0 0 0;
             0 mu_B*gx*sqrt(2)*Bx 1/4*D+mu_B*gz*1/2*(Bz-dBz) mu_B*gx*3/2*Bx 0 0;
             0 0 mu_B*gx*3/2*Bx 1/4*D-mu_B*gz*1/2*(Bz-dBz) mu_B*gx*sqrt(2)*Bx 0;
             0 0 0 mu_B*gx*sqrt(2)*Bx 9/4*D-mu_B*gz*3/2*(Bz-dBz) mu_B*gx*sqrt(5)/2*Bx;
             0 0 0 0 mu_B*gx*sqrt(5)/2*Bx 25/4*D-mu_B*gz*5/2*(Bz-dBz)]; 


%Diagonalize H to find energy eigenvalues (i=1..6) 
    Ex2=eig(Hx2);      %Ei(i=1..6)
    Ex1=eig(Hx1);
    Ez2=eig(Hz2);
    Ez1=eig(Hz1);


Zx2=squeeze(sum(exp(-Ex2/T),1));     
Zx1=squeeze(sum(exp(-Ex1/T),1));     
Zz2=squeeze(sum(exp(-Ez2/T),1));     
Zz1=squeeze(sum(exp(-Ez1/T),1));     

dlnZ_dBx=(log(Zx2)-log(Zx1))/dBx/2;
dlnZ_dBz=(log(Zz2)-log(Zz1))/dBx/2;

%%
%Magnetization in the polar reference frame in mu_B units 
Mx=k_B*T*dlnZ_dBx/9.27e-24; %Mx(B,Theta)
Mz=k_B*T*dlnZ_dBz/9.27e-24; %Mz(B,Theta)



OUTPUT=[Mx;Mz];


%%

%M_z=Mz(:,1); %Hard axis magnetization Mz(B,theta=0)
%M_x=Mx(:,Nth); %Easy-plane magnetization Mx(B,theta=90°)
%M_av=2/pi*(trapz(Theta,(Mx.*sin(Theta)+Mz.*cos(Theta)),2)); %averaging to the polar angle


%Bsym=sort([-B(2:end) 0 B(2:end)]);
%M_av_sym=[-flipud(M_av(2:end));0;M_av(2:end)];
%M_av_int=interp1(Bsym,M_av_sym,xdata);


