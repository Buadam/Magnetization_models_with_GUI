function [Mxx_T,Mzz_T,Mav_T]=Anis_2D_mf_mT(B0,param,T,Theta)

if nargin<4
    Theta=linspace(0,pi/2,61);
end
if nargin<3
    T=linspace(1,300,300);
end
if nargin<2
    D=8.3;
    gx=2;
    gz=2;
    J0=0.1;
else
D=param(1); gx=param(2); gz=param(3); J0=param(4);
end
if nargin<1
    T0=2;
end

NT=length(T); Tmax=max(T);
NTh=length(Theta); Thmax=max(Theta);

Mx_T=zeros(NT,NTh);
Mz_T=zeros(NT,NTh);   

parfor i=1:NTh*NT
     Th=linspace(0,Thmax,NTh);
     Tv=linspace(1,Tmax,NT);
     [k,l]=ind2sub([NT,NTh],i);
     Bx=B0*cos(Th(l));
     Bz=B0*sin(Th(l));
     M0=Anispm_cart_single(Tv(k),[D,gx,gz],[],[0,0,0],Bx,Bz); %Mx(B,Theta), Mz(B,Theta); B:-10..10, Theta: 0...2pi ('look-up table' for Mx and Mz in a polar grid)
     
     %% - Solving self-consistent equation to find magnetization
     fun=@(M)selfcons_anis_cart_single(Tv(k),[D,gx,gz],Bx,Bz,J0,M);
     options = optimoptions('fsolve','Display','off');
     M=fsolve(fun,M0,options);
     
     Mx_T(i)=M(1);
     Mz_T(i)=M(2);
end
            
Mxx_T=Mx_T(:,1);
Mzz_T=Mz_T(:,end);
Mav_T=2/pi*trapz(Theta,(Mx_T.*cos(Theta)+Mz_T.*sin(Theta)),2);  
          
end