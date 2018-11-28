function [Mxx,Mzz,Mav]=Anis_2D_mH(T0,param,B,Theta)

if nargin<4
    Theta=linspace(0,pi/2,61);
end
if nargin<3
    B=linspace(0,100,100);
end
if nargin<2
    D=8.3;
    gx=2;
    gz=2;
else
D=param(1); gx=param(2); gz=param(3);
end
if nargin<1
    T0=2;
end

NB=length(B); Bmax=max(B);
NTh=length(Theta); Thmax=max(Theta);

Mx=zeros(NB,NTh);
Mz=zeros(NB,NTh);

%% - CALCULATE Mx(B,Theta) and Mz(B,Theta) functions for a full polar grid
parfor i=1:NB*NTh
    Th=linspace(0,Thmax,NTh);
    B=linspace(0,Bmax,NB);   
    [k,l]=ind2sub([NB,NTh],i);
    Bx=B(k)*cos(Th(l));
    Bz=B(k)*sin(Th(l));
    M=Anis_2D_single(T0,[D,gx,gz],[],[0,0,0],Bx,Bz); %Mx(B,Theta), Mz(B,Theta); B:-10..10, Theta: 0...2pi ('look-up table' for Mx and Mz in a polar grid)
    
    Mx(i)=M(1);
    Mz(i)=M(2);
end

Mxx=Mx(:,1);
Mzz=Mz(:,end);
Mav=2/pi*trapz(Theta,(Mx.*cos(Theta)+Mz.*sin(Theta)),2);


figure(423)
cla
plot(B,Mxx,'g-')
hold on
plot(B,Mzz,'r')
plot(B,Mav,'b')

end