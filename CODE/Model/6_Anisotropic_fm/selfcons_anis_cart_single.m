function F=selfcons_anis_cart_single(T,param,Bx,Bz,J0,M)

if nargin<1
    T=2;
    gx=2;
    gz=2;
    D=10;
    Bx=5;
    Bz=5;
else
    D=param(1);
    gx=param(2);
    gz=param(3);
    
    Mx=M(1);
    Mz=M(2);
end

muB=0.6720345504; %In [K/T] units (muB/kB) -> J0 is defined in [K]


%% - CALCULATE Effective magnetic fields for all (B,Theta)
Beffx=Bx+2/gx^2*J0/muB*Mx;
Beffz=Bz+2/gz^2*J0/muB*Mz;

M2=Anispm_cart_single(T,[D,gx,gz],[],[0,0,0],Beffx,Beffz); %Mx(B,Theta), Mz(B,Theta); B:-10..10, Theta: 0...2pi ('look-up table' for Mx and Mz in a polar grid)
%Mx2=squeeze(M2(:,:,1));
%Mz2=squeeze(M2(:,:,2));

F=M-M2;

%{
figure(1532)
subplot(1,2,1)
cla
hold on
surf(BXg,BZg,Beffx)
subplot(1,2,2)
cla
hold on
surf(BXg,BZg,Beffz)
%[Bgx,Bgz]=meshgrid(Beffx,Beffz);
%[B]

%Mx2=interp2(Bg,Thg,Mx
%F=M2-M;

figure(4321122)
cla
hold on
plot(B,M2,'k-')
plot(B,M,'r-')
drawnow
%}