function Magn=Anisotropic_mf_single(T,param_fit,param_fixed,fixed)
tic
%% -Get parameters
if nargin<4
    fixed=[0,0,0];
end
if nargin<3
    param_fixed=[];
end
if nargin<2
    gx=2; gz=2;
    D=10; J0=0; T=2;
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
    J0=p(2);
end

%% - Define Cartesian grid
Bxmax=100;
Bzmax=100;
Nx=50;
Nz=100;
Bx=linspace(0,Bxmax,Nx);    
Bz=linspace(0,Bzmax,Nz);

[BXg,BZg]=meshgrid(Bx,Bz);
BXg=permute(BXg,[2,1]);
BZg=permute(BZg,[2,1]);

Mx=zeros(Nx,Nz);
Mz=zeros(Nx,Nz);

%% - CALCULATE Mx(B,Theta) and Mz(B,Theta) functions for a full polar grid
for k=1:length(Bx)
for l=1:length(Bz)
    %disp(['[',num2str(k),',',num2str(l)])
    M0=Anispm_cart_single(T,[D,gx,gz],[],[0,0,0],Bx(k),Bz(l)); %Mx(B,Theta), Mz(B,Theta); B:-10..10, Theta: 0...2pi ('look-up table' for Mx and Mz in a polar grid)
    %% - Solving self-consistent equation to find magnetization
    fun=@(M)selfcons_anis_cart_single(T,[D,gx,gz],Bx(k),Bz(l),J0,M);
    options = optimoptions('fsolve');
    M=fsolve(fun,M0,options);
    
    Mx(k,l)=M(1);
    Mz(k,l)=M(2);
end
end

figure(432)
plot(Bx,Mx(:,1))
hold on
plot(Bz,Mz(1,:))

toc
Magn=[Mx;Mz];
