function M=Isotropic_mf_MH(T0,param,B)

%% -Get parameters
if nargin<3
    B=linspace(0,5,100);
end
if nargin<2
    g=2; S=5/2; J0=0.1;
else
    g=param(1);
    S=param(2);
    J0=param(3);
end
if nargin<1
    T0=2;
end

M_Br=Brillouin(T0,[g,S],[],[0,0],B);
if J0<=0
    M0=M_Br;
else
    M0=ones(size(M_Br));
end

for b=1:length(B)
    Bakt=B(b);
    fun=@(M)selfcons(T0,g,S,Bakt,J0,M);
    options = optimoptions('fsolve');
    M(b)=fsolve(fun,M0(b),options);
end
M=M';
assignin('base','mH',M)