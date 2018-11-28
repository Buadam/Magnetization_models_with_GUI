function M=Isotropic_mf_MT(B0,param,T)

%% -Get parameters
if nargin<3
    T=linspace(1,300,300);
end
if nargin<2
    g=2; S=5/2; J0=0.1;
else
    g=param(1);
    S=param(2);
    J0=param(3);
end
if nargin<1
    B0=0.5;
end

%muB=0.6720345504;
%TC=2*J0*S*(S+1)/3
%M=g^2*muB*S*(S+1)./(T-TC)/3*B0;


M_Br_T=Brillouin_mT(B0,[g,S],T);
if J0<=0
    M0=0.9*M_Br_T;
else
    M0=0.1*ones(size(M_Br_T));
end


fun=@(M)selfcons(T,g,S,B0,J0,M);
options = optimoptions('fsolve');
M=fsolve(fun,M0,options);

assignin('base','mT',M)
%{
M=zeros(length(T));
for t=1:length(T)
    Takt=T(t);
    fun=@(M)selfcons(Takt,g,S,B0,J0,M);
    options = optimoptions('fsolve');
    M(t)=fsolve(fun,M0(t),options);
end
%}
end