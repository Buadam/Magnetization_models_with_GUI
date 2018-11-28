function OUTPUT = Brillouin_mH(T,param,B)

if nargin<3
    Bmax=10;
    NB=100;
    B=[0 logspace(-3,log10(Bmax),NB)];
end
%% -Get parameters
if nargin<2
    g=2; S=5/2;
else
    g=param(1);
    S=param(2);
end
if nargin<1
    T=2;
end



muB=0.6720345504;

x=g*muB*S/T*B;
Br=g*S*((2*S+1)/(2*S)*coth((2*S+1)/(2*S)*x)-1/(2*S)*coth(1/(2*S)*x));

OUTPUT=Br;

OUTPUT(find(x==0))=0;
