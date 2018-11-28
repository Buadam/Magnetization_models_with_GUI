function OUTPUT = Brillouin_mT(B,param,T)

if nargin<3
    T=linspace(1,300,300);
end
if param<2
    g=2; S=5/2;
else
    g=param(1);
    S=param(2);
end
if nargin<1
    B=0.5;
end



muB=0.6720345504;

x=g*muB*S*B./T;
Br=g*S*((2*S+1)/(2*S)*coth((2*S+1)/(2*S)*x)-1/(2*S)*coth(1/(2*S)*x));

OUTPUT=Br;

OUTPUT(find(x==0))=0;
