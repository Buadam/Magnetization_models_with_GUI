function OUTPUT = Brillouin(T,param_fit,param_fixed,fixed,xdata)

if nargin<5
    Bmax=10;
    NB=100;
    xdata=[0 logspace(-3,log10(Bmax),NB)];
end
if nargin<4
    fixed=[0,0];
end
if nargin<3
    param_fixed=[];
end
%% -Get parameters
if nargin<2
    g=2; S=5/2;
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
    g=p(1);
    S=p(2);
end

if nargin<1
    T=2;
end



muB=0.6720345504;

B=xdata;
x=g*muB*S./T.*B;
Br=g*S*((2*S+1)/(2*S)*coth((2*S+1)/(2*S)*x)-1/(2*S)*coth(1/(2*S)*x));

OUTPUT=Br;

OUTPUT(find(x==0))=0;
