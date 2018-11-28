function OUTPUT=fit_Brillouin(T0,B0,param_fit,param_fixed,fixed,xdata,sizes,mode)

if nargin<8
    mode=[1,0,0,0];
end
if nargin<7
    sizes=[];
end
if nargin<6
    xdata=linspace(0,5,100);
end
if nargin<5
    fixed=[0,0];
end
if nargin<4
    param_fixed=[];
end
%% -Get parameters
if nargin<3
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
if nargin<2
    B0=0.5;
end
if nargin<1
    T0=2;
end

OUTPUT=[];
if mode(1)==1
    OUTPUT=Brillouin_mH(T0,[g,S],xdata(sizes(1,1):sizes(1,2)));
end
if mode(2)==1
    OUTPUT2=Brillouin_mH(T0,[g,S],xdata(sizes(2,1):sizes(2,2)));
    OUTPUT=[OUTPUT;OUTPUT2];
end
if mode(3)==1
    OUTPUT3=Brillouin_mT(B0,[g,S],xdata(sizes(3,1):sizes(3,2)));
    OUTPUT=[OUTPUT;OUTPUT3];
end
if mode(4)==1
    OUTPUT4=Brillouin_mT(B0,[g,S],xdata(sizes(4,1):sizes(4,2)));
    OUTPUT=[OUTPUT;OUTPUT4];
end