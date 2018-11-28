function M=Isotropic_mf(T,param_fit,param_fixed,fixed,xdata)

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
    J0=p(3);
end


B1=xdata;
M_Br=Brillouin(T,[g,S],[],[0,0],B1);
if J0<=0
    M0=M_Br;
else
    M0=ones(size(M_Br));
end
tic
    fun=@(M)selfcons(T,g,S,B1,J0,M);
    options = optimoptions('fsolve');
    M=fsolve(fun,M0,options);

toc