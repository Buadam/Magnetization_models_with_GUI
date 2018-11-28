function M=Anisotropic_mf(T,param_fit,param_fixed,fixed,xdata)

%% -Get parameters
if nargin<2
    D=2; J0=0;
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

g=2;
S=5/2;

%% - Initial guess by Isotropic Brillouin function
B1=xdata;
M_Br=Brillouin(T,[g,S],[],[0,0],B1);
M0=M_Br;

%% - Solving self-consistent equation to find magnetization
fun=@(Mx,Mz)selfcons_anis(T,g,D,B1,J0,Mx,Mz);
options = optimoptions('fsolve');
M=fsolve(fun,M0,options);
