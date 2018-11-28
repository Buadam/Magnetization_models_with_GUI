function OUTPUT=Anisotropic_mf_single_pol(T,param_fit,param_fixed,fixed,B,xdata,sizes,mode)
tic
%% -Get parameters
if nargin<8
    mode='plot';
end
if nargin<7
   sizes=[];
end
if nargin<6
    xdata=[];
end
if nargin<5
    Bmax=100;
    NB=100;
    B=linspace(0.01,Bmax,NB);    
else
    NB=length(B);
    Bmax=max(B);
end

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
    D=p(1)
    gx=p(2)
    gz=p(3)
    J0=p(4)
end

%% - Define Polar grid
Thmax=pi/2;
NTh=61;

%[BXg,BZg]=meshgrid(Bx,Bz);
%BXg=permute(BXg,[2,1]);
%BZg=permute(BZg,[2,1]);
%[Bg,Thg]=meshgrid(B,Th);
%Bg=permute(Bg,[2,1]);
%Thg=permute(Thg,[2,1]);

Mx=zeros(NB,NTh);
Mz=zeros(NB,NTh);
Th=linspace(0,Thmax,NTh);
B=linspace(0.01,Bmax,NB); 

%% - CALCULATE Mx(B,Theta) and Mz(B,Theta) functions for a full polar grid
parfor i=1:NB*NTh
    Th=linspace(0,Thmax,NTh);
    B=linspace(0,Bmax,NB);   
    %k=mod(i-1,length(B))+1;
    %l=floor((i-1)/length(B))+1;
    [k,l]=ind2sub([NB,NTh],i);
    Bx=B(k)*cos(Th(l));
    Bz=B(k)*sin(Th(l));
    %disp(['[',num2str(k),',',num2str(l)])
    M0=Anispm_cart_single(T,[D,gx,gz],[],[0,0,0],Bx,Bz); %Mx(B,Theta), Mz(B,Theta); B:-10..10, Theta: 0...2pi ('look-up table' for Mx and Mz in a polar grid)
    
    %% - Solving self-consistent equation to find magnetization
    fun=@(M)selfcons_anis_cart_single(T,[D,gx,gz],Bx,Bz,J0,M);
    options = optimoptions('fsolve','Display','off');
    M=fsolve(fun,M0,options);
    
    Mx(i)=M(1);
    Mz(i)=M(2);

end



%Mx=reshape(Mx',length(B),length(Th));
%Mx=reshape(Mz',length(B),length(Th));

Mxx=Mx(:,1);
Mzz=Mz(:,end);
Mav=2/pi*trapz(Th,(Mx.*cos(Th)+Mz.*sin(Th)),2);
%{
figure(432)
cla
hold all
box on
ep=plot(B,Mxx,'g-');
ha=plot(B,Mzz,'r-');
av=plot(B,Mav,'k-');
legend([ep,ha,av],'Mx(Bx)','Mz(Bz)','Mav(B)')
%}

switch mode
    case 'plot'
        OUTPUT=[Mxx,Mzz,Mav];
    case 'plot_Tdep'
        NT=100;
        Tmax=300;
        Mx_T=zeros(NT,NTh);
        Mz_T=zeros(NT,NTh);
            parfor i=1:NTh*NT
                Th=linspace(0,Thmax,NTh);
                Tv=linspace(1,Tmax,NT);
                [k,l]=ind2sub([NT,NTh],i);
                Bx=0.5*cos(Th(l));
                Bz=0.5*sin(Th(l));
                M0=Anispm_cart_single(Tv(k),[D,gx,gz],[],[0,0,0],Bx,Bz); %Mx(B,Theta), Mz(B,Theta); B:-10..10, Theta: 0...2pi ('look-up table' for Mx and Mz in a polar grid)
                %% - Solving self-consistent equation to find magnetization
                fun=@(M)selfcons_anis_cart_single(Tv(k),[D,gx,gz],Bx,Bz,J0,M);
                options = optimoptions('fsolve','Display','off');
                M=fsolve(fun,M0,options);
                Mx_T(i)=M(1);
                Mz_T(i)=M(2);
            end
            
          Mxx_T=Mx_T(:,1);
          Mzz_T=Mz_T(:,end);
          Mav_T=2/pi*trapz(Th,(Mx_T.*cos(Th)+Mz_T.*sin(Th)),2);  
          
          OUTPUT={Mxx,Mzz,Mav,Mxx_T,Mzz_T,Mav_T};
          assignin('base','OUTPUT',OUTPUT)
     case 'fit_av'
        Mav_int=interp1(B,Mav,xdata);
        OUTPUT=Mav_int;
     case 'fit_av_ea'
        Mav_int=interp1(B,Mav,xdata(1:sizes(1)));
        Mxx_int=interp1(B,Mxx,xdata(sizes(1)+1:end));
        Mzz_int=interp1(B,Mzz,xdata(sizes(1)+1:end));
        OUTPUT=[Mav_int;Mzz_int]; 
        
     case 'fit_av_ea_MT'  
        NT=100;
        Tmax=300;
        Tv=linspace(1,Tmax,NT);
        Mx_T=zeros(NT,NTh);
        Mz_T=zeros(NT,NTh);
            parfor i=1:NTh*NT
                Th=linspace(0,Thmax,NTh);
                Tv=linspace(1,Tmax,NT);
                [k,l]=ind2sub([NT,NTh],i);
                Bx=0.5*cos(Th(l));
                Bz=0.5*sin(Th(l));
                M0=Anispm_cart_single(Tv(k),[D,gx,gz],[],[0,0,0],Bx,Bz); %Mx(B,Theta), Mz(B,Theta); B:-10..10, Theta: 0...2pi ('look-up table' for Mx and Mz in a polar grid)
                %% - Solving self-consistent equation to find magnetization
                fun=@(M)selfcons_anis_cart_single(Tv(k),[D,gx,gz],Bx,Bz,J0,M);
                options = optimoptions('fsolve','Display','off');
                M=fsolve(fun,M0,options);
                Mx_T(i)=M(1);
                Mz_T(i)=M(2);
            end
            
          Mxx_T=Mx_T(:,1);
          Mzz_T=Mz_T(:,end);
          Mav_T=2/pi*trapz(Th,(Mx_T.*cos(Th)+Mz_T.*sin(Th)),2);  
          
          
          Mav_int=interp1(B,Mav,xdata(1:sizes(1)));
          Mxx_int=interp1(B,Mxx,xdata(sizes(1)+1:sizes(1)+sizes(2)));
          Mav_T_int=interp1(Tv,Mav_T,xdata(sizes(1)+sizes(2)+1:sizes(1)+sizes(2)+sizes(3)));
          Mxx_T_int=interp1(Tv,Mxx_T,xdata(sizes(1)+sizes(2)+sizes(3)+1:end));
          
          OUTPUT=[Mav_int;Mxx_int;Mav_T_int;Mxx_T_int];
          assignin('base','OUTPUT',OUTPUT)
end
toc
end
  

