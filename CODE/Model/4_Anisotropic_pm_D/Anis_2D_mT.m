function [Mxx_T,Mzz_T,Mav_T]=Anis_2D_mT(B0,param,T,Theta)

D=param(1); gx=param(2); gz=param(3);
NT=length(T); Tmax=max(T);
NTh=length(Theta); Thmax=max(Theta);

Mx_T=zeros(NT,NTh);
Mz_T=zeros(NT,NTh);   

parfor i=1:NTh*NT
     Th=linspace(0,Thmax,NTh);
     Tv=linspace(1,Tmax,NT);
     [k,l]=ind2sub([NT,NTh],i);
     Bx=B0*cos(Th(l));
     Bz=B0*sin(Th(l));
     M=Anis_2D_single(Tv(k),[D,gx,gz],[],[0,0,0],Bx,Bz); %Mx(B,Theta), Mz(B,Theta); B:-10..10, Theta: 0...2pi ('look-up table' for Mx and Mz in a polar grid)
                
     Mx_T(i)=M(1);
     Mz_T(i)=M(2);
end
            
Mxx_T=Mx_T(:,1);
Mzz_T=Mz_T(:,end);
Mav_T=2/pi*trapz(Theta,(Mx_T.*cos(Theta)+Mz_T.*sin(Theta)),2);  
          
end