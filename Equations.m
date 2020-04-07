function F=Equations(x)
C_C2H6=reshape(x(1:Nz*Nr),Nz,Nr);
C_C2H4=reshape(x(Nz*Nr+1:2*Nz*Nr),Nz,Nr);
C_O2=reshape(x(2*Nz*Nr+1:3*Nz*Nr),Nz,Nr);
C_CO2=reshape(x(3*Nz*Nr+1:4*Nz*Nr),Nz,Nr);
C_CO=reshape(x(4*Nz*Nr+1:5*Nz*Nr),Nz,Nr);
C_H2O=reshape(x(5*Nz*Nr+1:6*Nz*Nr),Nz,Nr);
Cs_C2H6=reshape(x(6*Nz*Nr+1:7*Nz*Nr),Nz,Nr);
Cs_C2H4=reshape(x(7*Nz*Nr+1:8*Nz*Nr),Nz,Nr);
Cs_O2=reshape(x(8*Nz*Nr+1:9*Nz*Nr),Nz,Nr);
Cs_CO2=reshape(x(9*Nz*Nr+1:10*Nz*Nr),Nz,Nr);
Cs_CO=reshape(x(10*Nz*Nr+1:11*Nz*Nr),Nz,Nr);
Cs_H2O=reshape(x(11*Nz*Nr+1:12*Nz*Nr),Nz,Nr);
T=reshape(x(12*Nz*Nr+1:13*Nz*Nr),Nz,Nr);
Ts=reshape(x(13*Nz*Nr+1:14*Nz*Nr),Nz,Nr);

E_C_C2H6=zeros(Nz,Nr);
E_C_C2H4=zeros(Nz,Nr);
E_C_O2=zeros(Nz,Nr);
E_C_CO2=zeros(Nz,Nr);
E_C_CO=zeros(Nz,Nr);
E_C_H2O=zeros(Nz,Nr);
E_Cs_C2H6=zeros(Nz,Nr);
E_Cs_C2H4=zeros(Nz,Nr);
E_Cs_O2=zeros(Nz,Nr);
E_Cs_CO2=zeros(Nz,Nr);
E_Cs_CO=zeros(Nz,Nr);
E_Cs_H2O=zeros(Nz,Nr);
E_T=zeros(Nz,Nr);
E_Ts=zeros(Nz,Nr);

for i=1:Nz
    for k=1:Nr 
        E_C_C2H6(i,k)=
        
        E_C_C2H4(i,k)=
        
        E_C_O2(i,k)=
        
        E_C_CO2(i,k)=
        
        E_C_CO(i,k)=
        
        E_C_H2O(i,k)=
        
        E_Cs_C2H6(i,k)=
        
        E_Cs_C2H4(i,k)=
        
        E_Cs_O2(i,k)=
        
        E_Cs_CO2(i,k)=
        
        E_Cs_CO(i,k)=
        
        E_Cs_H2O(i,k)=
        
        E_T(i,k)=
        
        E_Ts(i,k)=
    end
end

F=[reshape(E_C_C2H6,Nz*Nr,1)  ;  reshape(E_C_C2H4,Nz*Nr,1)  ; ...
   reshape(E_C_O2,Nz*Nr,1)    ;  reshape(E_C_CO2,Nz*Nr,1)   ; ...
   reshape(E_C_CO,Nz*Nr,1)    ;  reshape(E_C_H2O,Nz*Nr,1)   ; ...
   reshape(E_Cs_C2H6,Nz*Nr,1) ;  reshape(E_Cs_C2H4,Nz*Nr,1) ; ...
   reshape(E_Cs_O2,Nz*Nr,1)   ;  reshape(E_Cs_CO2,Nz*Nr,1)  ; ...
   reshape(E_Cs_CO,Nz*Nr,1)   ;  reshape(E_Cs_H2O,Nz*Nr,1)  ; ...
   reshape(E_T,Nz*Nr,1)       ;  reshape(E_Ts,Nz*Nr,1)           ];
   
















