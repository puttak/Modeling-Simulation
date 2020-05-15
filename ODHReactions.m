function rxn = ODHReactions(C_solid,T,R,Pt,RxnKinetic,deltaS0,deltaH0,compnumber)
% This code is for modeling of ODH reaction kinetics

Ct_solid = sum(C_solid);   % total mole concentration in solid phase

% component order list: [C2H6 C2H4 O2 CO2 CO H2O N2]
P_solid = Pt*(C_solid/Ct_solid);

% component order list: [C2H6 C2H4 O2 CO2 CO H2O]
K = ones(1,6);
for n = 1:6
    K(n) = exp((deltaS0 - deltaH0*(1/T - 1/(25+273.15)) )/R );
end

Tetha_star = 1/(1 + K(1)*P_solid(1) + K(2)*P_solid(2) + sqrt(K(3)*P_solid(3)) + ...
                K(4)*P_solid(4) + K(5)*P_solid(5) + K(6)*P_solid(6) );
            
Tetha_O    = sqrt(K(3)*P_solid(3))/Tetha_star;
Tetha_C2H6 = K(1)*P_solid(1)/Tetha_star;
Tetha_C2H4 = K(2)*P_solid(2)/Tetha_star;

k = ones(1,5);
rxn = k(:);
for i = 1:5
    k(i)   = exp(RxnKinetic.Aprime(i) - RxnKinetic.EnergyA(i)/R *(1/T - 1/(25+273.15) ) );
    if i < 4    
        rxn(i) = k(i) * Tetha_O^(RxnKinetic.m) * Tetha_C2H6;
    else
        rxn(i) = k(i) * Tetha_O^(RxnKinetic.m) * Tetha_C2H4;
    end
end

% sum of reactions
rxn = rxn * RxnKinetic.vcoffrxn(:,compnumber) ;

end
