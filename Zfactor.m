function Z=Zfactor(Tr,Pr)
if Pr>0.2 && Pr<1.2 && Tr>1.05 && Tr<1.2
    Z=Pr*(1.6643*Tr-2.2114)-0.3647*Tr+1.4385;
elseif Pr>0.2 && Pr<1.2 && Tr>1.2 && Tr<1.4
    Z=Pr*(0.5222*Tr-0.8511)-0.0364*Tr+1.0490;
elseif Pr>0.2 && Pr<1.2 && Tr>1.4 && Tr<2z
    Z=Pr*(0.1391*Tr-0.2988)+0.0007*Tr+0.9969;
elseif Pr>0.2 && Pr<1.2 && Tr>2 && Tr<3
    Z=Pr*(0.0295*Tr-0.0825)+0.0009*Tr+0.9967;
elseif Pr>1.2 && Pr<2.8 && Tr>1.05 && Tr<1.2
    Z=Pr*(-1.3570*Tr+1.4942)+4.6315*Tr-4.7009;
elseif Pr>1.2 && Pr<2.8 && Tr>1.2 && Tr<1.4
    Z=Pr*(0.1717*Tr-0.3232)+0.5869*Tr+0.1229;
elseif Pr>1.2 && Pr<2.8 && Tr>1.4 && Tr<2
    Z=Pr*(0.0984*Tr-0.2053)+0.0621*Tr+0.8580;
elseif Pr>1.2 && Pr<2.8 && Tr>2 && Tr<3
    Z=Pr*(0.0211*Tr-0.0527)+0.0127*Tr+0.9549;
elseif Pr>2.8 && Pr<5.4 && Tr>1.05 && Tr<1.2
    Z=Pr*(-0.3278*Tr+0.4752)+1.8223*Tr-1.9036;
elseif Pr>2.8 && Pr<5.4 && Tr>1.2 && Tr<1.4
    Z=Pr*(-0.2521*Tr+0.3871)+1.6087*Tr-1.6635;
elseif Pr>2.8 && Pr<5.4 && Tr>1.4 && Tr<2
    Z=Pr*(-0.0284*Tr+0.0625)+0.4714*Tr-0.0011;
elseif Pr>2.8 && Pr<5.4 && Tr>2 && Tr<3
    Z=Pr*(0.0041*Tr+0.0039)+0.0607*Tr+0.7927;
elseif Pr>5.4 && Pr<15 && Tr>1.05 && Tr<3
    Z=Pr*(0.711+3.66*Tr)^-1.4667-(1.637/(0.319*Tr+0.522))+2.071;
else
    disp('Error');
    disp('Tr & Pr is:');
    Tr
    Pr
    Z=input('input Z factor from Standing-Katz chart');
end
end
