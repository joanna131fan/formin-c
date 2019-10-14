x=importdata('CAPUsingle');
NFil=2;
k_PAF=10;
c_PA=2.5;
k_poly1=0;
k_poly2=0;
p_occ1=[];
p_occ2=[];
Npp=5;
N_All=5;
kp1=[];
npp=[6 6 9 14 9];
for i=1:Npp
    p_occ1 = [p_occ1 x(16+2*(N_All+1)+7*(i-1))];
end

for n=1:Npp
    k_poly1=k_poly1+((1-p_occ1(n))*npp(n));
    kp1 = [kp1 (1-p_occ1(n))*npp(n)*k_PAF*c_PA];
end
k_poly1=k_poly1*k_PAF*c_PA
kp1 = [kp1; kp1];
fil=[1 2];
figure()
sgtitle('Capu Polymerization Rate Per Filament');
subplot(1,3,1);
bar(fil, kp1,0.5, 'stacked')
title('Single-Filament Capu');
xlim([0.5 1.5]);
xlabel('Filaments');
ylabel('k_p_o_l_y');
ylim([0 180]);

x=importdata('CAPUdouble');
NFil=2;
k_PAF=10;
c_PA=2.5;
k_poly21=0;
k_poly22=0;
p_occ1=[];
p_occ2=[];
Npp=5;
N_All=10;
kp1=[];
kp2=[];
kpboth = [];
npp=[6 6 9 14 9];
for i=1:Npp
    p_occ1 = [p_occ1 x(16+2*(N_All+1)+7*(i-1))];
    p_occ2 = [p_occ2 x(16+2*(N_All+1)+7*(i-1)+(6 + 7*Npp + 2 + NFil + NFil + 2*Npp))];
end

for n=1:Npp
    k_poly21=k_poly21+((1-p_occ1(n))*npp(n));
    kp1 = [kp1 (1-p_occ1(n))*npp(n)*k_PAF*c_PA];
end   
for n=1:Npp
    k_poly22=k_poly22+((1-p_occ2(n))*npp(n));
    kp2 = [kp2 (1-p_occ2(n))*npp(n)*k_PAF*c_PA];
    temp = kp1(n);
    kpboth = [kpboth temp+(1-p_occ2(n))*npp(n)*k_PAF*c_PA];
end    
k_poly21=k_poly21*k_PAF*c_PA
k_poly22=k_poly22*k_PAF*c_PA
kp=[];
kp=[kp1 ; kp2; kpboth];
fil=[1 2 3];
subplot(1,3,2);
bar(fil, kp ,0.5, 'stacked')
title('Double-Filament Capu');
xlabel('Filaments');
ylabel('k_{poly}');
xticklabels({1, 2, 'total'});
ylim([0 180]);

x=importdata('CAPUdimerized');
NFil=2;
k_PAF=10;
c_PA=2.5;
k_poly31=0;
k_poly32=0;
p_occ1=[];
p_occ2=[];
Npp=5;
N_All=10;
kp1=[];
kp2=[];
kpboth=[];
npp=[6 6 9 14 9];
for i=1:Npp
    p_occ1 = [p_occ1 x(16+2*(N_All+1)+7*(i-1))];
    p_occ2 = [p_occ2 x(16+2*(N_All+1)+7*(i-1)+(6 + 7*Npp + 2 + NFil + NFil + 2*Npp))];
end

for n=1:Npp
    k_poly31=k_poly31+((1-p_occ1(n))*npp(n));
    kp1 = [kp1 (1-p_occ1(n))*npp(n)*k_PAF*c_PA];
end
for n=1:Npp
    k_poly32=k_poly32+((1-p_occ2(n))*npp(n));
    kp2 = [kp2 (1-p_occ2(n))*npp(n)*k_PAF*c_PA];
    temp = kp1(n);
    kpboth = [kpboth temp+(1-p_occ2(n))*npp(n)*k_PAF*c_PA];
end    
k_poly31=k_poly31*k_PAF*c_PA
k_poly32=k_poly32*k_PAF*c_PA
kp=[kp1 ; kp2; kpboth];
fil=[1, 2, 3];
subplot(1,3,3);
bar(fil, kp ,0.5, 'stacked')
title('N-Dimerized Capu');
xlabel('Filaments');
ylabel('k_{poly}');
xticklabels({1, 2, 'total'});
ylim([0 180]);

saveas(gcf, 'Capu k-polymerization.pdf');