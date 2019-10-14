x=importdata('FHODsingle');
%Constants
k_PAF=10;
c_PA=2.5;
Npp=4;
N_All=4;
%Variables
k_poly1=0;
p_occ1=[];
kp1=[];
npp=[5 5 5 2];

for i=1:Npp
    p_occ1 = [p_occ1 x(16+2*(N_All+1)+7*(i-1))]; %loop through Npp to get probability
end

for n=1:Npp
    k_poly1=k_poly1+((1-p_occ1(n))*npp(n)); % add to poly rate
    kp1 = [kp1 (1-p_occ1(n))*npp(n)*k_PAF*c_PA];
end
k_poly1=k_poly1*k_PAF*c_PA
kp1 = [kp1; kp1];
fil=[1 2];
figure()
sgtitle('FHOD Polymerization Rate Per Filament');
subplot(1,3,1);
bar(fil, kp1,0.5, 'stacked')
title('Single-Filament FHOD');
xlim([0.5 1.5]);
xlabel('Filaments');
ylabel('k_p_o_l_y');
ylim([0 60]);

x=importdata('FHODdouble');
%Constants
NFil=2;
k_PAF=10;
c_PA=2.5;
Npp=4;
N_All=8;
%Variables
k_poly21=0;
k_poly22=0;
p_occ1=[];
p_occ2=[];
kp1 = [];
kp2 = [];
kpboth = [];
npp=[5 5 5 2];

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
kp=[kp1 ; kp2; kpboth];
fil=[1 2 3];
subplot(1,3,2);
bar(fil, kp ,0.5, 'stacked')
title('Double-Filament FHOD');
xlabel('Filaments');
ylabel('k_{poly}');
xticklabels({1, 2, 'total'});
ylim([0 60]);

x=importdata('FHODdimer');
%Constants
NFil=2;
k_PAF=10;
c_PA=2.5;
Npp=4;
N_All=8;
%Variables
k_poly31=0;
k_poly32=0;
p_occ1=[];
p_occ2=[];
kp1 = [];
kp2 = [];
kpboth = [];
npp=[5 5 5 2];

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
fil=[1 2 3];
subplot(1,3,3);
bar(fil, kp ,0.5, 'stacked')
title('N-Dimerized FHOD');
xlabel('Filaments');
ylabel('k_{poly}');
xticklabels({1, 2, 'total'});
ylim([0 60]);

saveas(gcf, 'FHOD k-polymerization.pdf');