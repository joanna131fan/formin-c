x=importdata('outputfile');
NFil=1;
k_PAF=10;
c_PA=2.5;
k_poly=0;
p_occ=[];
Npp=4;
N_All=4;
npp=[2, 5, 5, 5];
for j=1:NFil
    for i=1:Npp
        p_occ = [p_occ x(16+2*(N_All+1)+7*(i-1))];
    end
end
for n=1:Npp
    k_poly=k_poly+((1-p_occ(n))*npp(n));
end
k_poly=k_poly*k_PAF*c_PA