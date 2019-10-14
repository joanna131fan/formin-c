
m=importdata('outputfile');
m(1);
NFil = m(2);
isitetot = 9;
N_All = 9; %N_AlliSites
p_occlude=[];
iSite = [];
ree = [];
ree2 = [];
for j=1:NFil
    for i=1:isitetot
        iSite = [iSite m(15+2*(N_All+1)+7*(i-1))];
        p_occlude = [p_occlude m(16+2*(N_All+1)+7*(i-1))];
        ree = [ree m(26 + 2*(N_All +1) + 7*(isitetot-1) + 2*(NFil-1) + 2*(i-1))];
        ree2 = [ree2 m(27 + 2*(N_All +1) + 7*(isitetot-1) + 2*(NFil-1) + 2*(i-1))];

    end
end
figure()
plot(iSite, p_occlude)
figure()
plot(iSite, ree)
figure()
plot(iSite, ree2)
