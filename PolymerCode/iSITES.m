m=importdata('outputfileFHOD');
NFil = m(2);
isitetot = 4;
N_All = 8; %N_AlliSites
figure();
    p_occludea1=[];
    p_occludea2 =[];
    iSitea1 = [];
    iSitea2 = [];
    reea1 = [];
    reea2 = [];
    ree2a1 = [];
    ree2a2 =[];
    for i=1:isitetot
        iSitea1 = [iSitea1 m(15+2*(N_All+1)+7*(i-1))];
        p_occludea1 = [p_occludea1 m(16+2*(N_All+1)+7*(i-1))];
        reea1 = [reea1 m(26 + 2*(N_All +1) + 7*(isitetot-1) + 2*(NFil-1) + 2*(i-1))];
        ree2a1 = [ree2a1 m(27 + 2*(N_All +1) + 7*(isitetot-1) + 2*(NFil-1) + 2*(i-1))];
        iSitea2 = [iSitea2 m(15+2*(N_All+1)+7*(i-1)+ (6 + 7*isitetot + 2 + NFil + NFil + 2*isitetot))];
        p_occludea2 = [p_occludea2 m(16+2*(N_All+1)+7*(i-1)+ (6 + 7*isitetot + 2 + NFil + NFil + 2*isitetot))];
        reea2 = [reea2 m(26 + 2*(N_All +1) + 7*(isitetot-1) + 2*(NFil-1) + 2*(i-1)+ (6 + 7*isitetot + 2 + NFil + NFil + 2*isitetot))];
        ree2a2 = [ree2a2 m(27 + 2*(N_All +1) + 7*(isitetot-1) + 2*(NFil-1) + 2*(i-1)+ (6 + 7*isitetot + 2 + NFil + NFil + 2*isitetot))];
    end
    subplot(2, 3, 1);
    hold on;
    plot(iSitea1, p_occludea1, '-r')
    plot(iSitea2, p_occludea2, '-c')
    xlabel("Polyproline Binding Sites");
    ylabel("Probability of Occlusion");
    hold off;
    subplot(2, 3, 2);
    hold on;
    plot(iSitea1, reea1, '-r')
    plot(iSitea2, reea2, '-c')
    xlabel("Polyproline Binding Sites");
    ylabel("<r_e_e");
    hold off;
    
    subplot(2, 3, 3);
    hold on;
    plot(iSitea1, ree2a1, '-r')
    plot(iSitea2, ree2a2, '-c')
    xlabel("Polyproline Binding Sites");
    ylabel("<r_e_e^2>");
    hold off;

a=importdata('outputfileCAPU');
NFil = a(2);
isitetot = 5;
N_All = 10; %N_AlliSites
    p_occludeb1=[];
    p_occludeb2=[];
    iSiteb1 = [];
    iSiteb2=[];
    reeb1 = [];
    reeb2 = [];
    ree2b1 = [];
    ree2b2 =[];
    for i=1:isitetot
        iSiteb1 = [iSiteb1 a(15+2*(N_All+1)+7*(i-1))];
        p_occludeb1 = [p_occludeb1 a(16+2*(N_All+1)+7*(i-1))];
        reeb1 = [reeb1 a(26 + 2*(N_All +1) + 7*(isitetot-1) + 2*(NFil-1) + 2*(i-1))];
        ree2b1 = [ree2b1 a(27 + 2*(N_All +1) + 7*(isitetot-1) + 2*(NFil-1) + 2*(i-1))];
        iSiteb2 = [iSiteb2 a(15+2*(N_All+1)+7*(i-1)+ (6 + 7*isitetot + 2 + NFil + NFil + 2*isitetot))];
        p_occludeb2 = [p_occludeb2 a(16+2*(N_All+1)+7*(i-1)+ (6 + 7*isitetot + 2 + NFil + NFil + 2*isitetot))];
        reeb2 = [reeb2 a(26 + 2*(N_All +1) + 7*(isitetot-1) + 2*(NFil-1) + 2*(i-1)+ (6 + 7*isitetot + 2 + NFil + NFil + 2*isitetot))];
        ree2b2 = [ree2b2 a(27 + 2*(N_All +1) + 7*(isitetot-1) + 2*(NFil-1) + 2*(i-1)+ (6 + 7*isitetot + 2 + NFil + NFil + 2*isitetot))];
    end
    
    subplot(2, 3, 4);
    hold on;
    plot(iSiteb1, p_occludeb1, '-b');
    plot(iSiteb2, p_occludeb2, '-m');
    xlabel("Polyproline Binding Sites");
    ylabel("Probability of Occlusion");
    hold off;
    
    subplot(2, 3, 5);
    hold on;
    plot(iSiteb1, reeb1, '-b');
    plot(iSiteb2, reeb2, '-m');
    xlabel("Polyproline Binding Sites");
    ylabel("<r_e_e>");
    hold off;
    
    subplot(2, 3, 6);
    hold on;
    plot(iSiteb1, ree2b1, '-b');
    plot(iSiteb2, ree2b2, '-m');
    xlabel("Polyproline Binding Sites");
    ylabel("<r_e_e^2>");
    hold off;
saveas(gcf, 'doubleFHOD&CAPU.pdf');