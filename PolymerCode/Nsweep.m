m=dlmread('output1_25_300double_b.txt');

NFil=2;
%POcc1 = (-1)*ones(size(m,1) , 276);
p_occ_mat1=[];
p_occ_mat2=[];
p_occ_1=[];
p_occ_2=[];
N_mat =[];
N_array=1:25:101;
isite_mat=[];%shorthand?
N=1;
cnt=1;
%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
while(N<=126)
    N_All=NFil*N;
    isitetot=N;
    line=m(cnt, :);
    for isite=1:N
        temp1=line(16+2*(N_All+1)+7*(isite-1));
        if temp1>1
            temp1=1;
        end
        p_occ_1=[p_occ_1 temp1];
        temp2=line(16+2*(N_All+1)+7*(isite-1)+ (6 + 9*N + 2 + NFil + NFil));
        if temp2>1
            temp2=1;
        end
        p_occ_2=[p_occ_2 temp2];
        isite_mat=[isite_mat isite/N];
        N_mat=[N_mat N];
        %POcc1(iN, isite) =  m(16+2*(N_All+1)+7*(isite-1));
        %p_occludea2 = [p_occludea2 m(16+2*(N_All+1)+7*(isite-1)+ (6 + 7*isitetot + 2 + NFil + NFil + 2*isitetot))];
    end
    %POcc1(iN,isite) = [POcc1; p_occludea1];
    N=N+25;
    cnt=cnt+1;
end

figure()

 a = size(transpose(N_mat));   
 a = a(1);
 A = zeros(a,3);
 A(:,1) = transpose(N_mat);
 A(:,2) = transpose(isite_mat);
 A(:,3) = transpose(p_occ_1);
 

 disp(A);
 
 N_coarse = A(:,1);
 isite_coarse = A(:,2);
 p_occ_coarse = A(:,3);
 
 N_fine=linspace(min(N_coarse),max(N_coarse),130);
 isite_fine=linspace(min(isite_coarse),max(isite_coarse),130);
 [N_finemesh,isite_finemesh]=meshgrid(N_fine,isite_fine);
 F=scatteredInterpolant(N_coarse,isite_coarse,p_occ_coarse, 'linear', 'none');
 contourf(N_finemesh,isite_finemesh,F(N_finemesh,isite_finemesh), 100,'LineColor','none') 
 colorbar;
 colormap;
 xlabel('Length of FH1 domain')
 ylabel('Binding Site Location')
 caxis([0 1])
 

figure()

 a = size(transpose(N_mat));   
 a = a(1);
 A = zeros(a,3);
 A(:,1) = transpose(N_mat);
 A(:,2) = transpose(isite_mat);
 A(:,3) = transpose(p_occ_2);
 
 disp(A)
 
 N_coarse = A(:,1);
 isite_coarse = A(:,2);
 p_occ_coarse = A(:,3);

 N_fine=linspace(min(N_coarse),max(N_coarse),130);
 isite_fine=linspace(min(isite_coarse),max(isite_coarse),130);
 [N_finemesh,isite_finemesh]=meshgrid(N_fine,isite_fine);
 F=scatteredInterpolant(N_coarse,isite_coarse,p_occ_coarse, 'linear','none');
 contourf(N_finemesh,isite_finemesh,F(N_finemesh,isite_finemesh),100, 'LineColor','none') 
 colorbar;
 colormap;
 xlabel('Length of FH1 domain')
 ylabel('Binding Site Location')
 caxis([0 1])