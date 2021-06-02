function get_BUG_resultsmat()
% result from experiment and FEM model of CHARLIE
e_PER=[0.49, 0.55, 0.65, 0.6, 0.76, 0.8, 0.73, 0.9, 0.86];
e_PTE=[0.45, 0.48, 0.52, 0.5, 0.56, 0.57, 0.55, 0.6, 0.59];
e_Pplasma=[8.1, 17.28, 28.6, 9, 20.16, 31.35, 9.9, 21.6, 32.45];
f_PER=[6.32, 6.56, 6.68, 6.31, 6.82, 6.97, 6.9, 7.15, 7.61];
f_Rmetal=[0.38, 0.35, 0.38, 0.35, 0.34, 0.34, 0.35, 0.34, 0.33];
f_Ls=[7.47, 7.18, 7.52, 7.13, 6.98, 7.04, 7.18, 6.99, 6.87];

experiment.PER=zeros(3,1,3);
experiment.PTE=zeros(3,1,3);
experiment.Pplasma=zeros(3,1,3);
fem.PER=zeros(3,1,3);
fem.Rmetal=zeros(3,1,3);
fem.Ls=zeros(3,1,3);
i_data=0;
for i_p=1:3
    for i_Pin=1:3
        i_data=i_data+1;
        experiment.PER(i_p,1,i_Pin)=e_PER(i_data);
        experiment.PTE(i_p,1,i_Pin)=e_PTE(i_data);
        experiment.Pplasma(i_p,1,i_Pin)=e_Pplasma(i_data);
        fem.PER(i_p,1,i_Pin)=f_PER(i_data);
        fem.Rmetal(i_p,1,i_Pin)=f_Rmetal(i_data);
        fem.Ls(i_p,1,i_Pin)=f_Ls(i_data);
    end
end

fem.slitplasma=[[8.15, 0.3, 6.51];...
    [8.01, 0.28, 6.22];...
    [8.04, 0.28, 6.22]];

fem.magnetized_a=[[7.12, 0.48, 8.31];...
    [7.82, 0.48, 8.28];...
    [6.92, 0.49, 8.39]];

save('./others/for_paper210415/BUG_results.mat')
end