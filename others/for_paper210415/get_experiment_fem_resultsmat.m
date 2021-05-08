% result mat
clear
experiment.PER=[0.1, 0.35, 0.6, 0.51, 0.37, 0.17, 1.43, 2.38, 3.24, 5.35, 5.91, 5.23];
experiment.PTE=[0.51, 0.79, 0.87, 0.85, 0.81, 0.66, 0.85, 0.92, 0.94, 0.97, 0.97, 0.96];
experiment.PER=[experiment.PER(1:6)',experiment.PER(7:12)'];
experiment.PTE=[experiment.PTE(1:6)',experiment.PTE(7:12)'];

fem.dielectric_PER=[0.9, 1.29, 1.33, 1.01, 0.71, 0.38, 4.23, 5.32, 6.83, 8.33, 8.59, 7.5];
fem.dielectric_Rmetal=[0.047, 0.047, 0.047, 0.047, 0.047, 0.047, 0.095, 0.095, 0.095, 0.094, 0.092, 0.092];
fem.dielectric_Ls=[2.58, 2.55, 2.57, 2.61, 2.63, 2.64, 2.34, 2.28, 2.24, 2.33, 2.36, 2.47];
fem.conductor_PER=[1.24, 1.66, 1.56, 1.06, 0.73, 0.39, 13.88, 14.14, 14.15, 12.31, 11.22, 8.49];
fem.conductor_Rmetal=[0.047, 0.047, 0.047, 0.047, 0.047, 0.047, 0.095, 0.095, 0.095, 0.092, 0.092, 0.092];
fem.conductor_Ls=[2.61, 2.6, 2.61, 2.63, 2.64, 2.65, 2.24, 2.18, 2.18, 2.33, 2.39, 2.5];

field = fieldnames(fem); % cell
for i = 1:length(field)
    name_i = field{i};
    value_i = fem.(name_i); % 方式二
    fem.(name_i)=[value_i(1:6)',value_i(7:12)'];
end

fem.nonuniform_PER=[0.8, 0.4, 0.31, 0.27, 0.26, 0.25, 0.24];
fem.nonuniform_Rmetal=[0.047, 0.047, 0.045, 0.044, 0.044, 0.044, 0.042];
fem.nonuniform_Ls=[2.62, 2.64, 2.59, 2.56, 2.55, 2.56, 2.5];
fem.nonuniform_mesh_plasma=[2.62, 2.64, 2.59, 2.56, 2.55, 2.56, 2.5];

save('./others/for_paper210415/experiment_fem_results.mat')