function get_CHARLIE_resultsmat()
% result from experiment and FEM model of CHARLIE
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

fem.nonuniform_in4_PER=[0.8, 0.4, 0.31, 0.27, 0.26, 0.25];
fem.nonuniform_in4_Rmetal=[0.047, 0.047, 0.045, 0.044, 0.044, 0.044];
fem.nonuniform_in4_Ls=[2.62, 2.64, 2.59, 2.56, 2.55, 2.56];
fem.nonuniform_in2_PER=[0.73, 0.36, 0.28, 0.26, 0.25, 0.25];
fem.nonuniform_in2_Rmetal=[0.042, 0.042, 0.039, 0.041, 0.041, 0.042];
fem.nonuniform_in2_Ls=[2.5, 2.52, 2.43, 2.48, 2.49, 2.51];
fem.nonuniform_in4_ra_np5=[0.38, 0.044, 2.55];

fem.sweep_ne_PER=[0.04, 0.05, 0.09, 0.14, 0.21, 0.33, 0.38, 0.51, 0.77, 0.79, 1.2, 1.78, 2.46];
fem.sweep_ne_Rmetal=[0.046962502, 0.046972295, 0.046965189, 0.046956217, 0.046968535, 0.046959449, 0.046961613, 0.046926009, 0.046930069, 0.046941435, 0.046982031, 0.047337354, 0.047399779];
fem.sweep_ne_Ls=[2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.64, 2.64, 2.63, 2.63, 2.61, 2.58, 2.5];

save('./others/for_paper210415/CHARLIE_results.mat')
end