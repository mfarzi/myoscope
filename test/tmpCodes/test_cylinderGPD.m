% this script validate the a cylinder model developed in Matlab agains the
% model in camino
clc; clear all; close all; 

addpath('../../../camino/matlab_help_functions');
addpath('../models');
scheme = camino_read_scheme('myscheme.scheme');

S0 = 1;
D = 1.2e-9;
theta = pi/2;
phi = 0; 
R = 10e-6;
params = [S0; D; R; theta; phi];

%% camino
str = sprintf(['datasynth -synthmodel compartment 1 ' ...
               'CYLINDERGPD %d %d %d %d ' ...
               '-schemefile myscheme.scheme ',...
               '-voxels 1 -outputfile data_synth.Bfloat'], ...
               D, theta, phi, R); % 
           
system(str);     

fid= fopen('data_synth.Bfloat','r','b');
s_camino=fread(fid,'float');
fclose(fid);

%%
cylinder = cylinderGPD(params);
s_model = cylinder.synthesize(scheme);

%% plotting
thisDirection = 2;
colorCode = {'c', 'b', 'r', 'k', 'g'}; 
plot_handles = [];
thisCounter = 0;
figure('position', [100 100 700 650]);
hold on;
for thisDELTA = 0.01:0.01:0.05
    thisCounter = thisCounter +1;
     
    idx = scheme.DELTA==thisDELTA & scheme.G_dir==thisDirection;
    ref = scheme.DELTA==thisDELTA & scheme.G_dir==0;

    
    xplot = scheme.G_mag(idx);
    [~,sortID] = sort(xplot);
    xplot = xplot(sortID);
    
    yplot = s_camino(idx);
    yplot = yplot(sortID);
    
    s = s_model(idx);
    s = s(sortID);
    
    scatter(xplot.^2, yplot, ones(length(yplot), 1)*100, 'filled', colorCode{thisCounter});
    scatter(xplot.^2, s, ones(length(yplot), 1)*100,'*', 'k');
%     semilogy(xplot.^2, s, colorCode{thisCounter}, 'linewidth', 2);
end
ylim([0.05 1])
xlim([0 0.8])
set(gca, 'yscale', 'log')
xlabel('Squared diffusion gradient strength, (T/m)^2', 'FontSize', 20)
ylabel('Signal Attenuation', 'FontSize', 20)
legend(plot_handles, {'Delta=10ms', 'Delta=20ms', 'Delta=30ms', 'Delta=40ms', 'Delta=50ms'});
set(gca, 'LineWidth', 2)
set(gca, 'FontSize', 20)
axis square
box on 

%%
figure; scatter((s_model+s_camino)/2, s_model-s_camino)