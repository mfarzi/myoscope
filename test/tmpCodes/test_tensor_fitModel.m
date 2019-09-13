clc; clear all; close all; 

addpath('../../../camino/matlab_help_functions');
addpath('../models');
addpath('../optimizer');

%% read the scheme file and the data
scheme = camino_read_scheme('./data/scheme.scheme');
fid= fopen('./data/data.Bfloat','r','b');
data =fread(fid,'float');
fclose(fid);

%% fit a cylinder model
S0 = 1;
D = 1e-9*rand(3);

params = [S0; D(1,1); D(1,2); D(1,3); D(2,2); D(2,3); D(3, 3)];

idx = false(7,1); 
model = tensor(params);
model.fitModel(scheme, data);
%%
stickModel = ballStick([1; 0.2; 1e-9; 0.8; 1e-9; pi/2; 0]);
stickModel.fitModel(scheme, data);

%% fit a elliptical cylinder using camino
str = ['modelfit -inputfile ./data/data.Bfloat -schemefile ./data/scheme.scheme',...
           ' -model ldt -outputfile ./data/DT.Bdouble'];
system(str);
    
% rotate the scheme file around the diffusion tensor
fid= fopen('./data/DT.Bdouble','r','b');
d=fread(fid,'double');
fclose(fid);
DT = [[d(3), d(4), d(5)];...
      [d(4), d(6), d(7)];...
      [d(5), d(7), d(8)]];
camino = tensor([exp(d(2));d(3:8)]);
s_camino = stickModel.synthesize(scheme);

%%
s_model = model.synthesize(scheme);
%%
thisDirection = 6;
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
    
    yplot = data(idx);
    yplot = yplot(sortID);
    
    s = s_model(idx);
    s = s(sortID);
    
    sp = s_camino(idx);
    sp = sp(sortID);
    
    scatter(xplot.^2, yplot, ones(length(yplot), 1)*100, 'filled', colorCode{thisCounter});
    %scatter(xplot.^2, s, ones(length(yplot), 1)*100,'*', 'k');
    semilogy(xplot.^2, s, colorCode{thisCounter}, 'linewidth', 2);
    semilogy(xplot.^2, sp, '--k', 'linewidth', 2);
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