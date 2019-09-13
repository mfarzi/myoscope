clc; clear all; close all; 

addpath('../../../../camino/matlab_help_functions');
addpath('../../models');
addpath('../../');

%% read the scheme file and the data
load('../testData/experimental/scheme.mat');
dataID = 2;
load(sprintf('../testData/experimental/data_%02d.mat', dataID));
data = double(data);
% normalise for b0 signal
i = 0;
b0Signal = zeros(5, 1);
for thisDELTA = 0.01:0.01:0.05
    i = i + 1;
    idx = scheme.DELTA==thisDELTA;
    ref = scheme.DELTA==thisDELTA & scheme.G_dir==0;
    
    b0Signal(i) = mean(data(ref),2);
    data(idx) = data(idx)/b0Signal(i);
end
data = double(data);

%% fit models
model1 = tensor();
p1 = model1.fitMultiRun(scheme, data, 5);
S1 = model1.synthesize(scheme);

model2 = twoCompartmentModel(cylinder(), zeppelin());
p2 = model2.fitMultiRun(scheme, data, 10);
S2 = model2.synthesize(scheme);

model3 = threeCompartmentModel(ellipticalCylinder(), ball(), stick());
p3 = model3.fitMultiRun(scheme, data, 20);
S3 = model3.synthesize(scheme);

% model4 = threeCompartmentModel(...
%          ellipticalCylinder('alpha', model1.alpha, 'theta', model1.theta, 'phi', model1.phi), ...
%          ball('diff', 2.3e-9), ...
%          stick('diff', 2.3e-9, 'theta', model1.theta, 'phi', model1.phi));
% model4.comp1.fixParams('alpha', true, 'theta', true, 'phi', true);
% model4.comp2.fixParams('diff', true);
% model4.comp3.fixParams('diff', true, 'theta', true, 'phi', true);

model4 = threeCompartmentModel(cylinder(), ball(), stick());

p4 = model4.fitMultiRun(scheme, data, 20);
S4 = model4.synthesize(scheme);

%%
thisDirection = 3;
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
    
    s1 = S1(idx);
    s1 = s1(sortID);
    
    s2 = S2(idx);
    s2 = s2(sortID);
    
    s3 = S3(idx);
    s3 = s3(sortID);
    
    s4 = S4(idx);
    s4 = s4(sortID);
    
%     
    scatter(xplot.^2, yplot, ones(length(yplot), 1)*100, 'filled', colorCode{thisCounter});
    
    h1 = semilogy(xplot.^2, s1, 'k', 'linewidth', 2);
    h2 = semilogy(xplot.^2, s2, 'b', 'linewidth', 2);
    h3 = semilogy(xplot.^2, s3, 'r', 'linewidth', 2);
    h4 = semilogy(xplot.^2, s4, 'g', 'linewidth', 2);

end
ylim([0.05 1])
xlim([0 0.8])
set(gca, 'yscale', 'log')
xlabel('Squared diffusion gradient strength, (T/m)^2', 'FontSize', 20)
ylabel('Signal Attenuation', 'FontSize', 20)
legend([h1, h2, h3, h4], {model1.name, model2.name, model3.name, model4.name}, 'interpreter', 'none');%{'Delta=10ms', 'Delta=20ms', 'Delta=30ms', 'Delta=40ms', 'Delta=50ms'});
set(gca, 'LineWidth', 2)
set(gca, 'FontSize', 20)
axis square
box on 