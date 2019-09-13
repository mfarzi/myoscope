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
D = 1e-9;
theta = pi/2;
phi = 0; 
alpha = 0;
R1 = 10e-6;
R2 = 10e-6;
params = [S0; D; R1; R2; theta; phi; alpha];

idx = false(5,1); 
model = ellipticalCylinderGPD(params);

model.fitModel(scheme, data);

%% fit a elliptical cylinder using camino
nReps = 10;
str = ['modelfit -inputfile ./data/data.Bfloat -inputdatatype float ',...
       '-fitalgorithm MULTIRUNLM -samples ', num2str(nReps), ' ', ... %'-fitalgorithm LM ', ... 
       '-fitmodel ELLIPTICALCYLINDERONLY ', ...
       '-schemefile ./data/scheme.scheme ',...
       '-noisemodel Gaussian ', ...
       '-voxels 1 ',...
       '-outputfile ./data/params.Bdouble'];
system(str);

fid= fopen('./data/params.Bdouble','r','b');
d=fread(fid,'double');
fclose(fid);
D = reshape(d, [length(d)/nReps, nReps]);
if sum(D(1,:)==0) == 0
    [~,iMIN] = min(D(end,:));
else
    D(:,D(1,:)~=0)=[];
    [~,iMIN] = min(D(end,:));
end
params_camino = D(:,iMIN);

str = sprintf(['datasynth -synthmodel compartment 1 ' ...
               'ELLIPTICALCYLINDER %d %d %d %d %d %d ' ...
               '-schemefile ./data/scheme.scheme ',...
               '-voxels 1 -outputfile ./data/data_synth_camino.Bfloat'], ...
               params_camino(4), params_camino(5), params_camino(6), params_camino(7), params_camino(8), params_camino(9)); % 
           
system(str);     

fid= fopen('./data/data_synth_camino.Bfloat','r','b');
s_camino=fread(fid,'float');
fclose(fid);
    

%%
%cylinder = cylinderGPD(p);
s_model = model.synthesize(scheme);
%%
thisDirection = 10;
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