clc; clear all; close all;
files = dir('./data/synthetic/*.mat');

for i=5%:length(files)
    load(fullfile(files(i).folder, files(i).name), ...
        'scheme', 'params', 'hparams', 'sig');
    [~,modelName,~] = fileparts(files(i).name);
    model = compartment.str2model(modelName);
    model.addConstraint('s0>=0.8');
    model.addConstraint('s0<=1.2');
%     model.setHyperparams(hparams');
    N = size(sig, 2);
    P = ones(model.nParams, N);
    for n=1:N
        P(:,n) = model.fit(sig(:,n), scheme);
    end
%     if all(err<eps)
%         fprintf('%s synth passed\n', modelName);
%     else
%         fprintf('%s synth failed\n', modelName);
%     end
end
        