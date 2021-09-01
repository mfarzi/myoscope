clc; clear all; close all;
files = dir('./data/synthetic/*.mat');

for i=1:length(files)
    load(fullfile(files(i).folder, files(i).name), ...
        'scheme', 'params', 'hparams', 'sig');
    [~,modelName,~] = fileparts(files(i).name);
    model = compartment.str2model(modelName);
    model.setHyperparams(hparams);
    N = size(sig, 2);
    err = ones(N, 1);
    for n=1:N
        s = model.synthesize(params(:,n), scheme);
        err(n) = norm(sig(:,n)-s)/norm(sig(:,n));
    end
    if all(err<eps)
        fprintf('%s synth passed\n', modelName);
    else
        fprintf('%s synth failed\n', modelName);
    end
end
        