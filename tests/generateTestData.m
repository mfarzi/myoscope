clc; clear all;
% modelNames = {'tensor', 'ball', 'pancake', 'zeppelin', 'stick', ...
%               'cylinder', 'cylinderECS', 'cylinderGDR', 'cylinderBDA',...
%               'zeppelinBDA', 'stickBDA'};
modelNames = {'zeppelinBDA'};
scheme = myoscopeReadScheme('./data.scheme');

for i=1:length(modelNames)
    for n=1:100
        model = compartment.str2model(modelNames{i});
        model.addConstraint('s0>=0.8');
        model.addConstraint('s0<=1.2');
        model.setHyperparams(5810);
        p = model.randomInit();
        hparams = model.getHyperparams();
        s = model.synthesize(p, scheme);
        sig(:,n) = s;
        params(:,n) = p;
    end
    save(fullfile('./data/synthetic/', strcat(modelNames{i},'.mat')),...
        'params', 'hparams', 'scheme', 'sig');
    clear params sig hparams
end
        