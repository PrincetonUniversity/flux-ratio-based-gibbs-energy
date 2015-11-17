function dgr=dgf2dgr(model,dgf)
%% Code snippet from Component Contribution in Von Bertalanffy
% by Hulda Haraldsdóttir, Elad Noor, Ronan Fleming

% Define constants
R = 8.3144621e-3; % Gas constant in kJ/(K*mol)
T = model.T; % Temperature in K
F = 96485.3365e-3; % Faraday constant in kC/mol

% Estimate standard transformed reaction Gibbs energies
St = full(model.S);
hBool = strcmp(model.metFormulas,'H');
St(hBool,:) = 0; % Set proton coefficients to 0
dgr = St' * dgf;

% Configure model.cellCompartments
model.cellCompartments = reshape(model.cellCompartments,length(model.cellCompartments),1);
if ischar(model.cellCompartments)
    model.cellCompartments = strtrim(cellstr(model.cellCompartments));
end
if isnumeric(model.cellCompartments)
    model.cellCompartments = strtrim(cellstr(num2str(model.cellCompartments)));
end

% Adjust DrGt0 for transport across membranes
metCompartmentBool = strcmp(repmat(model.metCompartments,1,length(model.cellCompartments)),repmat(model.cellCompartments',length(model.metCompartments),1));

model_nHs = zeros(size(model.mets));
for i = 1:length(model.mets)
    model_nHs(i) = numAtomsOfElementInFormula(model.metFormulas{i},'H');
end
deltaPH = full(model.S') * diag(model_nHs) * metCompartmentBool * -(R * T * log(10) * model.ph); % Adjustment due to compartmental differences in pH

model_zs = model.metCharges;
deltaCHI = full(model.S') * diag(model_zs) * metCompartmentBool * F * model.chi/1000; % Adjustment due to compartmental differences in electrical potential

dgr = dgr + deltaPH + deltaCHI;
