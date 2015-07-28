function [ rxn_index ] = getRxnIndex(model, rxn)
%getRxnIndex Finds the index of the reaction specified as an argument in the
%metabolic model.
%   The index is the nth 1 of the 7th column of model.fluxes for the rxn
%   specified as an argument. The row of this reaction is its corresponding
%   column in model.rxns. (The 7th column of model.fluxes corresonds to the
%   free net fluxes, so this is in effect finding the nth free net flux.)

% Number of reactions in model
numRxns = length(model.rxns);

%Initialize index into the specified rxn name
rxn_row = 0;

% Determine rxn row for rxn in model.fluxes
for i=1:numRxns
if strcmp(model.rxns{i}, rxn) == 1
    rxn_row = i;
    break
end
end

% Return errors if rxn is not found or it is not a free flux
if rxn_row == 0
    error('Reaction name not found. Please check syntax and try again.')
end
if model.fluxes(rxn_row, 7) ~= 1
    error('Specified Reaction is not a free flux.')
end

% Determine the net free flux index of the reaction of interest
rxn_index = 0;
for j =1:rxn_row
    if model.fluxes(j, 7) == 1
        rxn_index = rxn_index + 1;
    end
end
end

