function m = latMetrics(latR, latT)
% -------------------------------------------------------------------------
% Calculate metrics to quantify the differences between a reference (true)
% and a reconstructed local activation time (LAT) map.
% -------------------------------------------------------------------------
%
% Inputs: latR: reconstructed LAT map (NxL)
%         latT: true (reference) LAT map (NxL)
%         N is the number of nodes; L is the max number of LATs per node
%         Unused entries must be NaNs
%
% Output: m: struct containing the fields:
%            TP:        number of true positives
%            FN:        number of false negatives
%            FP:        number of false positives
%            FNR:       false negative rate in %
%            FPR:       false positive rate in %
%            DSC:       Dice similarity coefficient
%            AEmean:    mean of absolute errors
%            AEstd:     standard deviation of absolute errors
%            truePosT:  indices of true positives in latT (1xN cell)
%            truePosR:  indices of true positives in latR (1xN cell)
%            falseNegT: indices of false negatives in latT (1xN cell)
%            falsePosR: indices of false positives in latR (1xN cell)
%            absErr:    individual absolute errors (1xN cell)
%
% - True positives are defined as the closest correspondences between 
%   rows in latR and rows in latT.
% - False negatives are the elements in latT after removing true positives.
% - False positives are the elements in latR after removing true positives.
% - Absolute errors are calculated between only true positives.
%
% -------------------------------------------------------------------------
% Steffen Schuler, August 2017
% Institute of Biomedical Engineering
% Karlsruhe Institute of Technology
% www.ibt.kit.edu
% -------------------------------------------------------------------------

m.TP = 0;
m.FN = 0;
m.FP = 0;

for i = 1:size(latR,1)

    T = latT(i,:)'; T(isnan(T)) = [];
    R = latR(i,:)'; R(isnan(R)) = [];

    truePosT = []; truePosR = [];
    falseNegT = []; falsePosR = [];
    
    if(numel(T) && ~numel(R))
        falseNegT = T;
    
    elseif(numel(R) && ~numel(T))
        falsePosR = R;
        
    else
        % determine true positives
        % by finding closest correspondences between T and R
        truePosT = 1:numel(T);
        truePosR = 1:numel(R);
        kT = 0; kR = 1;
        while(~isequal(kT,kR))
            kT = dsearchn(T(truePosT), R(truePosR));
            kR = dsearchn(R(truePosR), T(truePosT));
            truePosT = truePosT(unique(kT));
            truePosR = truePosR(unique(kR));
        end

        % determine false negatives and false positives
        % by removing true positives from T and R, respectively
        falseNegT = 1:numel(T); falseNegT(truePosT) = [];
        falsePosR = 1:numel(R); falsePosR(truePosR) = [];
    end
    
    m.TP = m.TP + numel(truePosT);
    m.FN = m.FN + numel(falseNegT);
    m.FP = m.FP + numel(falsePosR);
    m.truePosT{i} = truePosT;
    m.truePosR{i} = truePosR;
    m.falseNegT{i} = falseNegT;
    m.falsePosR{i} = falsePosR;
    m.absErr{i} = abs(R(truePosR)-T(truePosT));

    % % For debugging
    % hold off
    % plot(T(truePosT), ones(numel(truePosT),1), 'go')
    % hold on
    % plot(R(truePosR), -ones(numel(truePosR),1), 'gx')
    % plot(T(falseNegT), ones(numel(falseNegT),1), 'o')
    % plot(R(falsePosR), -ones(numel(falsePosR),1), 'x')
    % xlim([min(latT(:)) max(latT(:))])
    % ylim([-20 20])
    % title(i)
    % waitforbuttonpress
    
end

m.FNR = 100 * m.FN / numel(find(~isnan(latT)));
m.FPR = 100 * m.FP / numel(find(~isnan(latR)));
m.DSC = 2*m.TP / (2*m.TP + m.FN + m.FP);
AE = cell2mat(m.absErr(:));
m.AEmean = mean(AE(:));
m.AEstd = std(AE(:));

end