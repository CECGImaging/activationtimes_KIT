function lat = latDetection(sig, Fs, varargin)
% -------------------------------------------------------------------------
% Detect local activation times (LATs) by finding the peaks in the
% (positive or negative) temporal derivative.
% -------------------------------------------------------------------------
%
% Inputs: sig: input signals (NxL)
%              N is the number of nodes; L is the signal length
%         Fs:  sampling frequency in Hz
%         name-value pairs:
%         'vertices': list of vertices coords (Nx3)
%         'faces':    list of triangle vertices (Mx3)
%         'Edge':     'rising' (default) or 'falling' (edge to be detected)
%         'Mode':     'temporal' or 'spatiotemporal' (default)
%                     spatiotemporal mode: the derivative is additionally
%                     weighted with the norm of the spatial gradient
%         'MinPeakHeight':     min peak height in units of std(sig(i,:))
%                              Default: 0.5
%         'MinPeakDistance':   min distance between peaks in ms
%                              Default: 100ms
%         'MinPeakProminence': see documentation of findpeaks()
%                              Default: 0.3
%
% Output: lat: local activation times (NxL)
%              L is the max number of LATs per node
%              Unused entries are filled with NaNs
%
% -------------------------------------------------------------------------
% Steffen Schuler, August 2017
% Institute of Biomedical Engineering
% Karlsruhe Institute of Technology
% www.ibt.kit.edu
% -------------------------------------------------------------------------

p = inputParser;
addOptional(p,'vertices',[]);
addOptional(p,'faces',[]);
addParameter(p,'Edge','rising');
addParameter(p,'Mode','spatiotemporal');
addParameter(p,'MinPeakHeight',0.5);
addParameter(p,'MinPeakDistance',100);
addParameter(p,'MinPeakProminence',0.3);
parse(p, varargin{:});

switch p.Results.Edge
    case 'rising'
        sign = 1;
    case 'falling'
        sign = -1;
    otherwise
        error('Parameter ''Edge'' must be one of {''rising'',''falling''}.');
end

switch p.Results.Mode
    case 'temporal'
        peakSig = sign .* diff(sig,1,2);
    case 'spatiotemporal'
        if(isempty(p.Results.vertices) || isempty(p.Results.faces))
            error('Vertices and faces must be given for ''spatiotemporal'' mode.');
        end
        mesh = PrepareTriangleMesh(p.Results.vertices, p.Results.faces);
        [Gx,Gy,Gz] = Gradient(mesh);
        gradSig = sqrt((Gx*sig).^2+(Gy*sig).^2+(Gz*sig).^2);
        peakSig = sign * gradSig(:,1:end-1) .* diff(sig,1,2);
    otherwise
        error('Parameter ''Mode'' must be one of {''temporal'',''spatiotemporal''}.');
end

Ts = 1e3 / Fs; % sampling period in ms
minPeakDist = round(p.Results.MinPeakDistance / Ts);

lat = NaN(size(sig));
for i = 1:size(sig,1)

    s = peakSig(i,:);
    s = max(p.Results.MinPeakHeight, s/std(s));
    
    [pks,locs] = findpeaks(s, 'MinPeakDistance',minPeakDist, 'MinPeakProminence',p.Results.MinPeakProminence);
    lat(i,1:numel(locs)) = Ts * (locs-1);

    % % For debugging
    % hold off
    % findpeaks(s, 'MinPeakDistance',minPeakDist, 'MinPeakProminence',p.Results.MinPeakProminence)
    % hold on
    % plot((sig(i,:)-mean(sig(i,:))+2*std(sig(i,:)))./std(sig(i,:)))
    % text(locs+.02,pks,num2str((1:numel(pks))'))
    % title(i)
    % waitforbuttonpress

end
firstNanCol = find(isnan(min(lat,[],1)), 1, 'first');
lat(:,firstNanCol:end) = [];

end