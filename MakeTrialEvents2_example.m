function MakeTrialEvents2_example(cellid, TE, Events_TTL, Events_TS, varargin)
%MAKETRIALEVENTS2   Synchronize trial events to recording times. 
%	MAKETRIALEVENTS2 (cellid, TEEvents_TTL, Events_TS, varargin) loads events from recording system and adjusts
%	trial event times (trial-based behavioral data) to the recorded time stamps. This
%	way the neural recordings and behavioral time stamps are in register.
%	Stimulus time TTL pulses are used for synchronization. The synchronized
%	trial events structure is saved under the name 'TrialEvents.mat'. This
%	file becomes the primary store of behavioral data for a particular
%	session; it is retrieved by LOADCB via CELLID2FNAMES. This default
%	file name is one of the preference settings of CellBase - type 
%   getpref('cellbase','session_filename');
%
%   INPUT ARGUMENTS
%   cellid     - cellid- can be any cellid of a session-MakeTrialEvents2
%                only needs to be called once per session
%   TE         - Trial Events, non-aligned
%   Events_TTL - vector with TTL values from recording system
%   Events_TS - vector with timestamps from recording system (same time
%               format as spikes)

%   OUTPUT
%   is a TrialEvents.mat in the session folder, which contains a copy of TE
%   as well as a field TrialStart with timestamps for trial start
%   from the recording system. Also saves AlignedRecEvents.mat.

%   optional pairwise arguments
%   TrialStart_TTL - double (1) - TTL pulse corresponding to trial start
%   Events_TS_conversion - double (10^6) - conversion factor to convert
%                                        time stamps from recording system to seconds.
%   
%
%   See also MAKETRIALEVENTS2_GONOGO and LOADCB.

% This default function was based on example function such as
% MAKETRIALEVENTS2_GONOGO.

% TO, CSHL, 2018

% Parse input arguments
default_args={...
    'TrialStart_TTL',      2;...
    'conversion',          10^6,...
    };
g = parse_args(default_args,varargin{:});

% Create cell aray to group TTLs and their timestamps by trial. 
inxx=find(Events_TTL==g.TrialStart_TTL);
firsttarget = inxx(1);

% Break Rec events into a cell array by trials
nTrials = 0;
nEvents = 0;
RecEvents = cell(1,1);
for x = firsttarget:length(Events_TTL)  %
    if Events_TTL(x) == g.TrialStart_TTL %trialStart
        nTrials = nTrials + 1;
        nEvents = 0;
        RecEvents{1,nTrials} = [];
        RecEvents{2,nTrials} = [];
    end
    nEvents = nEvents + 1;
    RecEvents{1,nTrials} = [RecEvents{1,nTrials} Events_TTL(x)];
    RecEvents{2,nTrials} = [RecEvents{2,nTrials} Events_TS(x)];
end

% Set TTL alignement state for first state
idx=g.TrialStart_TTL; 

% create trial segmented events
AlignedRecEvents=cell(3,length(RecEvents));

for i=1:length(RecEvents)
    AlignedRecEvents{1,i}=RecEvents{1,i};
    size(RecEvents{2,i}(AlignedRecEvents{1,i}==idx));
    AlignedRecEvents{2,i}=RecEvents{2,i}-RecEvents{2,i}(AlignedRecEvents{1,i}==idx);
    AlignedRecEvents{3,i}=RecEvents{2,i}(AlignedRecEvents{1,i}==idx);
end


% Synchronization
son = Events_TTL==idx; 

TE2 = TE;
son2 = Events_TS(son)/g.conversion;   % Trial start time recorded by the recording system (Neuralynx)
ts = TE2.TrialStartTimeStamp;   % Trial start in absolut time recorded by Bpod


% Match timestamps - in case of mismatch, try to fix
if ~ismatch(ts,son2)
    % note: obsolete due the introduction of TTL parsing
    son2 = clearttls(son2); % eliminate recorded TTL's within 0.5s from each other - broken TTL pulse
    if ~ismatch(ts,son2)
        son2 = trytomatch(ts,son2);  % try to match time series by shifting
        if ~ismatch(ts,son2)
            son2 = tryinterp(ts,son2); % interpolate missing TTL's or delete superfluous TTL's up to 10 erroneous TTl's
            if ~ismatch(ts,son2)  % TTL matching failure
                error('MakeTrialEvents:TTLmatch','Matching TTLs failed.')
            else
                warning('MakeTrialEvents:TTLmatch','Missing TTL interpolated.')
            end
        else
            warning('MakeTrialEvents:TTLmatch','Shifted TTL series.')
        end
    else
        warning('MakeTrialEvents:TTLmatch','Broken TTLs cleared.')
    end
end

% Eliminate last TTL's recorded in only one system
sto = TE2.TrialStartTimeStamp;
if length(son2) > length(ts)   % time not saved in behavior file (likely reason: autosave was used)
    son2 = son2(1:length(ts));
elseif length(son2) < length(ts)  % time not recorded on Recording (likely reason: recording stopped)
    shinx = 1:length(son2);
    ts = ts(shinx);
    sto = sto(shinx);
    TE2 = shortenTE(TE2,shinx);
    warning('MakeTrialEvents: Recording events shorter than behavioral session - will be cut.')
end

TE2.TrialStart = son2;

[subject,session] = cellid2tags(cellid);

% Save synchronized 'TrialEvents' file
if ~isempty(TE2.TrialStartTimeStamp)
    save(fullfile(getpref('cellbase','datapath'),subject,session,getpref('cellbase','session_filename')),'-struct','TE2')
else
    error('MakeTrialEvents:noOutput','Synchronization process failed.');
end

if ~isempty(AlignedRecEvents)
    save(fullfile(getpref('cellbase','datapath'),subject,session, 'AlignedRecEvents.mat'),'AlignedRecEvents')
end

% -------------------------------------------------------------------------
function I = ismatch(ts,son2)

% Check if the two time series match notwithstanding a constant drift
clen = min(length(ts),length(son2));
I = abs(max(diff(ts(1:clen)-son2(1:clen)))) < 3;  % the difference between the timestamps on 2 systems may have a constant drift, but it's derivative should still be ~0

% note: abs o max is OK, the derivative is usually a small neg. number due
% to drift of the timestamps; max o abs would require a higher tolerance
% taking the drift into account (if 2 event time stamps are far, the drift
% between them can be large)

% -------------------------------------------------------------------------
function son2 = tryinterp(ts,son2)

% Interpolate missing TTL's or delete superfluous TTL's up to 10 erroneous TTl's
for k = 1:10
    if ~ismatch(ts,son2)
        son3 = son2 - son2(1) + ts(1);
        adt = diff(ts(1:min(length(ts),length(son2)))-son2(1:min(length(ts),length(son2))));
        badinx = find(abs(adt)>0.1,1,'first') + 1;  % find problematic index
        if adt(badinx-1) < 0    % interploate
            ins = ts(badinx) - linterp([ts(badinx-1) ts(badinx+1)],[ts(badinx-1)-son3(badinx-1) ts(badinx+1)-son3(badinx)],ts(badinx));
            son2 = [son2(1:badinx-1) ins+son2(1)-ts(1) son2(badinx:end)];
        else
%             ins = son3(badinx) - linterp([son3(badinx-1) son3(badinx+1)],[son3(badinx-1)-ts(badinx-1) son3(badinx+1)-ts(badinx)],son3(badinx));
%             ts = [ts(1:badinx-1) ins ts(badinx:end)];
            son2(badinx) = [];   % delete
        end
    end
end

% -------------------------------------------------------------------------
function son2 = trytomatch(ts,son2)

% Try to match time series by shifting
len = length(son2) - 15;
minx = nan(1,len);
for k = 1:len
    minx(k) = max(diff(ts(1:15)-son2(k:k+14)));  % calculate difference in the function of shift
end
mn = min(abs(minx));
minx2 = find(abs(minx)==mn);
minx2 = minx2(1);   % find minimal difference = optimal shift
son2 = son2(minx2:min(minx2+length(ts)-1,length(son2)));

% -------------------------------------------------------------------------
function son2 = clearttls(son2)

% Eliminate recorded TTL's within 0.5s from each other
inx = [];
for k = 1:length(son2)-1
    s1 = son2(k);
    s2 = son2(k+1);
    if s2 - s1 < 0.5
        inx = [inx k+1]; %#ok<AGROW>
    end
end
son2(inx) = [];

% -------------------------------------------------------------------------
function TE2 = shortenTE(TE2,shinx)

% Eliminate behavioral trials
fnm = fieldnames(TE2);
for k = 1:length(fieldnames(TE2))
    TE2.(fnm{k}) = TE2.(fnm{k})(shinx);
end



