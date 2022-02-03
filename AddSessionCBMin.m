%% READ THE HELP!!!!
%cellbase is actually quite well documented.
%cellbase provides tons of function and helper functions, which can be
%easily overlooked. It is worth spending some time at least scrolling over
%the function reference and reading function headers once in a while


%% Initialize cellbase
% initcb
% deletecb
% whosecb
% choosecb

%% Get cellbase preferences
getpref('cellbase')
datapath = getpref('cellbase','datapath');

%% find all cells irrespective if they are added to cellbase
cells = findallcells();

% alternative usage e.g. for adding a new session
% cells = findallcells(path_to_session_directory);

%cf findcells()

%extract subject and session name
[subject, session] = cellid2tags(cells{1});

%% Extract TTLs from neuralynx file
% this example is specific to nlx
% events = getRawTTLs(fullfile(datapath,subject,session,'Events.nev'));

%there is no common convention for naming these files - we could think
%about one
% EventTTL=events(:,2)';
% EventTimestamps=events(:,1)';

%saving optional
% save(fullfile(datapath,subject,session,'EventTTL_task.mat'),'EventTTL');
% save(fullfile(datapath,subject,session,'EventTimestamps.mat'),'EventTimestamps');

load(fullfile(datapath,subject,session,'EventTTL_task.mat'));
load(fullfile(datapath,subject,session,'EventTimestamps.mat'));

%% Load the behavioral file
[~,BehaviorFilesFullPath] = getDir(fullfile(datapath,subject,session),'file','Dual2AFCv3');
load(BehaviorFilesFullPath{1});

TE=SessionData; 

%here it would be a good idea to add anything you want to have in
%TrialEvents struct later. For this example we need stimulus onset and
%offset times for each trial.
for tr = 1:length(TE.CorrectChoice)
TE.StimulusOnset(tr) = TE.RawEvents.Trial{tr}.States.DeliverStimulus(1);
TE.StimulusOffset(tr) = TE.RawEvents.Trial{tr}.States.DeliverStimulus(2);
end

%% Align to NLX data and save trial events file
% the alignment function *could* be fairly general I believe but right now
% there is no general function to do this. I tried to write one.

%this will save a TrialEvents.mat in session folder
% with THE ADDITIONAL FIELD ***TRIALSTART*** cointaining the trial start
% timestamps from the recording system for each trial
MakeTrialEvents2_example(cells{1}, TE, EventTTL, EventTimestamps)


%% Define EventsEpochs function
[events,epochs] = defineEventsEpochs_example; %for demo, no need to call it extra in general

%% Prealign the spikes
prealignSpikes(cells,'FUNdefineEventsEpochs',@defineEventsEpochs_example,'filetype','event','events',[],'epochs',[],'writing_behavior','overwrite','ifsave',1);


%% Create the single trial PSTH arrays
% This is something I create for every cell to save me calculation time and
% redunant calls later. But there is no general function for it. See my
% branch?
% CreateSingleTrialPSTH(cells)

%% user demo
clear all

loadcb
% CELLIDLIST

%or listtag('cells')

cellid = CELLIDLIST{1};
TE = loadcb(cellid,'TrialEvents');
SE = loadcb(cellid,'EVENTSPIKES');

%% plot/get psth functions
% viewcell2b
% ultimate_psth

%viewcell2b logic
% --> gets stimes (EVENTSPIKES)
% --> calls stimes2binraster
% --> calls binraster2psth or binraster2apsth (adaptive)
% --> plot_raster2a + plot_timecourse

%(bug: 'Partitions','all' (default) does not work because it makes too
%strong assumptions about first TE field in partition_trials.m)
%note: partition corresponds to field in TE file. has to be a double or int
%bug with legends in plot_timcourse.m
%slow (because of inefficiently plotting rasters! (with uistack))
tic
viewcell2b(cellid,'TriggerName','StimulusOnset','SortEvent','StimulusOffset','eventtype','behav','ShowEvents',{{'StimulusOffset'}},'window',[-1,1],'Partitions','#CorrectChoice','stack_events_bottom',false,'FigureNum',5);
toc

% ultimate_psth logic
% used not for plotting, but for returning psth values
% --> gets stimes (EVENTSPIKES)
% --> calls stimes2binraster
% --> calls binraster2psth or binraster2apsth or binraster2dapsth (adaptive)
% --> calls psth_stats
tic
[psth, spsth, spsth_se, tags, spt, stats] = ultimate_psth(cellid,'trial','StimulusOnset',[-1,1],'parts','#CorrectChoice','display',true);
toc
figure,plot(spsth')

% useful figure export function
%writefigs(gcf,path_to_fig) %writes to pdf and will append if pdf exists

%% what Paul and I do
% calculate single trial psth and save for different alignments
% calculate trial average psth for conditions from there
% faster but lots of customized functions

%% Add cells to Cell Base / Compute the analysis
addcell(cells); %will compute added analyses
%cf addnewcells() addnewsessions() 

%% calculate analyses
% results go in ANALYSES and TheMatrix
% suited for any analysis which assigns a number/category to each cell
% you can additionally specify input to functions, which outputs to use etc
addanalysis(@myanalysis,'property_names',{'myvalue','mypvalue'})
addanalysis(@meanrate,'property_names',{'rate'})

%some functions to deal with analyses
delanalysis(@myanalysis)

findanalysis(@myanalysis)
findprop('mypvalue')

getvalue('mypvalue')
