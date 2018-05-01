function [out1, out2] = myanalysis(cellid,varargin)
%
%  MYANALYSIS     Calculates examaple analysis for a cell
%
%
%   x = meanrate(cellid,{isi_threshold})
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Analysis functions should always take the {cellid} as
%  the first argument and can possibly take a variable argument
%  list.
%  The output can be a *single* scalar or a vector of numerical values.
%  (i.e. don't give multiple outputs, put it into a vector)
%
%  Calling stimes = loadcell(cellid,'caller'); initially will load
%  the necessary data.
%
%  To support self-describing analysis functions return a cell array
%  of property descriptor strings back if cellid is 'default'.
%

TE = loadcb(cellid,'TrialEvents');
SE = loadcb(cellid,'EVENTSPIKES');


stim_rates = SE.epoch_rates{2};

outcome_index = (nanmean(stim_rates(TE.CorrectChoice==1)) - nanmean(stim_rates(TE.CorrectChoice==0))) ./ nanmean(stim_rates);
p = ranksum(stim_rates(TE.CorrectChoice==1),stim_rates(TE.CorrectChoice==0));

out1 = outcome_index;
out2 = p;
