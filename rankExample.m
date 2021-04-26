% please start CellNetAnalyzer before running this script
if ~exist('cnan','var')
    error('Please start CellNetAnalyzer first (''startcna(1)'')');
end

%% loading variables, defining path (also for saving returned values)
if ~exist('filename','var')
    filename = 'iJO1366-1-4-Butanediol';
    %filename = 'iJO1366-L-Methionine';
end
if ~exist('ext','var')
    ext = '.mat';
end
[path,filename,ext] = fileparts(which([filename ext]));
path = [path '/'];
load([filename ext],'mcs','cnap','T','D','t','d','coreReacs');
if ~exist('solver','var')
    solver = 0;
end

%% defining weighting coefficients
weightCoeff = [     1.5 , ...   #intv
                    1, ...      max growth
                    1, ...      min Y @ max growth
                    1, ...      min Y
                    0.5, ...    O2
                    0.5, ...    alt prods
                    0, ...      access met
                    0.5, ...    optMDF
                    0.5, ...    overlap
                    0.25]; %    core model

%% preparing structure that contains the the most important reaction indices
idx.subs   = find(strcmp(cellstr(cnap.reacID), 'EX_glc__D_e'));
idx.prod   = find(xor(T(1,:)',full(sparse(idx.subs,1,1,cnap.numr,1))));
idx.atpm   = find(strcmp(cellstr(cnap.reacID), 'ATPM'));
idx.o2     = find(strcmp(cellstr(cnap.reacID), 'EX_o2_e'));
idx.bm     = find(D(1,:));

%% translate cut sets into lower and upper bounds
mcs = mcs(1:4,:); %% For testing only 4 MCS! If you want to rank all MCS then outcomment this line!
mcs = logical([full(sparse(1,idx.o2,1,1,cnap.numr)); mcs]);
IS_to_rank = 2:size(mcs,1); % 1st is anaerobic reference
% upper and lower bounds for each model
IS_lb = repmat({cnap.reacMin},1,size(mcs,1));
IS_ub = repmat({cnap.reacMax},1,size(mcs,1));
for j = 1:size(mcs,1)
    IS_lb{j}(mcs(j,:)) = 0;
    IS_ub{j}(mcs(j,:)) = 0;
end
%% preparing parameters for: P06 accessible Metabolites
cytMet = regexp(cellstr(cnap.specID), '.*_c$', 'match');
cytMet = find(~cellfun(@isempty,cytMet));
idx.atp   = find(strcmp(cellstr(cnap.specID), 'atp_c'));
idx.adp   = find(strcmp(cellstr(cnap.specID), 'adp_c'));
idx.amp   = find(strcmp(cellstr(cnap.specID), 'amp_c'));
idx.pi    = find(strcmp(cellstr(cnap.specID), 'pi_c'));
idx.h2o   = find(strcmp(cellstr(cnap.specID), 'h2o_c'));
idx.h     = find(strcmp(cellstr(cnap.specID), 'h_c'));
idx.nad   = find(strcmp(cellstr(cnap.specID), 'nad_c'));
idx.nadh  = find(strcmp(cellstr(cnap.specID), 'nadh_c'));
idx.nadp  = find(strcmp(cellstr(cnap.specID), 'nadp_c'));
idx.nadph = find(strcmp(cellstr(cnap.specID), 'nadph_c'));

%% preparing parameters for: P08 MDF
idx.glc_e = find(strcmp(cellstr(cnap.specID), 'glc__D_e'));
idx.co2_e = find(strcmp(cellstr(cnap.specID), 'co2_c'));
mdfParam.G0      = cell2mat(CNAgetGenericReactionData_as_array(cnap,'deltaGR_0'));
mdfParam.uncert  = cell2mat(CNAgetGenericReactionData_as_array(cnap,'uncertGR_0'));
if size(mdfParam.G0,1) < cnap.numr % delta G is NaN for reactions that were added without thermd. data
    mdfParam.G0(end:cnap.numr) = NaN;
end
if size(mdfParam.uncert,1) < cnap.numr
    mdfParam.uncert(end:cnap.numr) = NaN;
end
mdfParam.Cmin    = 1e-6*ones(cnap.nums,1);
mdfParam.Cmin(idx.glc_e) = 1e-6;
mdfParam.Cmax = 0.02*ones(cnap.nums,1);
mdfParam.Cmax(idx.co2_e) = 1e-4;
mdfParam.Cmax(idx.glc_e) = 0.055557;
mdfParam.fixed_ratios(1,1:3) = [idx.atp   idx.adp   10];
mdfParam.fixed_ratios(2,1:3) = [idx.adp   idx.amp    1];
mdfParam.fixed_ratios(3,1:3) = [idx.nad   idx.nadh  10];
mdfParam.fixed_ratios(4,1:3) = [idx.nadph idx.nadp  10];
mdfParam.RT = 8.31446*300/1000; % Computation of MDF in kJ
mdfParam.bottlenecks = 0; % change to 1 to compute thermodynamic bottlenecks

%% preparing parameters for: P10 reactions that are active in core model
lbCore = zeros(cnap.numr,1);
lbCore(coreReacs) = cnap.reacMin(coreReacs);
ubCore = zeros(cnap.numr,1);
ubCore(coreReacs) = cnap.reacMax(coreReacs);
ubCore(2584:end) = 0;
lbCore(2584:end) = 0;
%% start characterization and ranking function
[IS_rankingStruct, IS_rankingTable] = characterizeIS(cnap,IS_lb,IS_ub,IS_to_rank,idx,cytMet,D,d,T,t,mdfParam,lbCore,ubCore,0:10,weightCoeff,solver);
%% save results into original .mat file
save([path filename '_ranking.mat'],'IS_rankingStruct', 'IS_rankingTable');
%% save ranking table to xls file
cell2csv([path filename '_ranking.xls'],IS_rankingTable,char(9));
