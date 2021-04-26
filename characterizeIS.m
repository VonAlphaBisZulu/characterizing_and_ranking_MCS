function [IS_rankingStruct, IS_rankingTable]...
= characterizeIS( cnapRef , IS_lb, IS_ub, IS_to_rank, idx, cytMet, D, d, T, t, mdfParam, lbCore, ubCore, crit, wCoeff,solver)
% CellNetAnalyzer API function 'characterizeIS'
% -----------------------------------------------
% Characterize a set of intervention strategies (IS) among different criteria, 
% group it into equivalence classes and generate a ranking.
% The characterization criteria are:
% 1) number of interventions, 2) max growth rate, 3) minimum Yield at maximum growth rate, 
% 4) minimum Yield, 5) oxygen requirement, 6) number of alternative Products, 
% 7) number of accessible metabolites, 8) optMDF, 9) overlap score, 
% 10) feasibility in core model
% The criteria to be assessed can be selected beforehand. The equivalence
% class grouping and the ranking is also optional.
% 
% Usage:    [IS_rankingStruct, IS_rankingTable] = characterizeIS( cnapRef , IS_lb, IS_ub, IS_to_rank, idx, cytMet, D, d, T, t, mdfParam, lbCore, ubCore, crit, wCoeff)
%
% Input:
%  =========== mandatory ===========
%    cnap                               : CellNetAnalyzer project structure containing
%                                          a stoichiometric model of the WT host
%    IS_lb, IS_ub (cell<n double<r>>)   : Upper and lower reaction bounds of each intervention
%           (n: number of IS)              strategy mutant. Flux bounds are r-dimentioal (double)
%           (r:number of reactions)        vectors, nested in cell arrays, each cell representing 
%                                          one of the n mutant strains / intervention stragies.
%    IS_to_rank  (double<m>)            : Index vector that specifies which mutants should be
%        (m: number of IS to be ranked)    included in the ranking. Contains integers between 
%                (m<=n)                    1 and n.
%    idx (struct)                       : Contains a set of importent reaction and metabolite indices
%        .subs (double)                 :  index of substrate source reaction
%        .prod (double)                 :  index of product export reaction
%        .atpm (double)                 :  index of ATP Maintanance reaction
%        .o2 (double)                   :  index of o2 uptake reaction
%        .bm (double)                   :  index of biomass growth reaction
%     ---- optional ----
%        .atp_c , .adp_c , .amp_c, 
%        .pi_c , .h2o_c , .h_c , .nad_c,:  indices of important cytosolic metabolites
%        .nadh_c , .nadp_c , .nadph_c
%
%  =========== optional ===========                                      
%   cytMet     (double<q>)              : Indices of all cytosolic metabolites
%    (q: number of cytosolic metabolites)
%   D, d   (double<k,r>, double<k,1>)   : Contains Desired flux vector
%   T, t   (double<l,r>, double<l,1>)   : Contains Target flux vector
%   mdfParam  (struct)                  : Hold the necessary parameters to do the OptMDF computation
%       .G0         (double<r>)         :  contains the Gibbs free energy of each reaction
%                                           (NaN for the reaction with non-specified Gibbs energies)
%       .uncert     (double<r>)         :  contains uncertainties of the Gibbs free energies
%       .Cmin       (double<s>)         :  minimum metabolite concentration
%       .Cmax       (double<s>)         :  maximum metabolite concentration
%           (s: number of species)
%       .fixed_ratios   (double<b,3>)   :  specifies fixed ratios of metabolites. First two columns
%                                           contain species indices (between 1 and s), third column
%                                           represents the ratio between the two
%       .RT                             :  Gas constant (R)  x  Temperature (T)
%       .bottlenecks: (logical)         :  switch, if thermodynamic bottlenecks should be computed
%   lbCore, ubCore   (double<r>)        : Boundaries to constrain the original model to a core model
%   crit             (double<c>)        : Specifies what criteria shall be assessed 
%     (c: number of assessed criteria)     (between 0 and 10, whereas criterion 0 is the equivalence 
%                                           class grouping)
%   wCoeff           (double<c>)        : Contains weighting factors for each criteria, used to
%                                         generate the final score
%   solver           (double)           : Solver (0: GLPK, 2: CPLEX)

if ~exist('solver','var')
    solver = 0;
    disp('No solver specified. Using GLPK solver.');
else
    switch solver
        case 0 % GLPK
            disp('using GLPK solver.');
        case 1 % linprog
            disp('using linprog solver.');
        case 2 % CPLEX
            disp('using CPLEX solver.');
        otherwise
            error('choose valid solver (GLPK: 0, CPLEX: 2).');
    end
end
% add reference bounds to set of strategies to analyze
IS_rankingStruct = repmat(struct,length(IS_to_rank),1); % return structure
IS_to_rank = IS_to_rank+1;
IS_lb = [{cnapRef.reacMin}, IS_lb];
IS_ub = [{cnapRef.reacMax}, IS_ub];
disp([num2str(length(IS_lb)-length(IS_to_rank)) ' reference model(s) and ' num2str(length(IS_to_rank)) ' strain design(s)']);
% init variables (double)
[   p0_eqClass,     p1_numInt,      p2_maxGrowth, ...   
    p3_minYatMaxGr, p4_minY,        p5_O2req, ...      
    p6_numAltProd,  p7_accessMet,   p8_optMDF, ...
    p9_overlap,     p10_feasCore,   numcutsCore, ...
    rCoreBMmax,     YCorePAtBMmax ] = deal(nan(length(IS_lb),1));
% init variables (cell)
[   itv, intervs, vMDFlimR, lb_fva, ub_fva ] = deal(repmat({nan},length(IS_lb),1));
% init variables (struct)
AuxReac = repmat(struct('ATP_coupled',nan,...
                        'NAD_red_coupled',nan,...
                        'NADP_red_coupled',nan,...
                        'NADH_ox_coupled',nan,...
                        'NADPH_ox_coupled',nan),length(IS_lb),1);

if ~exist('crit','var')
    crit = [0,1,9];
    if exist('idx','var')
        if ~isempty(idx)
            crit = [crit,2,3,4,5];
        end
    end
    if exist('cytMet','var')
        if ~isempty(cytMet)
            if exist('T','var') && exist('t','var')
                if ~isempty(T) && ~isempty(t)
                    crit = sort([crit,6]);
                end
            end
            crit = sort([crit,7]);
        end
    end
    if exist('mdfParam','var')
        if ~isempty(mdfParam) && solver == 2
            if all(arrayfun(@(x) isfield(mdfParam,x),{'G0', 'uncert','Cmin','Cmax','fixed_ratios','RT'}))
                crit = sort([crit,8]);
            end
        end
    end
    if exist('lbCore','var') && exist('ubCore','var') && exist('D','var') && exist('d','var') && exist('T','var') && exist('t','var')
        if ~isempty(lbCore) && ~isempty(ubCore) && ~isempty(D) && ~isempty(d) && ~isempty(T) && ~isempty(t)
            crit = [crit,10];
        end
    end
    crit = sort(crit);
end


% weighting
if exist('wCoeff','var')
    if ~isempty(wCoeff)
        wCoeff = full(sparse(crit(crit~=0),1,wCoeff,10,1));
    else
        wCoeff = zeros(10,1);
    end
else
    wCoeff = zeros(10,1);
end

% prepare criterion 8 - MDF
if ismember(8,crit)
    if isfield(mdfParam,'bottlenecks')
        bottlenecks = mdfParam.bottlenecks;
    else
        bottlenecks = 0;
    end
    jpath = javaclasspath('-dynamic');; % to be able to find CPLEX-java-library
    if solver ~= 2
        warning('OptMDF analysis needs CPLEX solver to be functional. If characterizeIS fails, consider to either exclude criterion 8 or install and configure CPLEX.');
    end
else, bottlenecks=[]; jpath = {''};
end

% prepare FVA / grouping MCS in families / calcuate zero rates to speed up future calculations
if ismember(0,crit)
    disp('preparational FVA for equivalence class determination...')
    [unknown_fluxes,zero_fluxes] = prep0_FVA_eq(cnapRef,min(cell2mat(IS_lb),[],2),max(cell2mat(IS_ub),[],2),solver);
else, unknown_fluxes=[]; zero_fluxes=[];
end

% prepare criterion 10 - reduced model - if existing
if ismember(10,crit)
    [cnapCore, Dfeas, Tfeas] = prep10_redModel(cnapRef,lbCore,ubCore,D,d,T,t,solver);
    if ~Dfeas || ~Tfeas
        % losen constraints to those that can occor in the mutants through
        % additional reactions
        rWiderBounds = any([(cell2mat(IS_lb)' < cnapRef.reacMin');(cell2mat(IS_ub)' > cnapRef.reacMax')]);
        lbc = lbCore;
        ubc = ubCore;
        lbmin = min(cell2mat(IS_lb),[],2);
        ubmax = max(cell2mat(IS_ub),[],2);
        lbc(rWiderBounds) = lbmin(rWiderBounds);
        ubc(rWiderBounds) = ubmax(rWiderBounds);
        warning('Desired or Target mode not feasible in wild type core model, checking with broader constraints');
        [~, Dfeas, Tfeas] = prep10_redModel(cnapRef,lbc,ubc,D,d,T,t,solver);
        if ~Dfeas || ~Tfeas
            error('Desired or Target mode still infeasible in core model');
            % none of the strain designs will work in the constrained model
        end
    end
else
    crit = crit(crit~=10);
    cnapCore = struct.empty();
end

% Criterion 1: Number of cuts
if ismember(1,crit),       p1_numInt = c1_numCuts(cnapRef.reacMin,cnapRef.reacMax,IS_lb,IS_ub,idx);
end

% Criterion 5: Is anaerobic
if ismember(5,crit),       p5_O2req = c5_isAnaerobic(IS_lb,idx);
end

% Criterion 9: Overlapping score
if ismember(9,crit),       p9_overlap = c9_overlapScore(cnapRef.reacMin,cnapRef.reacMax,IS_lb,IS_ub,idx,IS_to_rank);
end

%% use parallel pool or not? Approximation. Depends on size of the problem and number of iterations.
if license('test','Distrib_Computing_Toolbox') && ~isempty(ver('distcomp'))  && isempty(getCurrentTask()) % also check if on main or worker thread
    parforArg = inf;
    if length(IS_lb) > 1
        if ~verLessThan('matlab', '8.1') && isempty(gcp('nocreate'))
                parpool();
        elseif(verLessThan('matlab', '8.1') && ~matlabpool('size'))
                matlabpool(); %#ok<DPOOL> % suppress compatibility warning
        end
    else    % smaller problems
        parforArg = 0;
    end
    % suppress name warning
    if ~verLessThan('matlab', '8.1')
    	parfevalOnAll(@() warning('off', 'MATLAB:mir_nr_id_will_not_be_shared_variable'),1);
    end
    if ismember(8,crit) || solver == 2
        %parfevalOnAll(@() javaaddpath(jpath(~ismember(jpath,javaclasspath))),1); % add missing java libraries to parallel workspaces
    end
else
    parforArg = 0;
end

%% loop through the intervention strategies
rowNames = repmat({''},length(IS_lb),1);
disp('Parallel pool is used. Order of computation may vary, order in final output is conserved.');
parfor (i = 1:length(IS_lb), parforArg)
% for i = 1:length(IS_lb)
    if ismember(i,IS_to_rank)
        rowNames{i} = ['IS_' num2str(find(ismember(IS_to_rank,i)),['%0' num2str(ceil(log10(length(IS_to_rank)+1))),'d'])]; % leading zeros
        disp(['characterizing strain design ' num2str(find(ismember(IS_to_rank,i))) ' out of ' num2str(length(IS_to_rank)) ' (please ignore the following output)']);
    else
        rowNames{i} = ['ref_' num2str(num2str(find(ismember(find(~ismember(1:length(IS_lb),IS_to_rank)),i))),['%0' num2str(ceil(log10(sum(~ismember(1:length(IS_lb),IS_to_rank))))),'d'])]; % leading zeros
        disp(['characterizing reference model ' num2str(find(ismember(find(~ismember(1:length(IS_lb),IS_to_rank)),i))) ' out of ' ...
              num2str(length(IS_lb)-length(IS_to_rank)) ' (please ignore the following output)']);
    end
    % get interventions
    itv{i} = [-find(or(and(cnapRef.reacMin ~= IS_lb{i}, cnapRef.reacMin ~= 0) , and(cnapRef.reacMax ~= IS_ub{i}, cnapRef.reacMax ~= 0))); ... knock-outs or regulatory interventions
            find(or(and(cnapRef.reacMin ~= IS_lb{i}, cnapRef.reacMin == 0) , and(cnapRef.reacMax ~= IS_ub{i}, cnapRef.reacMax == 0)))]; % new reactions
    intervs(i) = cellstr(strtrim(strjoin(cellstr(num2str(itv{i},'%+d')))));
    
    cnap = cnapRef;
    cnap.reacMin = IS_lb{i};
    cnap.reacMax = IS_ub{i};
    cnap.epsilon = 1e-7;
    fv = nan(cnapRef.numr,1);
    
    % Get upper and lower bounds in FVA to identify equivalence classes later
    if ismember(0,crit),     [lb_fva{i}, ub_fva{i}] = c0_getFVA(cnap,unknown_fluxes,zero_fluxes,solver); end
    % Criterion 1: Number of cuts
        % calculated before
    % Criterion 2: Maximum growth rate
    if any(ismember([2,3],crit)),     [fv,p2_maxGrowth(i,1)]   = c2_maxGrowthrate(cnap,idx,solver);  end
    % Criterion 3: Minimum product yield at maximum growth rate
    if ismember(3,crit),                  p3_minYatMaxGr(i,1) = c3_minYatMaxMue(cnap,fv,idx,solver);  end
    % Criterion 4: Minimum product yield
    if ismember(4,crit),                  p4_minY(i,1)        = c4_minProdY(cnap,idx,solver);  end
    % Criterion 5: Dependency on oxygen
        % determined before
    % Criterion 6: Alternative products that could disrupt the coupling strategy
    if ismember(6,crit)
        [IdxProdCandidates, lpModel] = prep67_cpxModel(cnap,cytMet,solver);
        [p6_numAltProd(i,1),AuxReac(i,1)] = c6_disruptiveMetabolites(cnap,lpModel,IdxProdCandidates,T,t,idx,solver);
    end
    % Criterion 7: Metabolites that are still accessible
    if ismember(7,crit)
        [IdxProdCandidates, lpModel] = prep67_cpxModel(cnap,cytMet,solver);
        p7_accessMet(i,1) = c7_numMetAccessible(cnap,lpModel,IdxProdCandidates,solver);
    end
    % Criterion 8: MDF
    if ismember(8,crit)
        javaaddpath(jpath(~ismember(jpath,javaclasspath))); % add missing java libraries to parallel workspaces
        [p8_optMDF(i,1), vMDFlimR(i,1)] = c8_computeMDF(cnap, mdfParam.RT, ...
                                                              mdfParam.G0, ...
                                                              mdfParam.uncert, ...
                                                              mdfParam.Cmin, ...
                                                              mdfParam.Cmax, ... 
                                                              mdfParam.fixed_ratios, ...
                                                              D, d, bottlenecks );
    end
    % Criterion 9: Overlapping
    % determined before
    % Criterion 10: Reduced model
    if ismember(10,crit)
        [numcutsCore(i,1),  rCoreBMmax(i,1),  YCorePAtBMmax(i,1)] = ...
            c10_checkInCoreModel(cnapCore, lbCore, ubCore, IS_lb{i}, IS_ub{i}, itv{i}, idx, solver);
        p10_feasCore(i,1) = rCoreBMmax(i,1)~=0;
    end
end
cnap = struct.empty(); % necessary to make parfor-loop work
%% Write output table
if length(IS_lb)>1   
    IS_rankingTable = [[cellstr(['production reaction: ' num2str(idx.prod)]);cellstr(rowNames)], [{'Interventions'}; intervs]];
    
    %% property values
    IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'itv',y),1:length(IS_to_rank),itv(IS_to_rank)');
    if ismember(0,crit)
        p0_eqClass  = c0_getEC(lb_fva,ub_fva,1e-5);
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'eqClass',y),1:length(IS_to_rank),p0_eqClass(IS_to_rank)');
        IS_rankingTable = [IS_rankingTable, [{'equivalence_class'}; num2cell(p0_eqClass)]];
    end
    if ismember(1,crit)
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'p1_numInt',y),1:length(IS_to_rank),p1_numInt(IS_to_rank)');
        IS_rankingTable = [IS_rankingTable, [{'P1_num_Intv'}; num2cell(p1_numInt)]]; 
    end
    if ismember(2,crit)
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'p2_maxGrowth',y),1:length(IS_to_rank),p2_maxGrowth(IS_to_rank)');
        IS_rankingTable = [IS_rankingTable, [{'P2_max_growth'};num2cell(roundDec(p2_maxGrowth,3))]]; 
    end
    if ismember(3,crit)
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'p3_minYatMaxGr',y),1:length(IS_to_rank),p3_minYatMaxGr(IS_to_rank)');
        IS_rankingTable = [IS_rankingTable, [{'P3_min_Y_at_max_growth'};num2cell(roundDec(p3_minYatMaxGr,3))]]; 
    end
    if ismember(4,crit)
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'p4_minY',y),1:length(IS_to_rank),p4_minY(IS_to_rank)');
        IS_rankingTable = [IS_rankingTable, [{'P4_min_Y'};num2cell(roundDec(p4_minY,3))]];
    end    
    if ismember(5,crit)
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'p5_O2req',y),1:length(IS_to_rank),p5_O2req(IS_to_rank)');
        IS_rankingTable = [IS_rankingTable, [{'P5_anaerobic'}; num2cell(p5_O2req)]];
    end
    if ismember(6,crit)
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'p6_numAltProd',y),1:length(IS_to_rank),p6_numAltProd(IS_to_rank)');
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'ATP_coupled',y),1:length(IS_to_rank),[AuxReac(IS_to_rank).ATP_coupled]);
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'NAD_red_coupled',y),1:length(IS_to_rank),[AuxReac(IS_to_rank).NAD_red_coupled]);
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'NADP_red_coupled',y),1:length(IS_to_rank),[AuxReac(IS_to_rank).NADP_red_coupled]);
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'NADH_ox_coupled',y),1:length(IS_to_rank),[AuxReac(IS_to_rank).NADH_ox_coupled]);
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'NADPH_ox_coupled',y),1:length(IS_to_rank),[AuxReac(IS_to_rank).NADPH_ox_coupled]);
        IS_rankingTable = [IS_rankingTable, [{'P6_num_coupling_disrupting_metabolites'}; num2cell(p6_numAltProd)], ...
            [{'ATP_coupled'};{AuxReac(:).ATP_coupled}'], ...
            [{'NAD_reduction_disrupts_coupling'};{AuxReac(:).NAD_red_coupled}'], ...
            [{'NADP_reduction_disrupts_coupling'};{AuxReac(:).NADP_red_coupled}']...
            [{'NADH_oxidation_disrupts_coupling'};{AuxReac(:).NADH_ox_coupled}'], ...
            [{'NADPH_oxidation_disrupts_coupling'};{AuxReac(:).NADPH_ox_coupled}']];
    end  
    if ismember(7,crit)
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'p7_accessMet',y),1:length(IS_to_rank),p7_accessMet(IS_to_rank)');
        IS_rankingTable = [IS_rankingTable, [{'P7_num_accessible_metabolites'}; num2cell(p7_accessMet)]];
    end
    if ismember(8,crit)
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'p8_optMDF',y),1:length(IS_to_rank),p8_optMDF(IS_to_rank)');
        IS_rankingTable = [IS_rankingTable, [{'P8_optMDF'}; num2cell(p8_optMDF,3)]];
        if bottlenecks
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'vMDFlimR',y),1:length(IS_to_rank),vMDFlimR(IS_to_rank)');
            IS_rankingTable = [IS_rankingTable, [{'vMDFlimR'};  vMDFlimR]];
        end
    end
    if ismember(9,crit)
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'p9_overlap',y),1:length(IS_to_rank),p9_overlap(IS_to_rank)');
        IS_rankingTable = [IS_rankingTable, [{'P9_overlap_Score'}; num2cell(roundDec(p9_overlap,3))]];
    end
    if ismember(10,crit)
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'p10_feasCore',y),1:length(IS_to_rank),p10_feasCore(IS_to_rank)');
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'numcutsCore',y),1:length(IS_to_rank),numcutsCore(IS_to_rank)');
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'rCoreBMmax',y),1:length(IS_to_rank),rCoreBMmax(IS_to_rank)');
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'YCorePAtBMmax',y),1:length(IS_to_rank),YCorePAtBMmax(IS_to_rank)');
        IS_rankingTable = [IS_rankingTable, [[{'P10_feas_core'};num2cell(roundDec(p10_feasCore,3))], ...
            [{'num_intvs_core'}; num2cell(numcutsCore,3)], ...
            [{'rBM_max_core'};num2cell(roundDec(rCoreBMmax,3))], ...
            [{'Ymin_rBM_max_core'};num2cell(roundDec(YCorePAtBMmax,3))]]];
    end
    %% ranking
    if any(wCoeff) && length(IS_to_rank)>1
        % initialize
        [s1_numInt, s2_maxGrowth, s3_minYatMaxGr, s4_minY, s5_O2req, s6_numAltProd, s7_accessMet, s8_optMDF, s9_overlap, s10_feasCore] = deal(zeros(length(IS_lb),1));
        % individual scores (only if properties have a range greater than zero)
        if ismember(1,crit) && range(p1_numInt(IS_to_rank)),        s1_numInt(IS_to_rank)      = (max(p1_numInt(IS_to_rank))-p1_numInt(IS_to_rank))            /(max(p1_numInt(IS_to_rank))-min(p1_numInt(IS_to_rank))); end
        if ismember(2,crit) && range(p2_maxGrowth(IS_to_rank)),     s2_maxGrowth(IS_to_rank)   = (p2_maxGrowth(IS_to_rank)-min(p2_maxGrowth(IS_to_rank)))      /(max(p2_maxGrowth(IS_to_rank))-min(p2_maxGrowth(IS_to_rank))); end
        if ismember(3,crit) && range(p3_minYatMaxGr(IS_to_rank)),   s3_minYatMaxGr(IS_to_rank) = (p3_minYatMaxGr(IS_to_rank)-min(p3_minYatMaxGr(IS_to_rank)))  /(max(p3_minYatMaxGr(IS_to_rank))-min(p3_minYatMaxGr(IS_to_rank))); end
        if ismember(4,crit) && range(p4_minY(IS_to_rank)),          s4_minY(IS_to_rank)        = (p4_minY(IS_to_rank)-min(p4_minY(IS_to_rank)))                /(max(p4_minY(IS_to_rank))-min(p4_minY(IS_to_rank))); end
        if ismember(5,crit) && range(p5_O2req(IS_to_rank)),         s5_O2req(IS_to_rank)       = p5_O2req(IS_to_rank); end
        if ismember(6,crit) && range(p6_numAltProd(IS_to_rank)),    s6_numAltProd(IS_to_rank)  = (max(p6_numAltProd(IS_to_rank))-p6_numAltProd(IS_to_rank))    /(max(p6_numAltProd(IS_to_rank))-min(p6_numAltProd(IS_to_rank))); end
        if ismember(7,crit) && range(p7_accessMet(IS_to_rank)),     s7_accessMet(IS_to_rank)   = (max(p7_accessMet(IS_to_rank))-p7_accessMet(IS_to_rank))      /(max(p7_accessMet(IS_to_rank))-min(p7_accessMet(IS_to_rank))); end
        if ismember(8,crit) && range(p8_optMDF(IS_to_rank)),        s8_optMDF(IS_to_rank)      = max((p8_optMDF(IS_to_rank)-max(0,min(p8_optMDF(IS_to_rank)))) /(max(p8_optMDF(IS_to_rank))-max(0,min(p8_optMDF(IS_to_rank)))),0); end
        if ismember(9,crit) && range(p9_overlap(IS_to_rank)),       s9_overlap(IS_to_rank)     = (p9_overlap(IS_to_rank)-min(p9_overlap(IS_to_rank)))          /(max(p9_overlap(IS_to_rank))-min(p9_overlap(IS_to_rank))); end
        if ismember(10,crit)&& range(p10_feasCore(IS_to_rank)),     s10_feasCore(IS_to_rank)   = p10_feasCore(IS_to_rank); end
        totalScore = sum([  wCoeff( 1)*s1_numInt, ...
                            wCoeff( 2)*s2_maxGrowth, ...
                            wCoeff( 3)*s3_minYatMaxGr, ...
                            wCoeff( 4)*s4_minY, ...
                            wCoeff( 5)*s5_O2req, ...
                            wCoeff( 6)*s6_numAltProd, ...
                            wCoeff( 7)*s7_accessMet, ...
                            wCoeff( 8)*s8_optMDF, ...
                            wCoeff( 9)*s9_overlap, ...
                            wCoeff(10)*s10_feasCore],2);
        repstv = logical(iv(length(IS_lb),IS_to_rank));
        rank = zeros(length(IS_lb),1);
        if ismember(0,crit) && ismember(7,crit) %% Select representatives and rank
            for c = unique(p0_eqClass(IS_to_rank))'
                repstv(p0_eqClass == c) = 0;
                fewest_accMet = (p7_accessMet== min(p7_accessMet(p0_eqClass == c))) & p0_eqClass == c;
                highest_score = totalScore.*fewest_accMet == max(totalScore.*fewest_accMet);
                repstv(find(highest_score,1)) = 1;
            end
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'is_representative',y),1:length(IS_to_rank),repstv(IS_to_rank)');
        end
        [score,pos] = sort(repstv.*totalScore,'descend');
        rank(pos) = [1:sum(repstv) zeros(1,sum(repstv==0))];
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'rank',y),1:length(IS_to_rank),rank(IS_to_rank)');
        % write in struct and table
        if ismember(1,crit)
            IS_rankingTable = [IS_rankingTable, [{'S1'}; num2cell(s1_numInt)]]; 
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'S1',y),1:length(IS_to_rank),s1_numInt(IS_to_rank)'); 
        end
        if ismember(2,crit)
            IS_rankingTable = [IS_rankingTable, [{'S2'}; num2cell(roundDec(s2_maxGrowth,3))]]; 
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'S2' ,y),1:length(IS_to_rank),s2_maxGrowth(IS_to_rank)'); 
        end
        if ismember(3,crit)
            IS_rankingTable = [IS_rankingTable, [{'S3'}; num2cell(roundDec(s3_minYatMaxGr,3))]]; 
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'S3' ,y),1:length(IS_to_rank),s3_minYatMaxGr(IS_to_rank)'); 
        end
        if ismember(4,crit)
            IS_rankingTable = [IS_rankingTable, [{'S4'}; num2cell(roundDec(s4_minY,3))]]; 
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'S4' ,y),1:length(IS_to_rank),s4_minY(IS_to_rank)'); 
        end
        if ismember(5,crit)
            IS_rankingTable = [IS_rankingTable, [{'S5'}; num2cell(s5_O2req)]]; 
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'S5' ,y),1:length(IS_to_rank),s5_O2req(IS_to_rank)'); 
        end
        if ismember(6,crit)
            IS_rankingTable = [IS_rankingTable, [{'S6'}; num2cell(s6_numAltProd)]]; 
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'S6' ,y),1:length(IS_to_rank),s6_numAltProd(IS_to_rank)'); 
        end
        if ismember(7,crit)
            IS_rankingTable = [IS_rankingTable, [{'S7'}; num2cell(s7_accessMet)]]; 
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'S7' ,y),1:length(IS_to_rank),s7_accessMet(IS_to_rank)'); 
        end
        if ismember(8,crit)
            IS_rankingTable = [IS_rankingTable, [{'S8'}; num2cell(roundDec(s8_optMDF,3))]]; 
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'S8' ,y),1:length(IS_to_rank),s8_optMDF(IS_to_rank)'); 
        end
        if ismember(9,crit)
            IS_rankingTable = [IS_rankingTable, [{'S9'}; num2cell(roundDec(s9_overlap,3))]]; 
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'S9' ,y),1:length(IS_to_rank),s9_overlap(IS_to_rank)'); 
        end
        if ismember(10,crit)
            IS_rankingTable = [IS_rankingTable, [{'S10'}; num2cell(s10_feasCore)]]; 
            IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'S10',y),1:length(IS_to_rank),s10_feasCore(IS_to_rank)'); 
        end
        IS_rankingTable = [IS_rankingTable, [{'Total'}; num2cell(roundDec(totalScore,3))]];
        IS_rankingStruct = arrayfun(@(x,y) setfield(IS_rankingStruct(x),'totalScore',y),1:length(IS_to_rank),totalScore(IS_to_rank)'); 
        if ismember(0,crit) && ismember(7,crit)
            IS_rankingTable = [IS_rankingTable, [{'is_representative'}; num2cell(repstv)]];
        end
        IS_rankingTable = [IS_rankingTable, [{'rank'};num2cell(rank)]];
    end
    % prepare output as table (return evalTable instead of evalCell if output as table is prefered)
%     C = IS_rankingTable(2:end,2:end);
%     IS_rankingTableTb = cell2table(C);
%     IS_rankingTableTb.Properties.VariableNames = IS_rankingTable(1,2:end);
%     IS_rankingTableTb.Properties.RowNames      = IS_rankingTable(2:end,1);
else
    warning(['no intervention strategy was passed ' char(cnap.reacID(idx.prod,:))]);
    IS_rankingTable = cell.empty(0,0);
end
end

function [unknown_fluxes,zero_fluxes] = prep0_FVA_eq(cnap,lb_min,ub_max,solver)
%% preparational - equivalence classes: get fixed rates to speed up future FVA computations
% returns vector
% Initial flux variability analysis to speed up identifying equivalence classes later
unknown_fluxes = find(lb_min ~= ub_max);

cnap.reacMin   = lb_min;
cnap.reacMax   = ub_max;
cnap.epsilon   = 1e-6;

[lb,ub] = CNAfluxVariability(cnap,[],[],solver,unknown_fluxes);
zero_fluxes = unknown_fluxes(lb==ub & ub==0);
unknown_fluxes = unknown_fluxes(lb~=ub);
end

function [lb_fva, ub_fva] = c0_getFVA(cnap,unknown_fluxes,zero_fluxes,solver)
%% IS equivalent classes: get ub and lb vectors of the fva for each MCS
% single value
cnap.reacMin(zero_fluxes) = 0;
cnap.reacMax(zero_fluxes) = 0;
[lb_fva,ub_fva] = CNAfluxVariability(cnap,[],[],solver,unknown_fluxes);
end

function eqClass = c0_getEC(lb_fva,ub_fva,epsilon)
%% IS equivalent classes: get ub and lb vectors of the fva for each MCS
% single value
eqClassMat = full(sparse(1,1,1,1,length(lb_fva)));
for i = 2:length(lb_fva)
    for j = 1:size(eqClassMat,1)
        if all(abs(lb_fva{i} - lb_fva{find(eqClassMat(j,:),1)}) < epsilon) && all(abs(ub_fva{i} - ub_fva{find(eqClassMat(j,:),1)}) < epsilon)
            eqClassMat(j,i) = 1;
        end
    end
    switch sum(eqClassMat(:,i))
        case 1
        case 0
            eqClassMat(j+1,i) = 1;
        otherwise
            warning('ambiguos equivalence class association');
    end
end
[eqClass,~] = find(eqClassMat);
end

function numIntv = c1_numCuts(lbref,ubref,lb,ub,idx)
%% criterion 1: number of interventions
% returns vector with the number of necessary interventions per strategy
% interventions are defined as number of boundaries that are changed from
% original reference
numIntv = zeros(length(lb),1);
for i = 1:length(lb)
    % number of interventions w/o o2
    numIntv(i) = length(setdiff([find(ubref ~= ub{i}) ; find(lbref ~= lb{i})],idx.o2));
end
end

function [fv,maxGrowthRate] = c2_maxGrowthrate(cnap, idx,solver)
%% criterion 2: maximum possible growth rate
% returns single value
cnap.objFunc = -iv(cnap.numr,idx.bm);
[fv,~,~,maxGrowthRate] = CNAoptimizeFlux(cnap,[],[],solver);
maxGrowthRate = -maxGrowthRate;
end

function minProdYAtMueMax  = c3_minYatMaxMue(cnap,fv,idx,solver)
%% criterion 3: minimum product yield at maximum growth
% returns single value
if solver == 1 % 1e-6 is used as tolerance for the priorily calculated
               % growth rate. Otherwise the optimize Yield function fails with
    cnap.epsilon = 1e-6;
end
minProdYAtMueMax = -CNAoptimizeYield(cnap,-iv(cnap.numr,idx.prod)',-iv(cnap.numr,idx.subs)',fv(idx.bm)*iv(cnap.numr,idx.bm,nan)-cnap.epsilon*10,[],solver);
end

function minY = c4_minProdY(cnap,idx,solver)
%% criterion 4: minimum product yield
% returns single value
cnap.objFunc = iv(cnap.numr,idx.prod);
minY = -CNAoptimizeYield(cnap,-iv(cnap.numr,idx.prod)',-iv(cnap.numr,idx.subs)',[],[],solver);
end

function isAnaerobic = c5_isAnaerobic(lb,idx)
%% criterion 5: is mcs anaerobic?
% returns column vector
isAnaerobic = cellfun(@(x) x(idx.o2)>=0,lb)';
end

function [IdxProdCandSinkReacs, lpModel] = prep67_cpxModel(cnap, cytMet, solver)
%% preparation for criterion 6 and/or 7
% one run to get possibly reachable cytoplasmic metabolites
% prebuild Cplex model
% - this preselection speeds up later computations

reacsCytEx = full(-sparse(cytMet,1:length(cytMet),1,cnap.numis,length(cytMet))); % export reactions for cyt Metabolites

lpModel.N = initsmat(cnap.stoichMat,cnap.mue,cnap.macroComposition,cnap.local.c_makro,cnap.specInternal);
lpModel.N(:,cnap.numr+1:cnap.numr+length(cytMet)) = reacsCytEx;
N = -lpModel.N; % Needed to turn constraint into Nr >= 0
lpModel.b = zeros(1,cnap.numis);
lpModel.lb = zeros(cnap.numr+length(cytMet),1);
lpModel.ub = [zeros(cnap.numr,1)    ; ones(length(cytMet),1)];
lpModel.lb(cnap.reacMin<0) = -Inf;
lpModel.lb(cnap.reacMin>0) = cnap.reacMin(cnap.reacMin>0);
lpModel.ub(cnap.reacMax>0) = Inf;
lpModel.lb(cnap.reacMax<0) = cnap.reacMax(cnap.reacMax<0);
lpModel.obj = -[zeros(cnap.numr,1); ones(length(cytMet),1)];
lpModel.glpk_ctype = repmat('U',cnap.numis,1); % only bound by rhs
lpModel.glpk_options.round = 1;
if ~verLessThan('matlab', '8.1') % use different functions to suppress linprog verbose depending on MATLAB version
    lpmodel.linprog_options = optimoptions('linprog','Display','off');
else
    lpmodel.linprog_options = optimset('Display','off');
end
lpModel.vartype = repmat('C',size(N,2),1);
switch solver
    case 0 % GLPK
        [x, ~,exitflag] =  glpk(lpModel.obj,N,lpModel.b,lpModel.lb,lpModel.ub,lpModel.glpk_ctype);
        % exitflag ==  2: optimal % exitflag ==  5: feasible
    case 1 % linprog
        [x, ~, exitflag] = linprog(lpModel.obj, N , lpModel.b , [], [], lpModel.lb, lpModel.ub,lpModel.linprog_options);
    case 2 % CPLEX
        [x, ~, exitflag] = cplexlp(lpModel.obj, N , lpModel.b , [], [], lpModel.lb, lpModel.ub);
        % exitflag ==  1: optimal, feasible  % exitflag == -2: infeasible
end
if (exitflag == 1 && (solver == 2 || solver == 1)) ... CPLEX & Linprog
        || ((exitflag == 2 || exitflag == 5) && solver == 0) % GLPK
    x(abs(x)<=1e-10) = 0;
else
    error('problem infeasible');
end
IdxProdCandSinkReacs = find(any(N(:,cnap.numr+find(x(cnap.numr+1:end))),2)); % find producible candidates
IdxReachableReac     = cnap.numr+find(x(cnap.numr+1:end)); % find producible candidates
lpModel.N         = lpModel.N(:, [1:cnap.numr IdxReachableReac']);
lpModel.lb        = [cnap.reacMin; zeros(length(IdxReachableReac),1)];
lpModel.ub        = [cnap.reacMax; zeros(length(IdxReachableReac),1)];
lpModel.obj       = [cnap.objFunc; zeros(length(IdxReachableReac),1)];
lpModel.vartype   = lpModel.vartype([1:cnap.numr IdxReachableReac']);
end

function [numDisruptMet,AuxReac] = c6_disruptiveMetabolites(cnap,lpModel,IdxProdCandidates,T,t,idx,solver)
%% criterion 6: check disruptive metabolites
numDisruptMet=zeros(1,length(IdxProdCandidates));
T = [T zeros(size(T,1),length(IdxProdCandidates))]; % append T with zeros

if ~verLessThan('matlab', '8.1') % use different functions to suppress linprog verbose depending on MATLAB version
    linprog_options = optimoptions('linprog','Display','off');
else
    linprog_options = optimset('Display','off');
end
glpk_ctype = [repmat('S',1,cnap.numis), repmat('U',1,length(t))];

for j = 1:length(IdxProdCandidates) % parfor
    exitflag = 0;
    % open export for a certain cytoplasmic metabolite
    ub = lpModel.ub;
    ub(cnap.numr+j) = Inf;
    % check if target modes become feasible
    switch solver
        case 0 % GLPK       % exitflag ==  2: optimal  % exitflag ==  5: feasible
            [~, ~,exitflag] =  glpk(lpModel.obj,[lpModel.N; T],[lpModel.b, t'],lpModel.lb,ub,glpk_ctype);
        case 1 % linprog    % exitflag ==  1: optimal, feasible
            [~, ~, exitflag] = linprog(lpModel.obj, T, t, lpModel.N, lpModel.b', lpModel.lb, ub, linprog_options);
        case 2 % CPLEX      % exitflag ==  1: optimal, feasible  % exitflag == -2: infeasible
            [~, ~, exitflag] = cplexlp(lpModel.obj, T, t, lpModel.N, lpModel.b', lpModel.lb, ub);
    end
    if (exitflag == 1 && (solver == 2 || solver == 1)) ... CPLEX & Linprog
            || ((exitflag == 2 || exitflag == 5) && solver == 0) % GLPK
        numDisruptMet(j) = 1; % Target mode became feasible again -> mark as disruptive
    end
end
numDisruptMet = sum(numDisruptMet);

%% check also if ATP recovery or NADPH or NADH recovery could
% disrupt coupling
numRows = size(lpModel.N,1);
numCols = size(lpModel.N,2);
% if provided, check for reduction equivalent regeneration
if all(arrayfun(@(x) isfield(idx,x),{'atp', 'h2o','adp','h','pi','nadh','nad','nadph','nadp'}))
    R_ATP_recover   = iv(numRows,[idx.atp idx.h2o]) - iv(numRows,[idx.adp idx.h idx.pi]);
    R_NAD_red       = iv(numRows,[idx.nadh idx.h])  - iv(numRows,idx.nad);
    R_NADP_red      = iv(numRows,[idx.nadph idx.h]) - iv(numRows,idx.nadp);
    R_NADH_ox       = iv(numRows,idx.nad)           - iv(numRows,[idx.nadh idx.h]);
    R_NADPH_ox      = iv(numRows,idx.nadp)          - iv(numRows,[idx.nadph idx.h]);

    lpModel.N = [lpModel.N, ...
        R_ATP_recover,...
        R_NAD_red,...
        R_NADP_red,...
        R_NADH_ox,...
        R_NADPH_ox];

    lpModel.lb  = [lpModel.lb;  zeros(5,1)];
    lpModel.ub  = [lpModel.ub;  zeros(5,1)];
    lpModel.obj   = [lpModel.obj ;  zeros(5,1)];
    T = [T zeros(size(T,1),5)];
end

% check feasabilitiy of target vectors if recovery reactions are added
% ATP
numCouplingMech=zeros(1,5);
for j = 1:5 % parfor
    exitflag = 0;
    ub = lpModel.ub;
    ub(numCols+j) = Inf;
    switch solver
        case 0 % GLPK       % exitflag ==  2: optimal  % exitflag ==  5: feasible
            [~, ~,exitflag] =  glpk(lpModel.obj,[lpModel.N; T],[lpModel.b, t'],lpModel.lb,ub,glpk_ctype);
        case 1 % linprog    % exitflag ==  1: optimal, feasible
            [~, ~, exitflag] = linprog(lpModel.obj, T, t, lpModel.N, lpModel.b', lpModel.lb, ub, linprog_options);
        case 2 % CPLEX      % exitflag ==  1: optimal, feasible  % exitflag == -2: infeasible
            [~, ~, exitflag] = cplexlp(lpModel.obj, T, t, lpModel.N, lpModel.b', lpModel.lb, ub);
    end
    if (exitflag == 1 && (solver == 2 || solver == 1)) ... CPLEX & Linprog
            || ((exitflag == 2 || exitflag == 5) && solver == 0) % GLPK
        numCouplingMech(j) = 1; % if feasible -> 1 -> coupling is disrupted / T is feasible again
    end
end
 AuxReac.ATP_coupled = numCouplingMech(1); % TRUE if problem is not 'infeasible'
 AuxReac.NAD_red_coupled = numCouplingMech(2);
 AuxReac.NADP_red_coupled = numCouplingMech(3);
 AuxReac.NADH_ox_coupled = numCouplingMech(4);
 AuxReac.NADPH_ox_coupled = numCouplingMech(5);
end

function numMetAccessible = c7_numMetAccessible(cnap,lpModel,IdxProdCandidates,solver)
%% criterion 7: Number of accessible metabolites
N = -lpModel.N; % To mime Nr >= 0 constraint (-Nr <= 0)
if ~verLessThan('matlab', '8.1') % use different functions to suppress linprog verbose depending on MATLAB version
    linprog_options = optimoptions('linprog','Display','off');
else
    linprog_options = optimset('Display','off');
end
glpk_ctype = repmat('U',1,cnap.numis); % -Nr <= 0 (that's why upper bound is used)
numMetAccessible=zeros(1,5);
for j = 1:length(IdxProdCandidates) % parfor
    exitflag = 0;
    r = zeros(cnap.numr,1);
    obj = zeros(size(N,2),1);
    ub  = lpModel.ub;
    obj(cnap.numr+j) = -1;
    ub(cnap.numr+j)  = 1;
    switch solver
        case 0 % GLPK       % exitflag ==  2: optimal  % exitflag ==  5: feasible
            [r, ~,exitflag] =  glpk(   obj, N, lpModel.b,          lpModel.lb, ub, glpk_ctype);
        case 1 % linprog    % exitflag ==  1: optimal, feasible
            [r, ~, exitflag] = linprog(obj, N, lpModel.b', [], [], lpModel.lb, ub, linprog_options);
        case 2 % CPLEX      % exitflag ==  1: optimal, feasible  % exitflag == -2: infeasible
            [r, ~, exitflag] = cplexlp(obj, N, lpModel.b', [], [], lpModel.lb, ub);
    end
    if (exitflag == 1 && (solver == 2 || solver == 1)) ... CPLEX & Linprog
            || ((exitflag == 2 || exitflag == 5) && solver == 0) % GLPK
        numMetAccessible(j) = r(cnap.numr+j)>1e-10;
    end
end
numMetAccessible = sum(numMetAccessible);
end

function [mdf, limReacID]  = c8_computeMDF(cnap, RT, G0, uncert, Cmin, Cmax, fixed_ratios, D, d, bottlenecks)
%% if MDF data are available: Test if Output reaction is in best flux vector
mdf= CNAcomputeOptMDFpathway(cnap, RT, [G0 uncert], Cmin, Cmax,[],[],fixed_ratios); % Add D and d for also verifying desired mode
% Find thermodynamic bottlenecks if required
    if bottlenecks
        mdf2 = nan(cnap.numr,1);
        for i = find(~isnan(G0))'
            G0_1    = G0;
            G0_1(i) = G0(i)-abs(G0(i)*0.01);
            mdf2(i) = CNAcomputeOptMDFpathway(cnap, RT, [G0_1 uncert], Cmin, Cmax,[],[],fixed_ratios);
        end
        limReacID = logical(mdf2>mdf);
        limReacID = {strjoin(cellstr(cnap.reacID(limReacID,:)),', ')};
    else
        limReacID = {''};
    end
end

function overlapScore = c9_overlapScore(lbref,ubref,lb,ub,idx,index_IS)
% overlapping score (count occurencies of cuts, sum them up for each mcs and
% devide by the number of cuts)
overlapScore = zeros(length(lb),1);
intv = (repmat(ubref,1,length(index_IS)) ~= [ub{index_IS}] | repmat(lbref,1,length(index_IS)) ~= [lb{index_IS}])';
ivWoO2 = intv(:,[1:idx.o2-1 idx.o2+1:end]);
[~, reacInd] = find(ivWoO2);
reacInd = reshape(reacInd,length(reacInd),1); % compatibility
[n,reacInd] = hist(reacInd,unique(reacInd));
intvWeighted = double(ivWoO2);
for r = reacInd'
    intvWeighted(ivWoO2(:,r)==1,r) = n(r==reacInd');
end
overlapScore(index_IS) = sum(intvWeighted,2)./sum(ivWoO2,2);
end

function [cnapCore,Dfeas,Tfeas] = prep10_redModel(cnap,lbCore,ubCore,D,d,T,t,solver)
cnapCore = cnap;
cnapCore.reacMin = lbCore;
cnapCore.reacMax = ubCore;
cnapCore.epsilon = 1e-7;
% Check if target and desired vectors are feasible in core model
Dfeas = CNAcheckDesired(cnapCore,D,d,solver);
Tfeas = CNAcheckDesired(cnapCore,T,t,solver);
end

function [numcutsCore,  rCoreBMmax,  YCorePAtBMmax] = ...
            c10_checkInCoreModel(cnapCore, lbCore, ubCore, lb, ub, itv, idx,solver)
% test in Core model with reactions added or deleted
interv = setdiff(abs(itv),idx.o2); % interventions (except for oxygen uptake)
numcutsCore = sum(any([lbCore(interv) ~= lb(interv) , ubCore(interv) ~= ub(interv)],2));
cnapCore.reacMin(abs(itv)) = lb(abs(itv));
cnapCore.reacMax(abs(itv)) = ub(abs(itv));
cnapCore.objFunc = -iv(cnapCore.numr, idx.bm);
fv = CNAoptimizeFlux(cnapCore,[],[],solver);
rCoreBMmax = fv(idx.bm);
YCorePAtBMmax = -CNAoptimizeYield(cnapCore,-iv(cnapCore.numr,idx.prod)',-iv(cnapCore.numr,idx.subs)',fv(idx.bm)*iv(cnapCore.numr,idx.bm,nan)-cnapCore.epsilon,[],solver);
end

% creates a column vector of *length, with ones at *indices
function ive = iv( len, indices,nan) % if third parameter is passed, all zeros are replaced by nan
% len: length of vector
% indices: indices of positions that should be 1
% row or column vector
ive = zeros(len,1);
ive(indices) = 1;
if nargin > 2
    ive(ive==0) = nan();
end
end

function f = roundDec( f , decimals )
switch nargin
    case 1
        f = round(f);
    case 2
        f = round(f*(10^(decimals)))/(10^(decimals));
    otherwise
        error('Invalid number of arguments');
end
end

function s = strjoin(terms, delimiter)
assert(iscellstr(terms), 'strjoin:invalidarg', ...
    'The first argument should be a cell array of strings.');

if nargin < 2
    d = ' ';
else
    d = delimiter;
    assert(ischar(d) && ndims(d) == 2 && size(d,1) <= 1, ...
        'strjoin:invalidarg', ...
        'The delimiter should be a char string.');
end
n = numel(terms);
if n == 0
    s = '';
elseif n == 1
    s = terms{1};
else
    ss = cell(1, 2*n-1);
    ss(1:2:end) = terms;
    [ss{2:2:end}] = deal(d);
    s = [ss{:}];
end
end

function [feasible, r] = CNAcheckDesired(cnap,D,d,solver)
if nargout >= 2 % do parsimonious FBA if two outputs are requested
    parsimonious = 1;
else
    parsimonious = 0;
end
if nargin < 4
    solver = 0;
end
if nargin == 1
    D = double.empty(0,cnap.numr);
    d = double.empty(0,1);
end
r   = zeros(cnap.numr,1);
N   = initsmat(cnap.stoichMat,cnap.mue,cnap.macroComposition,cnap.macroDefault,cnap.specInternal);
ub  = cnap.reacMax;
lb  = cnap.reacMin;
b   = zeros(cnap.numis, 1);
obj = zeros(cnap.numr, 1);
if parsimonious == 1  % Max rates from D vector (workaround for ratios - maximize rates)
    for i = 1:length(d)
        if d(i) ~= 0
            obj(D(i,:)~=0) = -D(i,D(i,:)~=0)*sign(d(i));
        else
            obj(D(i,:)~=0) = -D(i,D(i,:)>0);
            break;
        end
    end
end

switch solver
    case 0 % GLPK       % exitflag ==  2: optimal  % exitflag ==  5: feasible
        ctype = [repmat('S',1,cnap.numis), repmat('U',1,length(d))];
        [x, ~,exitflag] = glpk(obj,[N; D],[b; d],lb,ub,ctype);
    case 1 % linprog    % exitflag ==  1: optimal, feasible
	if ~verLessThan('matlab', '8.1') % use different functions to suppress linprog verbose depending on MATLAB version
    		linprog_options = optimoptions('linprog','Display','off');
	else
    		linprog_options = optimset('Display','off');
	end
        [x, ~, exitflag] = linprog(obj, D, d, N, b, lb, ub,linprog_options);
    case 2 % CPLEX      % exitflag ==  1: optimal, feasible  % exitflag == -2: infeasible
        [x, ~, exitflag] = cplexlp(obj, D, d, N, b, lb, ub);
    otherwise
        error('define a valid solver.');
end

if (exitflag == 1 && (solver == 2 || solver == 1)) ... CPLEX & Linprog
        || ((exitflag == 2 || exitflag == 5) && solver == 0) % GLPK
    feasible = 1;
    r = x;
else
    feasible = 0;
    return
end

if parsimonious == 1
fixed_fluxes = find(sum(logical(D),1));

% split network into forward and backward reactios
N2  = [N,-N];
ub2 = [ub;max([zeros(size(lb,1),1),-lb]')'];
lb2 = [max([zeros(size(lb,1),1),lb]')';-min([zeros(size(lb,1),1),ub]')'];
% fixate the already optimized fluxes
ub2(fixed_fluxes) = max([zeros(length(fixed_fluxes),1) r(fixed_fluxes)]')';
lb2(fixed_fluxes) = ub2(fixed_fluxes);
ub2(fixed_fluxes+cnap.numr) = -min([zeros(length(fixed_fluxes),1) r(fixed_fluxes)]')';
lb2(fixed_fluxes+cnap.numr) = ub2(fixed_fluxes+cnap.numr);
obj2 = ones(2*cnap.numr, 1);

switch solver
    case 0 % GLPK
        ctype = repmat('S',1,cnap.numis);
        [x2, ~, exitflag] =  glpk(  obj2,         N2, b, lb2, ub2, ctype);
    case 1 % linprog
        [x2, ~, exitflag] = linprog(obj2, [], [], N2, b, lb2, ub2, linprog_options);
    case 2 % CPLEX
        [x2, ~, exitflag] = cplexlp(obj2, [], [], N2, b, lb2, ub2);
    otherwise
        error('define a valid solver.');
end
if (exitflag == 1 && (solver == 2 || solver == 1)) ... CPLEX & Linprog
        || ((exitflag == 2 || exitflag == 5) && solver == 0) % GLPK
    r = x2(1:cnap.numr) - x2(cnap.numr+1:end);
else
    error('numerical error');
end
end
end
