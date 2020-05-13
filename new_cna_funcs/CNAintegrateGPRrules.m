function [ gcnap, rmap, gmap, gkoCost, gkiCost, gT, gD, rType, sType, gpr_rules ] = ...
    CNAintegrateGPRrules( cnap, gpr_rules, koCost, kiCost, T, D, gkoCost, gkiCost, expand_notknockable )
%
% ------------------------------------------------
% CellNetAnalyzer API function 'CNAintegrateGPRrules'
% ------------------------------------------------
% --> Extends a given stoichiometric network with pseudoreactions and 
%     pseudometabolites for genes and enzyme pools.
%
% Usage: [ gcnap, rmap, gmap, gkoCost, gkiCost, gT, gD, rType, sType, gpr_rules ] = ...
%     CNAintegrateGPRrules( cnap, gpr_rules, koCost, kiCost, T, D, gkoCost, gkiCost, expand_notknockable )
%
% r1 (g1 or g2) translates to pseudoreactions:
%  g1:    -> E1 
%  p_g1,r1:  E1 -> Q1
%  g2:    -> E2
%  p_g2,r1:  E2 -> Q1
%  Q1 is consumed in r1
%
% r2 (g1 and g2) translates to pseudoreactions:
%  g1:       -> E1 
%  g2:       -> E2
%  p_g1,g2,r2:  E1 + E2 -> Q2
%  Q2 is consumed in r2
%
% Reversible reactions are split in forward and reverse direction. Reactions 
% that are catalyzed by multiple enzymes (containing 'or'-rules) are not copied.
%
% ------------------------------------------------
% Input:
%   cnap            : CellNetAnalyzer mass-flow project. (with gene rules in reaction notes)
%  ---- optional ----
%   gpr_rules       : A struct-array that contains the fiels reaction(s) [double] 
%                      and genes [double] that represent the gene-enzyme-reaction
%                      associations.
%   koCost, kiCost  : Vectors that specify the knockable reactions.
%                      Length must be identical to the number of reactions
%                      of the network. 'nan' -> reaction is not knockable.
%                      Knock-Ins 'override' knock-outs.
%                      If reactions with no gene association are marked
%                      as knockable (such as substrate or O2 uptake), 
%                      they will keep their 'knockability'.
%                      In case of enzyme associated reactions, all underlying
%                      genes (or enzymes) will inherit the KO/KI costs and
%                      the reaction itself will become not-knockable.
%                      In case of conflicts, when two reactions share the
%                      same gene, but one is knockable and the other is not,
%                      the gene will be set 'not-knockable' / 'nan'.
%   T, D            : Target and Desired vectors that should be reshaped.
%   gkoCost, gkiCost: Vectors with knockout and knock-In costs for each gene. 
%                      The length of these vectors must match the number of genes 
%                      occur in all rules. A way to get a complete list of all genes
%                      is via using [~,~,genes,~] = CNAgenerateGPRrules(cnap)
%                      or           [~,~,genes,~] = CNAgenerateGPRrules(cnap,grRules);
%                      Every entry represents the individual intervention costs.
%                      Notknockables carry 'nan'. Knock-Ins 'override' knock-outs.
%   expand_notknockable:Add enzyme and gene-association, even if their reactions
%                       were marked as notknockable 
%                      (koCost(i) = nan) (default: 1).
% ------------------------------------------------
% Output:   
%   gcnap            : CellNetAnalyzer mass-flow project, extended by
%                       1. source reactions for every gene product.
%                       2. conversion reactions from genes product to enzyme
%                          pool pseudometabolites defined by the reaction rules.
%                       3. pseudo-metabolite consumption is added to catalyzed reactions.
%   rmap             : Matrix that maps reactions from the extended to the original network
%                       This map is needed because some reactions are split and the
%                       original reaction indices are not maintained.
%   gmap             : Matrix that maps the list of genes on the gene source reactions in
%                       the new stoichiometric matrix.
%   gkoCost, gkiCost : Knockout and Knock-In costs for all reactions and
%                       gene-pseudoreactions in the extended network. These
%                       vectors are derived from the input vectors of knockout-costs.
%   gT,gD            : Target and Desired constraints mapped on the extended network.
%   rType            : indicates whether reactions in extended network are genes ('g'),
%                       GPR pseudoreactions ('p') or metabolic reactions ('r').
%   sType            : indicates whether species in extended network are genes ('g'),
%                       reaction-pseudometabolites ('p') or actual metabolites ('m').
%   gpr_rules        : Either same as entered or derived from reaction notes of CNA project.
%

%
% This file is part of CellNetAnalyzer. Please visit
% http://www.mpi-magdeburg.mpg.de/projects/cna/cna.html
% for more information and the latest version of CellNetAnalyzer.
%
% Copyright (C) 2000-2020 by Steffen Klamt and Axel von Kamp,
% Max Planck Institute for Dynamics of Complex Technical Systems, Magdeburg, Germany.
%
% Contributors are listed in CONTRIBUTORS.txt.
%
% This software can be used under the terms of our CellNetAnalyzer License.
% A copy of the license agreement is provided in the file named "LICENSE.txt"
% included with this software distribution. The license is also available online at
% http://www2.mpi-magdeburg.mpg.de/projects/cna/license.html
%
% For questions please contact: cellnetanalyzer@mpi-magdeburg.mpg.de
%

if nargin < 2 || isempty(gpr_rules)
    % if enzymes don't have names, generate names
    % if genes don't have names, generate names
    % check if cnap has enzyme
    if ~isfield(cnap,'enzymes')
        [cnap,~,~,gpr_rules] = CNAgenerateGPRrules(cnap);
    else
        gpr_rules = cnap.gpr_rules;
    end
end
numRules = length(gpr_rules);
if nargin < 3 || isempty(koCost)
    koCost = ones(cnap.numr,1);
end
if nargin < 4 || isempty(kiCost)
    kiCost = nan(cnap.numr,1);
end
if nargin < 5
    D = {};
    T = {};
end
if nargin < 7 || isempty(gkoCost) % number of genes is defined by 
    numGenes = max([gpr_rules(:).genes]);
    gkoCost = ones(numGenes,1);
    gkoCost_provided = 0;
else
    numGenes = length(gkoCost);
    gkoCost_provided = 1;  
end
if nargin < 8 || isempty(gkiCost)
    gkiCost = nan(numGenes,1);
    gkiCost_provided = 0;   
else, gkiCost_provided = 1;
end
if nargin < 9
    expand_notknockable = 1;
end

%% 0 initialize
koCost = double(koCost(:)'); % reshape as row (double)
kiCost = double(kiCost(:)');
koCost(~isnan(kiCost)) = nan; % knock-ins 'override' knock-outs

gkoCost = double(gkoCost(:)'); % reshape as row (double)
gkiCost = double(gkiCost(:)');
gkoCost(~isnan(gkiCost)) = nan; % knock-ins 'override' knock-outs

ruleMat = repmat({[]},1,cnap.numr);
for i = 1:cnap.numr
    rule_idx = find(ismember([gpr_rules(:).reaction],i));
    rule_i = arrayfun(@(x) gpr_rules(x).genes,rule_idx,'UniformOutput',0);
    ruleMat{i} = sparse(repelem_loc(1:length(rule_idx),cellfun(@length,rule_i)),cell2mat(rule_i),1,length(rule_idx),numGenes);
end

if ~expand_notknockable
    rule_removelist = zeros(1,length(gpr_rules));
    [gene_abund,genes] = hist([gpr_rules(:).genes],unique([gpr_rules(:).genes]));
    if ~gkoCost_provided
        for i = find(isnan(koCost) & ~cellfun(@isempty,ruleMat)) % for notknock reaction with gene reaction rule:
            % if a reaction is notknockable and there are genes that only occur with this reaction,
            % set those genes to notknockable, because it is never advantageous to knock them out.
            genes_1 = [gpr_rules([gpr_rules(:).reaction] == i).genes];
            for j = unique(genes_1)
                % if the number of occurrencies of this gene in the rules of one reaction is equal to
                % the total number of occurrencies of this gene in all rules.
                if sum(genes_1 == j) == gene_abund(genes == j)
                    gkoCost(j) = nan;
                end
            end
        end
    end
    if ~isempty(gkoCost)
        for i = 1:length(gpr_rules)
            % if all genes of one rule are notknockable, delete all rules for the same reaction
            % because the reaction can never be knocked out
            if all(isnan(gkoCost(gpr_rules(i).genes))) && all(isnan(gkiCost(gpr_rules(i).genes)))
                rule_removelist([gpr_rules(:).reaction] == gpr_rules(i).reaction) = 1;
            end
        end
    end
    gpr_rules = gpr_rules(~rule_removelist);
    numRules = length(gpr_rules);
end

if isempty(gpr_rules)
    error('Model extension with enzymes and genes not possible. Make sure that at least one enzyme can be knocked out.');
end

% set gene and enzyme names
if ~isfield(gpr_rules,'strGene')
    geneNames  = strcat('G',num2str(1:numGenes));
else
    [a,b] = unique([gpr_rules(:).genes]);
    gname = [gpr_rules(:).strGene]';
    geneNames = arrayfun(@(x) ['gene ' num2str(x)],1:numGenes,'UniformOutput',0);
    geneNames(a) = gname(b);
end
if ~isfield(gpr_rules,'name')
    ruleNames = strcat('PM',num2str(1:length(gpr_rules)));
else
    ruleNames = [gpr_rules(:).name]';
    if any(ismember(ruleNames,geneNames))
        ruleNames = strcat('Rule-',ruleNames);
    end
end

% split/duplicate reactions that are reversible and have a gene association
rules_per_reac = hist([gpr_rules(:).reaction],1:cnap.numr);
rev_reac_w_rule = rules_per_reac & (sign(cnap.reacMin) == -1)'; % subset of these reactions that are active in reverse direction
for_reac_w_rule = rules_per_reac & (sign(cnap.reacMax) ==  1)'; % subset of these reactions that are active in forward direction
% which reactions need to be split - create mapping
r_old2new = num2cell(1:cnap.numr);
for i = find(rev_reac_w_rule | for_reac_w_rule) % indicate which reactions are cloned (+/-)
    r_old2new{i} = [i*ones(1,for_reac_w_rule(i)) -i*ones(1,rev_reac_w_rule(i))]; % reactions
end
r_new2old = [r_old2new{:}];
rules_per_reac_new = rules_per_reac(abs(r_new2old));
rType(1:length(r_new2old)) = 'r'; % reaction Type: enzyme catalyzed reactions
sType(1:cnap.nums) = 'm'; % species Type: metabolites
rmap = full(sparse(abs(r_new2old),1:length(r_new2old),sign(r_new2old)));
 % reactions that stem back from reversible reactions with enzymes (and were split or inverted)
 % 0: reaction was not split, 1: positive sense reaction part, -1: reverse sense reaction part
r_was_rev = rev_reac_w_rule*rmap;

%% 1. extend stoichMat
%
% Take original model + gen-enzyme-reaction-association
% Split/duplicate reversible enzyme catalyzed reactions if necessary:
%   extend stoichMat with k reaction columns to split up reversible
%                          reactions

% Model: Stoichmat + bounds + objFunc
gcnap.stoichMat = cnap.stoichMat* rmap; % new Stoichmat
gcnap.reacMin = zeros(size(gcnap.stoichMat,2),1);
gcnap.reacMax = zeros(size(gcnap.stoichMat,2),1);
gcnap.reacMin(~r_was_rev) = cnap.reacMin(abs(r_new2old(~r_was_rev)));
gcnap.reacMin(~~r_was_rev) = 0;
gcnap.reacMax(r_was_rev >= 0) =  cnap.reacMax(abs(r_new2old(r_was_rev >= 0)));
gcnap.reacMax(r_was_rev <  0) = -cnap.reacMin(abs(r_new2old(r_was_rev <  0)));
gcnap.objFunc   = rmap'*cnap.objFunc;
% reacIDs
% Counter for split reactions (only for naming)
gcnap.reacNotes = cnap.reacNotes(abs(r_new2old));
gcnap.reacID    = cellstr(cnap.reacID(abs(r_new2old),:));
gcnap.reacID(rules_per_reac_new>0,:) = strcat(gcnap.reacID(rules_per_reac_new>0,:),'_gen');
gcnap.reacID(r_was_rev==-1,:) = strcat(gcnap.reacID(r_was_rev==-1,:),'_rev');
gcnap.reacID    = char(gcnap.reacID);
% specIDs
gcnap.specID    = cnap.specID;
gcnap.specNotes = cnap.specNotes;

%% 2. add reaction-pseudometabolites and genes
%   extend stoichMat with m1 rows for reaction-pseudometabolites
%                                          (all have values -1 below N
%                                           and 1 to the bottom right of N)
%   extend stoichMat with m2 rows for genes (all bottom right corner)
%   extend stoichMat with n1 = m1 columns for gene-to-reaction-pseudometabolite reactions
%   extend stoichMat with n2 = m2 columns for gene-producing reactions

% create vector that indicates whether reaction is gene expression, gene rule ('p') or reaction
reac_pseudo_met = unique(find(rules_per_reac)); % 1 pseudometabolite for each reaction
numReacPseudom = length(reac_pseudo_met);       %   (not 2 pseudo m. if reaction was split due to rev.)

% add pseudometabolites and pseudometabolites consumption to reactions
sType(end+1:end+numReacPseudom) = 'p';
pseudomet_consump_x = find(ismember(abs(r_new2old),reac_pseudo_met));
pseudomet_consump_y  = repelem_loc(1:numReacPseudom,(rev_reac_w_rule(reac_pseudo_met) & for_reac_w_rule(reac_pseudo_met))+1);
gcnap.stoichMat(sType=='p',rType=='r') = -sparse(pseudomet_consump_y,pseudomet_consump_x,1,numReacPseudom,size(gcnap.stoichMat,2));

% add pseudometabolite source terms (gene consumption is added later)
rType(end+1:end+numRules) = 'p';
pseudomet_source_x = 1:length(gpr_rules);
pseudomet_source_y = arrayfun(@(x) find(reac_pseudo_met == x), [gpr_rules(:).reaction]);
gcnap.stoichMat(sType=='p',rType=='p') = sparse(pseudomet_source_y,pseudomet_source_x,1,numReacPseudom,numRules); % add pseudometabolite source reaction
gcnap.reacMin(rType=='p') = 0;
gcnap.reacMax(rType=='p') = 1000;
gcnap.specID = char([cellstr(gcnap.specID) ; ...
                     strcat('psmet-',cellstr(cnap.reacID(reac_pseudo_met,:)))]);
gcnap.reacID = char([cellstr(gcnap.reacID) ; strcat('GPR-',ruleNames)]);

% add genes, gene consumption / gene rule part to pseudometabolite, also add gene source terms
sType(end+1:end+numGenes) = 'g';
rType(end+1:end+numGenes)   = 'g';
% change pseudoenzyme source terms into gene-pseudoenzyme reactions
pseudomet_gene_x = repelem_loc(1:length(gpr_rules),cellfun(@length,{gpr_rules.genes}));
pseudomet_gene_y = [gpr_rules.genes];
gcnap.stoichMat(sType=='g',rType=='p') = -sparse(pseudomet_gene_y,pseudomet_gene_x,1,numGenes,numRules);
gcnap.stoichMat(sType=='g',rType=='g') = -eye(numGenes); % add gene source pseudo-reaction
gcnap.reacMin(rType=='g') = -1000;
gcnap.reacMax(rType=='g') = 0;
unused_genes = setdiff(1:numGenes,[gpr_rules(:).genes]); % deactivate unused genes
gcnap.reacMin(find(rType=='g',1)-1+unused_genes) = 0;
gkoCost(unused_genes) = nan;
gkiCost(unused_genes) = nan;
gcnap.specID = char([cellstr(gcnap.specID) ; geneNames(:)]);
gcnap.reacID = char([cellstr(gcnap.reacID) ; strcat('GP-',geneNames(:))]);

gcnap.specNotes(size(gcnap.stoichMat,1)) = {''};
gcnap.reacNotes(size(gcnap.stoichMat,2)) = {''};
gcnap.objFunc(rType=='p' | rType=='g')   = 0;
rmap = [rmap, zeros(cnap.numr,numGenes+numRules)];

%% 3. adapt kiCost, koCost, D, d, T, t
% kiCost and koCost vectors will be derived from reaction-Cost vectors
% if reactions have associated enzymes/genes, those will be marked as knockable.
% Else, the reactions will stay knockable. Conflicts are treated as follows:
% When geneA catalyzes R1 (knockable) and R2 (notknockable), the corresponding geneA 
% will be knockable. If R1 has a lower cost than R2, lower cost will be taken into account.

gmap = [zeros(numGenes,sum(rType ~= 'g')), eye(numGenes)];
if gkoCost_provided
    gkoCost = [nan(1,sum(rType ~= 'g')) , gkoCost];
else
    gkoCost = nan(1,length(rType));
end
if gkiCost_provided
    gkiCost = [nan(1,sum(rType ~= 'g')) , gkiCost];
else
    gkiCost = nan(1,length(rType));
end
cgenes = {gpr_rules(:).genes};

for i = 1:size(gcnap.stoichMat,2)
    switch rType(i)
        case 'r' % reactions stay knockable if they were knockable before and don't have a gene-association
            if ~rules_per_reac_new(i)
                gkoCost(i) = koCost(r_new2old(i));
                gkiCost(i) = kiCost(r_new2old(i));
            end
        case 'p'
        case 'g'
            if ~gkoCost_provided || ~gkiCost_provided
                gene_id = find(find(rType=='g') == i);
                rules = cellfun(@(x) ismember(gene_id,x),cgenes);
                if any(~isnan(koCost([gpr_rules(rules).reaction]))) && ~gkoCost_provided
                    gkoCost(i) = min(koCost([gpr_rules(rules).reaction]));
                end
                if any(~isnan(kiCost([gpr_rules(rules).reaction]))) && ~gkiCost_provided
                    gkiCost(i) = min(kiCost([gpr_rules(rules).reaction]));
                end
            end
        otherwise
            error(['invalid rType. Should be ''r'' for reactions, ''e'''...
                   'for enzyme rules/sources and ''g'' for gene sources']);
    end
end
if any(~isnan(gkiCost) & ~isnan(gkoCost))
    warning('There are conflicts between KIable and KOable vector. Please check.');
end

if ~gkoCost_provided
    knockable_genes = gkoCost(rType=='g');
    not_knockable_warning = '';
    for i = find(isnan(koCost) & rules_per_reac > 0) % check if notknockable reactions still depend on at least one rule that is notknockable
        rules = find([gpr_rules(:).reaction] == i);
        independent = 0;
        for j = rules
            if all(isnan(knockable_genes(gpr_rules(j).genes)))
                independent = 1;
            end
        end
        if ~independent
                not_knockable_warning = [not_knockable_warning strtrim(cnap.reacID(i,:)) ...
            ' :RULE: ' num2str(rules) ...
            ' :GENE: ' num2str([gpr_rules(rules).genes])  '; ' newline];
        end
    end
    if ~isempty(not_knockable_warning)
        warning(['A reaction that was before not knockable has become knockable through the expansion. '...
            'The affected reactions, enzymes and genes are ' newline not_knockable_warning 'Please check for ' ...
            'genetic overlaps of knockable and notknockable reactions before expanding the model.']);
    end
end
% T and D
if ~isempty(T)
    for i = 1:length(T) % fill up to match reaction vector size
        gT{i}  = T{i}*rmap;
    end
else
    gT = {};
end
if ~isempty(D)
    for i = 1:length(D) % fill up to match reaction vector size
        gD{i}  = D{i}*rmap;
    end
else
    gD = {};
end

gcnap = CNAgenerateMFNetwork(gcnap);
%% 4. in case of gui, generate map
if cnap.has_gui
    gcnap.local.errval  = 0;
    gcnap.numr          = size(gcnap.stoichMat,2);
    gcnap.nums          = size(gcnap.stoichMat,1);
    gcnap.mue           = find(r_new2old == cnap.mue');
    gcnap.figs          = cnap.figs;
    gcnap.path          = cnap.path;
    gcnap.color1        = cnap.color1;
    gcnap.color2        = cnap.color2;
    gcnap.color3        = cnap.color3;
    gcnap.color4        = cnap.color4;
    gcnap.maps          = cnap.maps;
    gcnap.nummaps       = cnap.nummaps;
    gcnap.show_flux_format= cnap.show_flux_format;
    gcnap.net_var_name  = cnap.net_var_name;
    gcnap.reacDefault   = cnap.reacDefault(abs(r_new2old)).*sign(r_new2old);
    gcnap.reacDefault(r_was_rev == -1) = double(gcnap.reacDefault(r_was_rev == -1)>=0).*gcnap.reacDefault(r_was_rev == -1);
    gcnap.reacFontSize  = cnap.reacFontSize;
    gcnap.reacBoxWidth  = cnap.reacBoxWidth;
    gcnap.reacBoxHeight = cnap.reacBoxHeight;
    gcnap.specFontSize  = cnap.specFontSize;
    gcnap.specBoxWidth  = cnap.specBoxWidth;
    gcnap.specBoxHeight = cnap.specBoxHeight;
    gcnap.specBoxColor  = cnap.specBoxColor;
    %% extend reac Boxes vector
    gcnap.reacBoxes = cnap.reacBoxes(abs(r_new2old),:);
    gcnap.reacBoxes(r_was_rev == -1, 4) = 0; % reversible part
    gcnap.reacBoxes(r_was_rev == -1, 2) = gcnap.reacBoxes(r_was_rev == -1, 2) + 20;
    gcnap.reacBoxes(r_was_rev == -1, 3) = gcnap.reacBoxes(r_was_rev == -1, 3) + 3;
    disp('drawing new reaction boxes');
    for reacIndex = find(r_was_rev == -1)
        % declaration of handle
        cnan.open_projects.(gcnap.net_var_name).gui.handles= struct;
        gcnap = update_after_change(gcnap);
        currmap = gcnap.reacBoxes(reacIndex,5);
        % is visible/editable
        zw=gcnap.reacBoxes(reacIndex,6);
        fig = gcnap.figs(currmap,:);
        % make box
        zw1=uicontrol('Style', 'edit','Parent',fig(1), 'String', '###',...
            'Units','normalized','HorizontalAlignment','left','BackgroundColor',gcnap.color1,'ForegroundColor',gcnap.textColor,'TooltipString',gcnap.reacID(reacIndex,:));
        set(zw1, 'ButtonDownFcn', {@execute_callback, gcnap.net_var_name,...
            {'check_for_right_click', 'reaceditmask'}, {'reacenr', reacIndex}});
        % save handle
        gcnap.reacBoxes(reacIndex,4)=zw1;
        % adjust "zoom"
        place_box(fig,zw1,...
            gcnap.reacBoxes(reacIndex,2),...
            gcnap.reacBoxes(reacIndex,3),...
            gcnap.reacFontSize(currmap),...
            gcnap.reacBoxWidth(currmap),...
            gcnap.reacBoxHeight(currmap));
        % put default rate in textbox
        if isnan(gcnap.reacDefault(reacIndex))
            set(zw1,'String','#');
        else
            set(zw1,'String',num2str(gcnap.reacDefault(reacIndex)));
        end
        % is visible/editable
        if(zw==2)
            set(zw1,'Style', 'text');
        elseif(zw==3)  %%non-visible
            set(zw1,'Visible','off');
        end
    end
    disp('generating map for gene and pseudometabolite reactions');
    gcnap.reacBoxes(rType ~= 'r',5) = -1;
    gcnap = myCNAgenerateMap(gcnap,sum(rType ~= 'r')<=110);
end

% create gcnap Network
gcnap.epsilon = cnap.epsilon;
gcnap.rType = rType;
gcnap.sType = sType;
end

%% supplementary
function place_box(fig,handle,xp,yp,fontsize,box_width,box_height) % copied from zoom_single_box
    zzz=get(fig(1),'CurrentAxes');
    if(~fig(4))
    xxx=get(zzz,'XLIM');
    yyy=get(zzz,'YLIM');
    xxxl=xxx(2)-xxx(1);
    yyyl=yyy(2)-yyy(1);
    xp=(xp-xxx(1))/xxxl;
    yp=(yyy(2)-yp)/yyyl;
    box_width=box_width*fig(3)/xxxl;
    box_height=box_height*fig(2)/yyyl;
    fontsize=fontsize*fig(3)/xxxl;
    end
    pos=[xp yp box_width box_height];
    set(handle,'Position',pos,'FontSize',fontsize);
end

% local implementation of the MATLAB repelem function.
function U = repelem_loc(V,varargin)
% This function is used to guarantee compatibility with MATLAB2014b and below
    if exist('repelem','builtin')
        U = repelem(V,varargin{:});
        return;
    end
    % if V = 1xn or nx1 and N = 1xn
    if size(V,1) == 1 && length(varargin) == 1
        varargin{2} = 1;
        varargin = flip(varargin);
    elseif size(V,2) == 1 && length(varargin) == 1
        varargin{2} = 1;
    end
    % if V = 1xn and N = 1x1
    for i = find(cellfun(@length,varargin) == 1)
        varargin{i} = varargin{i}*ones(1,getelements(size(V)',i));
    end
	[reps{1:numel(varargin)}] = ndgrid(varargin{:});
	U = cell(cellfun(@length,varargin));
	U(:) = arrayfun(@(x) V(x)*ones(cellfun(@(y) y(x),reps)) ,1:numel(V),'UniformOutput',false);
	U = cell2mat(U);
end
