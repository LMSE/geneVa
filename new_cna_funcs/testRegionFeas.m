function feas = testRegionFeas(cnap,V,v,solver)
% Tests if region(s) are feasible (overlap with the flux space) by performing an FBA 
% with the respective additional constraints: V*r <= v.
% Returns feasibility status for all sets of constraints.
% A single set of V and v can be provided as matrix and vector, multiple sets must
% be nested in a cell array that contains these matrices and vectors for each field.
if nargin < 4
    LPavail=LP_solver_availability(true);
    if LPavail(3)
        solver = 2;
    else
        solver = 0;
    end
end
if isnumeric(V) && isnumeric(v) 
	V = {V};
	v = {v};
end
feas = nan(1,length(v));
for i = 1:length(v)
	feas(i) = any(~isnan(CNAoptimizeFlux(cnap, [], [], solver, -1, 0, V{i}, v{i})));
end
end
