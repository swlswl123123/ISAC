function CPM_LFM_Total = FreSum(CPM_LFM_Total, CPM_LFM)
%PlotFre - Description
%
% Syntax: output = PlotFre(CPM, LFM)
%
% Long description
    if isempty(CPM_LFM_Total)
        CPM_LFM_Total = CPM_LFM;
    else
        CPM_LFM_Total = CPM_LFM_Total + CPM_LFM;
    end
end