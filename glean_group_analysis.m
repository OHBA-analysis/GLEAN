%function glean_group_analysis(GLEAN,design)

res = glean_group_temporal_stats(GLEAN,design.ev,design.contrasts);

for i = 1:size(design.ev,2)
    reg = unique(design.grouplabels(design.ev(:,i) == 1))
end



for resname = fieldnames(res)'
    glean_plot_hmm_stats(res.(char(resname)),design.grouplabels,design.conlabels);
end