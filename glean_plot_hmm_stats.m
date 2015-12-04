function glean_plot_hmm_stats(res,group_labels,contrast_labels)
% Pretty plotting of HMM statistics:


num_states      = size(res.stats,2);
num_contrasts   = size(res.tstats,1);

if ~exist('group_labels','var')
    group_labels = repmat({' '},1,size(res.stats,1));
end


% Plot statistics for each group as violin plot
fig = figure('color','w');
for n = 1:num_states
    h = autosubplot(n,num_states);
    violinplot(res.stats(:,n),group_labels,20);
    set(h,'linewidth',2);
    set(h,'fontsize',16,'fontWeight','bold');
    ylabel(sprintf('%s (%s)',res.label,res.units),'fontsize',16,'fontWeight','bold');
    set(findobj(gca,'Type','text'),'fontsize',16,'fontWeight','bold')
    title(['state ' num2str(n)],'fontsize',16,'fontWeight','bold');
end


% Plot t-statistics for each contrast as bar plot and show significance
fig = figure('color','w');
cols = prettylines(num_contrasts);
h = axes('parent',fig);
b = bar(res.tstats');
colormap(cols);
set(h,'linewidth',2);
set(b,'linewidth',2);
set(h,'fontsize',16,'fontWeight','bold');
ylabel('t-statistic','fontsize',16,'fontWeight','bold');
xlabel('state','fontsize',16,'fontWeight','bold');
set(h,'xlim',[0 num_states+1])
legend(contrast_labels,'location','northoutside','orientation','horizontal')
legend boxoff

end

function h = autosubplot(n,N)
% [i j] = autosubplot(N,n)
% Creates a subplot at position n of N, automatically taking care of the
% plot arrangement
 
% Try to factor, then check aspect ratio is acceptable:
f(1) = max(factor(N));
f(2) = N / f(1);

if max(f)/min(f) <= 2
    i = max(f);
    j = min(f);
else
    i = ceil(sqrt(N));
    j = ceil(N/i);    
end

h = subplot(j,i,n);

end


function h = violinplot(data,group_labels,nbins)

cla, hold on

groups = unique(group_labels);

cols = prettylines(length(groups));

for i = 1:length(groups)
  
  Y = data(strcmp(group_labels,groups(i)));
  
  bins = linspace(nanmin(Y),nanmax(Y),nbins-2);
  bins = sort([bins bins(1)-mode(diff(bins)) bins(1)+mode(diff(bins))]);
  
  h = histc(Y,bins);
  h = 0.4*h./max(h);
  h = smooth(h,5);
  h = h(:)';
  patch([i-h i+h(end:-1:1)],[bins bins(end:-1:1)],cols(i,:),'EdgeColor','k','linewidth',2)
  
  plot(i,median(Y),'xk','MarkerSize',16,'LineWidth',2)  
  
end
legend(groups,'location','northoutside','orientation','horizontal')

end