function y = plotgrid(xcells,ycells)
ny = length(ycells);
nx = length(xcells);
greycolor = [211 211 211]/256;
lw = 0.1;
% plotar ny linhas horizontais
% plotar nx linhas verticais
y = figure;
hold on
xmin = min(xcells);
xmax = max(xcells);
ymin = min(ycells);
ymax = max(ycells);
for i=1:ny
   plot([xmin xmax],[ycells(i) ycells(i)],'color',greycolor,'linewidth',lw) 
end

for i=1:nx
   plot([xcells(i) xcells(i)],[ymin ymax],'color',greycolor,'linewidth',lw) 
end



end