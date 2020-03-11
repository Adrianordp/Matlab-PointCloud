function plotmatrix(A,level)
% axis ij -> aqui inverte a mostra
  [r,c] = find(A>level); %A > 0
  scatter(r,c,'filled','red')
  hold on
  [r,c] = find(A==-1); %A > 0
  scatter(r,c,'filled','green')
  
  grid on
  
end