function plotmatrix(A)
% axis ij -> aqui inverte a mostra
while any(A(:) > 0)
  [r,c] = find(A>0);
  scatter(r,c,'filled','black')
  A = A-1;
end
grid on
end