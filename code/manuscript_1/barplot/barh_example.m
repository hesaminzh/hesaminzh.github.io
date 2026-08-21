x = [1, 2, 3];
vals = [[2,3,6]; [11,23,26]];
b = barh(x,vals);

% Create Bar Labels Using text function
% display the values at the tips of the first series of bars
xtips1 = b(1).YEndPoints + 0.3;
ytips1 = b(1).XEndPoints;
labels = string(b(1).YData);
text(xtips1, ytips1, labels, 'VerticalAlignment','middle');


% display the values at the tips of the second series of bars
xtips2 = b(2).YEndPoints + 0.3;
ytips2 = b(2).XEndPoints;
labels2 = string(b(2).YData);
text(xtips2, ytips2, labels2,'VerticalAlignment','middle');
