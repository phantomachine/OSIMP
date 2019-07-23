function bar2stack(X,Y)

Yneg = Y;
Yneg(Yneg>0) = 0;
Ypos = Y;
Ypos(Ypos<0) = 0;
hold on;
bar(X,Yneg,'stack');
bar(X,Ypos,'stack');
hold off;