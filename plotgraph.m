subplot(3,1,1)
x1=1:istep;
y1=E(:,1);
plot(x1,y1);
Legend('Kinetic');

subplot(3,1,2)
x2=1:istep;
E(:,2)=abs(E(:,2));
y2=E(:,2);
plot(x2,y2,'color','r');
Legend('Potential');

subplot(3,1,3)
x3=1:istep;
E(:,3)=E(:,1)+E(:,2);
y3=E(:,3);
plot(x3,y3,'color',[0 0.6 0]);
Legend('Total');