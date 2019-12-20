data = All_data(1).CINE(1).mark(5).data;
x = data(1:3:30,:);
y = data(2:3:30,:);
z = data(3:3:30,:);

plot3(x,1:111,y);
grid on;
axis equal;

t = linspace(0,2*pi);
for i = 1 : 111
    dx(i)= std(x(:,i));
    dy(i)= std(y(:,i));
    c(i,:) = [mean(x(:,i)), mean(y(:,i))];
    
    X = dx(i)*cos(t) + c(i,1);
    Y = dy(i)*sin(t) + c(i,2);
    hold on
    plot3(X,repmat(i,length(X),1),Y)
end



