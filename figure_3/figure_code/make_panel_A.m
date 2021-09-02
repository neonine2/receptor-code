mu = 0;
kappa = 80;
vonmises = @(x) exp(kappa*cos(x-mu))./(2*pi*besseli(0,kappa));
angle = linspace(-pi,pi,1000)';
unitcircle = [cos(angle),sin(angle)];
receptor = 2.3*vonmises(angle)+10;

flat = [1,100,200];
for ii = 1:3
    coord = circshift(movmean(receptor,flat(ii)),1000/4) .* unitcircle;
    plot(coord(:,1),coord(:,2))
    hold on
end
hold off
box off
pbaspect([2,3,1])