clearvars
close all

% Weight matrix for the intermarriages.
W = [0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; % 1. Castellani
	1 0 1 0 0 0 1 0 0 0 0 0 0 0 0; % 2. Peruzzi
	1 1 0 0 1 0 1 0 0 0 0 0 0 0 0; % 3. Strozzi
	1 0 0 0 0 0 0 0 1 0 0 0 0 0 0; % 4. Barbadori
	0 0 1 0 0 0 0 0 1 0 0 0 0 0 0; % 5. Ridolfi
	0 0 0 0 0 0 0 0 1 0 0 0 0 0 0; % 6. Acciaiuoli
	0 1 1 0 0 0 0 0 0 1 0 0 0 0 0; % 7. Bischeri
	0 0 0 0 0 0 0 0 1 1 0 0 0 0 0; % 8. Tornabuoni
	0 0 0 1 1 1 0 1 0 0 1 1 0 0 0; % 9. Medici
	0 0 0 0 0 0 1 1 0 0 1 0 1 0 0; % 10. Guadagni
	0 0 0 0 0 0 0 0 1 1 0 0 0 1 0; % 11. Albizzi
	0 0 0 0 0 0 0 0 1 0 0 0 0 0 1; % 12. Salviati
	0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; % 13. Lamberteschi
	0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; % 14. Ginori
	0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]; % 15. Pazzi

% Calculating the invariant probability distribution P = D\W = diag(w)\W.
% w is the out-degree of each node, so the sum of each row
P = diag(sum(W,2))\W;

stubborn = [3 9];
nonstubborn = setdiff(1:15, stubborn);

% Creating matrices Q, R and S. P = [Q R; zeros S]
Q = P(nonstubborn,nonstubborn);
R = P(nonstubborn,stubborn);
S = P(stubborn,stubborn);

steps = 100;

x = zeros(15,steps);
x(stubborn,1) = [-1,1];

for k = 2:steps
    x(nonstubborn,k) = Q*x(nonstubborn,k-1) + R*x(stubborn,k-1);
    x(stubborn,k) = x(stubborn,k-1);
end

plot(x')
names = ["Castellani","Peruzzi","Strozzi","Barbadori","Ridolfi",...
    "Acciaiuoli","Bischeri","Tornabuoni","Medici","Guadagni",...
	"Albizzi","Salviati","Lamberteschi","Ginori","Pazzi"];
legend(names)

for k = 1:15
    fprintf('%12s \t %1.4f\n',names(k),x(k,end))
end





