%% Scaling law for Rayliegh taylor

clear

P(1) = 1;
P(2) = 1;
P(3) = 1;
P(4) = 1;

grad(1)=  0.0000217880391690281237306210659898653148;
grad(2)= -0.0000215651432274706561856292275081159460;
grad(3)= -0.0000194068852722766583261161665863170356;
grad(4)=  0.0000000900509960614182950205494518174099;

Q = 0.000021565336561051446213855762;

exp = (grad .* P) ./ (Q);
% STD Scaling
tem = 1;
for i=1:length(P)
    tem = tem* P(i).^exp(i);
end
prefactor = (Q)./tem;

display(sprintf('Scaling law = %4.4e * %4.4e^%4.4e',prefactor,P(1),exp(1)))
for i=2:length(P)
    display(sprintf('                         * %4.4e^%4.4e',P(i),exp(i)))
end

abc = 1;







