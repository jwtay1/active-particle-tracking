function expZZ = fngauss(params, xx, yy)

A = params(1);

x0 = params(2);
y0 = params(3);

sigmaX = params(4);
sigmaY = params(5);

offset = params(6);

expZZ = A * exp(- ( ((xx - x0).^2)./(2 * sigmaX^2) +  ((yy - y0).^2)./(2 * sigmaY^2) )) + offset;

end