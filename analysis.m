global k
ret=[];
delta = 0.001;
tolerance = 0.001;
krange = 10 : -0.01 : 1;
Ra = 7000;


for k = krange
    r = 1;
    while abs(r) > tolerance
        f = detA(Ra);
        d = (detA(Ra + delta) - detA(Ra - delta) ) / (2 * delta);
        r = - f / d;
        Ra = Ra + r;
    end
    ret = [ret Ra];
end
plot ( krange, ret )
xlabel('Wavenumber k') , ylabel('Rayleigh Number')


function ret = detA(Ra)
    global k;
    a0= 1i * k * (-1 + (Ra / k^4) ^ (1 / 3)) ^ (1 / 2);
    a1= k * (1 + (Ra / k^4) ^ (1/3) * (1 / 2 + 1i * sqrt(3) / 2))^(1 / 2);
    a2= k * (1 + (Ra / k^4) ^ (1/3) * (1 / 2 - 1i * sqrt(3) / 2))^(1 / 2);
    
    A = [1, 1, 1, 1, 1, 1
        exp(a0), exp(-a0), exp(a1), exp(-a1), exp(a2), exp(-a2)
        a0, -a0, a1, -a1, a2, -a2
        a0*exp(a0), -a0*exp(-a0), a1*exp(a1), -a1*exp(-a1), a2*exp(a2), -a2*exp(-a2),
        (a0^2-k^2)^2, (a0^2-k^2)^2, (a1^2-k^2)^2, (a1^2-k^2)^2, (a2^2-k^2)^2, (a2^2-k^2)^2
        (a0^2-k^2)^2*exp(a0), (a0^2-k^2)^2*exp(-a0), (a1^2-k^2)^2*exp(a1), (a1^2-k^2)^2*exp(-a1), (a2^2-k^2)^2*exp(a2), (a2^2-k^2)^2*exp(-a2) ];
    
    ret=det(A);
end

