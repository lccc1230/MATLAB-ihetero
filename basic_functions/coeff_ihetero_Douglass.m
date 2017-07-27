function c = coeff_ihetero_Douglass(alpha, K1, K2, l0, p10, p20)


% NOTE: alpha represents cooperativity as defined in the paper

a0 = (alpha-1)*alpha;

a1 = (alpha^2*(2*p10+l0+2*p20)+2*alpha*(K1+K2-p10-p20)-2*(K1+K2));

a2 = alpha^2*(p10^2+2*p10*(l0+2*p20)+p20*(2*l0+p20))+...
     alpha*(-p10*(l0+3*p20-2*K1-3*K2)+K2*(l0+2*p20+K1)-p10^2+l0^2-l0*p20-p20^2+l0*K1+3*p20*K1)+...
     (K2^2-2*K2*(p10+p20+K1)+K1*(K1-2*(p10+p20)))-...
     (K1-K2)^2/alpha;
 
a3 = alpha^2*(l0*p20^2+2*p10*p20*(2*l0+p20)+p10^2*(l0+2*p20))-...
    (p10*(l0+p20)+p20*(p20-l0-K1))*K1-...
    (p10^2+p20*(l0+K1)+p10*(p20-l0+K1))*K2+...
    alpha*(-p10^2*(l0+p20-K2)+...
            p10*(l0^2+l0*(K1+K2-2*p20)-p20^2+3*p20*(K1+K2)+K1*K2)+...
            p20*(l0^2+l0*(K1+K2-p20)+K1*(p20+K2)))+...
    p10*K2^2;

a4 = alpha*p10*p20*(l0^2+...
                      l0*(2*alpha*p20+K1+K2-p20)+...
                      p10*((2*alpha-1)*l0+alpha*p20+K2)+...
                      K1*(p20+K2));

a5 = alpha^2*p10^2*l0*p20^2;

c = [a0 -a1 a2 -a3 a4 -a5];