clear all
R = 8.314;
T = 350;
k1 = 2.20 * 10^6 * exp(-65.6*10^3 / (R*T) );
k2 = 2.32 * 10^4 * exp(-51.4*10^3 / (R*T) );
k3 = 1.14 * 10^5 * exp(-60.3*10^3 / (R*T) );
k4 = 3.58 * 10^3 * exp(-48.3*10^3 / (R*T) );
k5 = 3.36 * 10^2 * exp(-45.9*10^3 / (R*T) );
k6 = 4.80 * 10^1 * exp(-40.2*10^3 / (R*T) );
k7 = 6.49 * 10^0 * exp(-34.9*10^3 / (R*T) );
k8 = 4.37 * 10^0 * exp(-39.1*10^3 / (R*T) );
molar_density = 9.969 * 10^3;
xA(1) = 0.024028;
xB1(1) = 0.528713;
xcB2(1) = 0;
xtB2(1) = 0;
xBA(1) = 0.053561;
xIB(1) = 0.331857;
xIBA(1) = 0.056199;
xH2 = 0.04;
CA(1) = molar_density * xA(1);
CB1(1) = molar_density * xB1(1);
CcB2(1) = molar_density * xcB2(1);
CtB2(1) = molar_density * xtB2(1);
CBA(1) = molar_density * xBA(1);
CIB(1) = molar_density * xIB(1);
CIBA(1) = molar_density * xIBA(1);
xA_allowed(1) = 0.001;
xA_target(1) = 0.00025;
CH2 = molar_density * xH2;
residence_time =200;
step_count = 1000;
dt = residence_time / step_count;
time_count(1) = 0;
for i = 2:step_count;
    CA(i) =  CA(i-1) + dt*(-k1-k2-k3-k4)*CA(i-1)*(CH2)^0.5;
    xA(i) = CA(i) / molar_density;
    CB1(i) = CB1(i-1) + dt*(k1*CA(i-1)*(CH2)^0.5 - (k5+k6+k7)*CB1(i-1)*(CH2)^0.5);
    xB1(i) = CB1(i) / molar_density;
    CcB2(i) = CcB2(i-1) + dt*(k3*CA(i-1)*(CH2)^0.5 + k6*CB1(i-1)*(CH2)^0.5);
    xcB2(i) = CcB2(i) / molar_density;
    CtB2(i) = CtB2(i-1) + dt*(k2*CA(i-1)*(CH2)^0.5 + k5*CB1(i-1)*(CH2)^0.5);
    xtB2(i) = CtB2(i) / molar_density;
    CBA(i) = CBA(i-1) + dt*(k4*CA(i-1)*(CH2)^0.5 + k7*CB1(i-1)*(CH2)^0.5);
    xBA(i) = CBA(i) / molar_density;
    CIB(i) = CIB(i-1) - dt*(k8*CIB(i-1)*(CH2)^0.5);
    xIB(i) = CIB(i) / molar_density;
    CIBA(i) = CIBA(i-1) + dt*(k8*CIB(i-1)*(CH2)^0.5);
    xIBA(i) = CIBA(i) / molar_density;
    xA_allowed(i) = 0.001;
    xA_target(i) = 0.00025;
    time_count(i) = time_count(i-1) + dt;
end
plot(time_count, xA,'b','linewidth',1)
xlabel('residence time / s')
ylabel('mole fraction')
hold on
plot(time_count,xA_allowed, 'r','LineStyle','--')
plot(time_count,xA_target, 'g','LineStyle','--')
legend('1,2-butadiene', 'required concentration','target concentration')