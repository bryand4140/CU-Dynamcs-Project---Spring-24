clc; clear;

%ASEN 5022: Dynamics of Aerospace Structures

%Blade Params
m_b = 6; %Mass of Blade (kg)
len = 1.5; %Length of Blade (m)
k_b = 300*(10^6); %Spring Constant of Blade (Nm)

n = 16; % Number of Blades
i = 1:n; % Individual Blade Designation

%External Forces Params
omega = 262.85; %Rotational Speed (rad/s)
L_b = 720; %Lift Force (N)

% Blade Defect Params
d_b = 12; % Defect Blade Number
blade_work_percent = 0.5;

M_blade = (m_b)*eye(2*n);
M_blade(d_b,d_b) = M_blade(d_b,d_b)*blade_work_percent;
M_blade(n+d_b,n+d_b) = M_blade(n+d_b,n+d_b)*blade_work_percent;

M = trace(M_blade(1:n,1:n))*eye(2);
K = k_b*eye(2);
T = [eye(2) [cos((2*pi*(i-1))/n) -sin((2*pi*(i-1))/n);sin((2*pi*(i-1))/n) cos((2*pi*(i-1))/n)];[cos((2*pi*(i-1))/n)' sin((2*pi*(i-1))/n)';-sin((2*pi*(i-1))/n)' cos((2*pi*(i-1))/n)'] eye(2*n)];

F_blade = [M_blade(1:n,1:n)*ones(n,1)*(omega^2)*len;L_b*ones(n,1)];
F_blade(n+d_b,1) = F_blade(n+d_b,1)*(blade_work_percent^4);

F = T*[zeros(2,1);F_blade];
S = F(1:2,1);

%%Plot Params
t_end = 10;
del_t = 0.1;
t = linspace(0,t_end,t_end/del_t);
plot(t,S(1)*sin(omega*t),'r',t,S(2)*sin(omega*t),'k',LineWidth=2);
legend('F_{x_r}','F_{y_r}');
xlim([0 t_end]);
xlabel('Time (s)');
ylabel('Force (N)');
title(['Force on Rotor with del t = ' num2str(del_t)]);
hold off;

%%Calc Params
[V,E] = eig(K\M);
X_r_calc = V(:,1);
X_s_calc = V(:,2);
lambda_r_calc = E(1,1);
omega_r_calc = sqrt(lambda_r_calc);
lambda_s_calc = E(2,2);
omega_s_calc = sqrt(lambda_s_calc);
mu_r_calc = X_r_calc'*M*X_r_calc;
mu_s_calc = X_s_calc'*M*X_s_calc;

%%Duhamels Integral

%Initial Conditions
duhamels_r.time = 0;
duhamels_r.nu = 0;
duhamels_r.nu_dot_del_t = 0;

duhamels_s.time = 0;
duhamels_s.nu = 0;
duhamels_s.nu_dot_del_t = 0;

%Duhamels-r Integral
while true
    P_r = [S(1)*sin(omega*duhamels_r.time(end));S(2)*sin(omega*duhamels_r.time(end))];
    P_next_r = [S(1)*sin(omega*(duhamels_r.time(end)+del_t));S(2)*sin(omega*(duhamels_r.time(end)+del_t))];
    phi_r = (X_r_calc'*P_r)/mu_r_calc;
    phi_next_r = (X_r_calc'*P_next_r)/mu_r_calc;
    alpha_r = omega_r_calc*del_t;

    duhamels_r.nu(end+1) = (cos(alpha_r)*duhamels_r.nu(end)) + ((sin(alpha_r)/alpha_r)*duhamels_r.nu_dot_del_t(end)) + ((del_t/alpha_r)^2)*(((sin(alpha_r)/alpha_r)-cos(alpha_r))*phi_r + (1-(sin(alpha_r)/alpha_r))*phi_next_r);
    duhamels_r.nu_dot_del_t(end+1) = (-alpha_r*sin(alpha_r)*duhamels_r.nu(end)) + (cos(alpha_r)*duhamels_r.nu_dot_del_t(end)) + ((del_t/alpha_r)^2)*((alpha_r*sin(alpha_r)+cos(alpha_r)-1)*phi_r + (1-cos(alpha_r))*phi_next_r);
    duhamels_r.time(end+1) = duhamels_r.time(end)+del_t;
    
    if duhamels_r.time(end) >= t_end
        break;
    end
end

%Duhamels-s Integral
while true
    P_s = [S(1)*sin(omega*duhamels_s.time(end));S(2)*sin(omega*duhamels_s.time(end))];
    P_next_s = [S(1)*sin(omega*(duhamels_s.time(end)+del_t));S(2)*sin(omega*(duhamels_s.time(end)+del_t))];
    phi_s = (X_s_calc'*P_s)/mu_s_calc;
    phi_next_s = (X_s_calc'*P_next_s)/mu_s_calc;
    alpha_s = omega_s_calc*del_t;

    duhamels_s.nu(end+1) = (cos(alpha_s)*duhamels_s.nu(end)) + ((sin(alpha_s)/alpha_s)*duhamels_s.nu_dot_del_t(end)) + ((del_t/alpha_s)^2)*(((sin(alpha_s)/alpha_s)-cos(alpha_s))*phi_s + (1-(sin(alpha_s)/alpha_s))*phi_next_s);
    duhamels_s.nu_dot_del_t(end+1) = (-alpha_s*sin(alpha_s)*duhamels_s.nu(end)) + (cos(alpha_s)*duhamels_s.nu_dot_del_t(end)) + ((del_t/alpha_s)^2)*((alpha_s*sin(alpha_s)+cos(alpha_s)-1)*phi_s + (1-cos(alpha_s))*phi_next_s);
    duhamels_s.time(end+1) = duhamels_s.time(end)+del_t;

    if duhamels_s.time(end) >= t_end
        break;
    end
end

q_duhamels = duhamels_r.nu'*X_r_calc' + duhamels_s.nu'*X_s_calc';

%%Plot
plot(duhamels_r.time,q_duhamels(:,1),'r',duhamels_s.time,q_duhamels(:,2),'k',LineWidth=2);
legend('x_r','y_r');
xlim([0 t_end]);
xlabel('Time (s)');
ylabel('Amplitude (m)');
title(['Response of Rotor with del t = ' num2str(del_t)]);
