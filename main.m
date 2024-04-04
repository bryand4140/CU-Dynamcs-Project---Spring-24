clc; clear;

%==========================================================================
%                             ** CONTROLS **
show_plots = false;

%==========================================================================

%Blade Params
m_b = 6; %Mass of Blade (kg)
len = 1.5; %Length of Blade (m)
k_b = 298*(10^6); %Spring Constant of Blade (Nm)

%Rotor Params
m_r = 3500; %Mass of Rotor (kg)
k_r = 300*(10^6); %Spring Constant of Rotor (Nm)

n = 16; %Number of Blades
i = 1:n; %Individual Blade Designation

%External Forces Params
omega = 262.85; %Rotational Speed (rad/s)
g = 9.81; %Acceleration due to gravity (m/s^2)

%Aero Force Params
L_b = 700; %Force (N)

% Blade Defect Params
d_b = 1; % Defect Blade Number
blade_work_percent = 0.5;

%--------------------------------------------------------------------------
%Calculations:
M_blade = (m_b)*eye(2*n);
M_blade(d_b,d_b) = M_blade(d_b,d_b)*blade_work_percent;
M_blade(n+d_b,n+d_b) = M_blade(n+d_b,n+d_b)*blade_work_percent;

K_blade = [k_b*eye(n) zeros(n,n); zeros(n,n) k_b*eye(n)];
K_blade(d_b,d_b) = K_blade(d_b,d_b)/blade_work_percent;
K_blade(n+d_b,n+d_b) = K_blade(n+d_b,n+d_b)/blade_work_percent;

F_blade = [zeros(n,1);L_b*ones(n,1)];
F_blade(n+d_b,1) = F_blade(n+d_b,1)*(blade_work_percent^4);

%Mass Matrix
M = [0.5*trace(M_blade)*eye(2)+m_r*eye(2) [cos((2*pi*(i-1))/n) -sin((2*pi*(i-1))/n);sin((2*pi*(i-1))/n) cos((2*pi*(i-1))/n)]*M_blade;M_blade*[cos((2*pi*(i-1))/n)' sin((2*pi*(i-1))/n)';-sin((2*pi*(i-1))/n)' cos((2*pi*(i-1))/n)'] M_blade];

%Stiffness Matrix
K_centripetal = -(omega^2)*M;
K_spring = [0.5*trace(K_blade)*eye(2)+k_r*eye(2) [cos((2*pi*(i-1))/n) -sin((2*pi*(i-1))/n);sin((2*pi*(i-1))/n) cos((2*pi*(i-1))/n)]*K_blade;K_blade*[cos((2*pi*(i-1))/n)' sin((2*pi*(i-1))/n)';-sin((2*pi*(i-1))/n)' cos((2*pi*(i-1))/n)'] K_blade];
%K_spring = [k_r*eye(2) zeros(2,2*n);zeros(2*n,2) K_blade];
K = K_spring+K_centripetal;

%Force Matrix
F = [eye(2) [cos((2*pi*(i-1))/n) -sin((2*pi*(i-1))/n);sin((2*pi*(i-1))/n) cos((2*pi*(i-1))/n)];[cos((2*pi*(i-1))/n)' sin((2*pi*(i-1))/n)';-sin((2*pi*(i-1))/n)' cos((2*pi*(i-1))/n)'] eye(2*n)]*[zeros(2,1);F_blade];

%--------------------------------------------------------------------------
[EVec, Eval, NatFreq, ~, ~] = MDOF_Analysis(M,K);
%--------------------------------------------------------------------------

%Time Params
t_end = 10;
del_t = 0.1;
t = linspace(0,t_end,t_end/del_t);

[V,E]= eig(K,M);

sol = zeros(size(t,2)+1,2*n+2);

for i = 1:2*n+2
    duhamels_soln = duhamels(V(:,i),V(:,i)'*M*V(:,i),F,omega,sqrt(E(i,i)),t);
    sol = sol+duhamels_soln.nu'*V(:,i)';
end

sol = sol(1:end-1,:);

if show_plots
    plot(t,sol(:,1),'r',t,sol(:,2),'k',LineWidth=2);
    legend('x_r','y_r');
    xlim([0 t_end]);
    %ylim([-0.05 0.05]);
    xlabel('Time (s)');
    ylabel('Amplitude (m)');
    title(['Response on Rotor with Defect Blade Mass Loss of ' num2str(1 - blade_work_percent),'%']);
end

function duhamels = duhamels(X,mu,F,omega,omega_calc,t)

    duhamels.time = t(1);
    duhamels.nu = 0;
    duhamels.nu_dot_del_t = 0;
    del_t = t(2)-t(1);
    t_end = t(end);

    while true
        P = F*sin(omega*duhamels.time(end));
        P_next = F*sin(omega*(duhamels.time(end)+del_t));
        phi = (X'*P)/mu;
        phi_next = (X'*P_next)/mu;
        alpha = omega_calc*del_t;

        duhamels.nu(end+1) = (cos(alpha)*duhamels.nu(end)) + ((sin(alpha)/alpha)*duhamels.nu_dot_del_t(end)) + ((del_t/alpha)^2)*(((sin(alpha)/alpha)-cos(alpha))*phi + (1-(sin(alpha)/alpha))*phi_next);
        duhamels.nu_dot_del_t(end+1) = (-alpha*sin(alpha)*duhamels.nu(end)) + (cos(alpha)*duhamels.nu_dot_del_t(end)) + ((del_t/alpha)^2)*((alpha*sin(alpha)+cos(alpha)-1)*phi + (1-cos(alpha))*phi_next);
        duhamels.time(end+1) = duhamels.time(end)+del_t;
    
        if duhamels.time(end) >= t_end
            break;
        end
    end

end
