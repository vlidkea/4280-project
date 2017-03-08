function [output, te] = steadyState (V1e,V2e)
    tank_d = 4.445; % "diameter" of the tank in cm
    tank_in = [0.635 0.635; 0.47625 0.47625]; % diameter of the intakes for each tank in cm
    tank_out = [0.47625 0.47625; 0.47625 0.47625]; %diameter of the outlets for each tank in cm

    motor_constant = 3.3;
    gravitational_constant = 981;
    pi_4 = pi()/4;

    pump = motor_constant / (tank_d^2*pi_4) * [ tank_in(1,1)^2/(tank_in(1,1)^2+tank_in(2,2)^2), tank_in(1,2)^2/(tank_in(1,2)^2+tank_in(2,1)^2);
                                                tank_in(2,1)^2/(tank_in(1,2)^2+tank_in(2,1)^2), tank_in(2,2)^2/(tank_in(1,1)^2+tank_in(2,2)^2);];
    outlet = (2*gravitational_constant)^0.5 * tank_out^2 / tank_d^2;
    
    % state space values for tank 11
    E = pump(1,1); % input
    A = outlet(1,1)^2/(2*pump(1,1)*V1e); % tank state

    % state space values for tank 12
    F = pump(1,2); % input
    B = outlet(1,2)^2/(2*pump(1,2)*V1e); % tank state

    % state space values for tank 21
    dT21c = 2*pump(1,1)*V1e + 2*pump(2,1)*V2e; % constant from linearization
    G = (pump(1,1)^2*V1e+2*pump(1,1)*pump(2,1)*V2e)/dT21c; % V1 input
    H = (pump(2,1)^2*V2e+2*pump(1,1)*pump(2,1)*V1e)/dT21c; % V2 input
    C = outlet(2,1)^2/dT21c; % tank state

    % state space values for tank 22
    dT22c = 2*pump(1,2)*V2e + 2*pump(2,2)*V1e; % constant from linearization
    I = (pump(2,2)^2*V1e+2*pump(1,2)*pump(2,2)*V2e)/dT22c; % V1 input
    J = (pump(1,2)^2*V2e+2*pump(1,2)*pump(2,2)*V1e)/dT22c; % V2 input
    D = outlet(2,2)^2/dT22c; % tank state

    te = [pump(1,1)*V1e*outlet(1,1)^-2, pump(1,2)*V2e*outlet(1,2)^-2;
          (pump(2,1)*V2e+ pump(1,1)*V1e)*outlet(2,1)^-2, (pump(2,2)*V1e+ pump(1,2)*V2e)*outlet(2,2)^-2];
    
    output = [ E*V1e/A, F*V2e/B;
              0 , 0];
    output(2,1) = ((G-E)*V1e + H*V2e + A*(output(1,1)+te(1,1)))/C;
    output(2,2) = (I*V1e+(J-F)*V2e + B*(output(1,2)+te(1,2)))/D;
    
    output = output + te;
end