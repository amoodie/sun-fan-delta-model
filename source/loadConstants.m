function [con] = loadConstants(conset)
%load_conset load a set of environmental constants 
%   con = load_conset(conset) loads a set of environmental constants defined by
%   the environment set by the defining string conset. Environmental constants
%   at a minimum are g (gravitational acceleration), rho_f (fluid density),
%   rho_s (sediment density), and nu (fluid kinematic viscosity).
%
%   Acceptable values for conset are:
%       - earth-quartz-water
%       - mars-quartz-water
%       - mars-basalt-water
%
%   Use type load_cons to view the constants set by each conset.

    switch conset
        case 'earth-quartz-water'
            con.g = 9.81; % gravitational constant
            con.rho_f = 1000; % fluid density, kg/m^3
            con.rho_s = 2650; % particle density, kg/m^3
            con.nu = 1.004 * 1e-6; % fluid kinematic viscosity, m^2/s
        case 'mars-quartz-water'
            con.g = 3.7; % gravitational constant
            con.rho_f = 1000; % fluid density, kg/m^3
            con.rho_s = 2650; % particle density, kg/m^3
            con.nu = 1.004 * 1e-6; % fluid kinematic viscosity, m^2/s
        case 'mars-basalt-water'
            con.g = 3.7; % gravitational constant
            con.rho_f = 1000; % fluid density, kg/m^3
            con.rho_s = 2900; % particle density, kg/m^3
            con.nu = 1.004 * 1e-6; % fluid kinematic viscosity, m^2/s
    end
    con.R = (con.rho_s - con.rho_f) / con.rho_f;
end
