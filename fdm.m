clc; clear variables, close all

format long
%% time setting
t = 1000 %integration time
dt = 0.0001; %time step
    
%% Cunstract initial domain and grid 
  n_regs = 2;
  n_el = 2;
  n_ph = 2;
  n_el_regs = [2, 2];
  n_ph_regs = [1, 1];
  reg1 = [0, 1e-9] ;
  reg2 = [1e-9, 1e-6]; % domain specification: [left, position of interface, right]
  np_regs = [5, 5]; % number of points in the domain [phase1, phase2]
  grid1 = linspace(reg1(1), reg1(2), np_regs(1)); %% grid[tstp, region]
  grid2 = linspace(reg2(1), reg2(2), np_regs(2)); %% grid[tstp, region]    
%% Initial LE ad BC
%[alpha/Beta: xC-xFe
% Beta/alpha: xC-xFe ]
xLE1 = [ 7.91841E-04, 9.99208E-01 ];
xLE2 = [ 2.83247E-02, 9.71675E-01 ];
%[alpha/Beta: uC-uFe
% Beta/alpha: uC-uFe ]
uLE1 =  xLE1(1)/xLE1(2);%, xLE1(2)/xLE1(2)];
uLE2 =  xLE2(1)/xLE2(2);%, xLE2(2)/xLE2(2)];
  
xBC1 = [ 7.91841E-04, 9.99208E-01 ];  
xBC2 = [ 7.91841E-04, 9.99208E-01 ];
xBC3 = [ 2.83247E-02, 9.71675E-01 ];
xBC4 = [ 9.23190E-03, 9.90770E-01 ];
     
uBC1 =  xBC1(1)/xBC1(2);%, xBC1(2)/xBC1(2) ]; 
uBC2 =  xBC2(1)/xBC2(2);%, xBC2(2)/xBC2(2) ]; 
uBC3 =  xBC3(1)/xBC3(2);%, xBC3(2)/xBC3(2) ];
uBC4 =  xBC4(1)/xBC4(2);%, xBC4(2)/xBC4(2) ]; 
  
%% Initial concentration profiles
    %% x-u-[tstp, region, phase, element]
x1(1,:) = grid1 * 0 + xBC1(1);%% x[tstp_1,reg_1,Phase_1:alpha,element_1:C]
x1(2,:) = grid1 * 0 + xBC1(2);%% x[tstp_1,reg_1,Phase_1:alpha,element_1:Fe]

x2(1,:) = grid2 * 0 + xBC4(1);%% x[tstp_1,reg_2,phase_2:alpha,elemet_1:C]
x2(2,:) = grid2 * 0 + xBC4(2);%% x[tstp_1,reg_2,phase_2:alpha,elemet_1:Fe]

u1(1,:) = x1(1,:) ./ x1(2,:);%% u[tstp_1,reg_1,Phase_1:alpha,element_1:C]

u2(1,:) = x2(1,:) ./ x2(2,:);%% u[tstp_1,reg_2,phase_2:alpha,elemet_1:C]


%% apply boundary condition
%% region 1
u1(1,1) = uBC1;
u1(1,end) = uBC2;
%% region 2
u2(1,1) = uBC3;
u2(1,end) =uBC4;
%% timestep loop
for tstp = 1 : t/dt;
    %% update domain, grid amd position of interface
    if tstp < 2
    else
    %% calculate u frations / initial Local equilibrium
    D = 2.5e-10 * u2 + 5.9e-13;
    dotu = (diff(u2) ./ diff( grid2 )); %D(tstp,2,2,1) .*
    j = 0;
    for i = 1 : size(dotu(1))
     j = j - D(i) * dotu(i);
     u2(i) = u2(i) - D(i) * dotu(i)
    end
    %% solve flux equations
    v = dt * j / (uBC3 - uBC2);
    %% integrate concentration values 
    
    end
end