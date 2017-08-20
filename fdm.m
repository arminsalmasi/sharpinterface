clc; clear variables, close all

format long

% time setting
t = 3600; %integration time
dt = 0.1; %time step
    
%% timestep loop
for tstp = 1 : floor(t/dt);   
    if tstp < 2
      % Initiate concentration
      [er,u1,u2,grid1,grid2, reg1,reg2, np_regs, uBC1] = do_init();
      %plot( [grid1,grid2], [u1(1,:),u2(1,:)] );
      
    else 
      % solve PDE to calculate flux
      D21 = 2.5e-10 * u2(1,:) + 5.9e-13; % D reg=2 el=1
      [j] = do_calc_flux(u2(1,:), grid2, D21);
    
      % solve flux equations, calculate intevace speed and displacement
      v = j(1) / (uBC1(3) - uBC1(2));
      dx = v * dt;
      
      % update concentration values 
      u2_new = u2;
      for i = 2 : size(u2,2)
        u2_new(1,i) = u2(1,i) + j(i-1); 
      end
      %hold on
      %% reappply boundry condition
      %% update grid
      if reg1(2)+dx>=reg2(2) 
          er = 'larger than doamin'
          tstp
          tstp*dt
          %break
      else
        grid1_new = linspace(reg1(1), reg1(2)+dx, np_regs(1)); %% grid[tstp, region]
        grid2_new = linspace(reg2(1)+dx, reg2(2), np_regs(2)); %% grid[tstp, region]    
      end
      hold on
      if (tstp==2 || dt*tstp==t )%mod(t,tstp*dt*100)==0)
        tstp
        plot( [grid1,grid2], [u1(1,:),u2(1,:)] );
        plot( [grid1_new,grid2_new], [u1(1,:),u2_new(1,:)] );
      end
      grid1 = grid1_new;
      grid2 = grid2_new;
      u2 = u2_new;

    end
    %% update keepers
end
%% EOF
%%
function [j]=do_calc_flux(u,grd,D)
    [dotu] = do_fdm_diff(u,grd);
    for i = 1 : size(dotu,2)
        j(i) = -1 * D(i) * dotu(i);
    end
    return
end
%%
function [dotf_z] = do_fdm_diff(f,z)
%fdm first derivitive
    for i = 1 : size(f,2)-1
        dotf_z(i) = (f(1,i+1)- f(1,i)) / (z(i+1) + z(i)) ; 
    end
    return
end
%%
function [u1,u2]=do_apply_BC(uBC,u1,u2)
    %% region 1
    u1(1,1) = uBC(1);
    u1(1,end) = uBC(2);
    %% region 2
    u2(1,1) = uBC(3);
    u2(1,end) =uBC(4);
end
%%  
function [er, u1, u2, grid1, grid2, reg1,reg2, np_regs, uBC1] = do_init()
    format long
    er = true;
    n_regs = 2;
    n_el = 2;
    n_ph = 2;
    n_el_regs = [2, 2];
    n_ph_regs = [1, 1];
    reg1 = [0, 1e-9] ;
    reg2 = [1e-9, 5e-5]; % domain specification: [left, position of interface, right]
    np_regs = [10, 1000]; % number of points in the domain [phase1, phase2]
    grid1 = linspace(reg1(1), reg1(2), np_regs(1)); %% grid[tstp, region]
    grid2 = linspace(reg2(1), reg2(2), np_regs(2)); %% grid[tstp, region]    
    % Initial LE ad BC
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
    
    uBC1 =  [ xBC1(1)/xBC1(2), ...%, xBC1(2)/xBC1(2) ]; 
              xBC2(1)/xBC2(2), ...%, xBC2(2)/xBC2(2) ]; 
              xBC3(1)/xBC3(2), ...%, xBC3(2)/xBC3(2) ];
              xBC4(1)/xBC4(2)];%, xBC4(2)/xBC4(2) ]; 

    % Initial concentration profiles
    % x-u-[tstp, region, phase, element]
    x1(1,:) = grid1 * 0 + xBC1(1);%% x_reg_1[Phase_1:alpha,element_1:C]
    x1(2,:) = grid1 * 0 + xBC1(2);%% x_reg_1[Phase_1:alpha,element_1:Fe]
    x2(1,:) = grid2 * 0 + xBC4(1);%% x_reg_2[phase_2:alpha,elemet_1:C]
    x2(2,:) = grid2 * 0 + xBC4(2);%% x_reg_2[phase_2:alpha,elemet_1:Fe]
    u1(1,:) = x1(1,:) ./ x1(2,:);%% u_reg_1[Phase_1:alpha,element_1:C]
    u1(2,:) = x1(2,:) ./ x1(2,:);%% u_reg_1[Phase_1:alpha,element_1:Fe]
    u2(1,:) = x2(1,:) ./ x2(2,:);%% u_reg_2[phase_2:alpha,elemet_1:C]
    u2(2,:) = x2(2,:) ./ x2(2,:);%% u_reg_2[Phase_1:alpha,element_1:Fe]
    % apply boundary condition
    [u1,u2]=do_apply_BC(uBC1,u1,u2);
end
%%
