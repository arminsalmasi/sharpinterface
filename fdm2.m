clc; clear variables, close all

%% time setting
t = 5200; %integration time
dt = 0.001; %time step
plot_t = [1, 100, 200, 500, 1000, 3600, 5200];

%% domain setting  
rgs_geo = [0   , 1e-9 ...    % geometery of region 1
          ;1e-9, 5e-5 ];     % geometery of region 2
np_rgs = [10, 1000];   % number of points in each region
grid_types = [1, 1]; % linear in both regions;
geo_facts = [1,1]; % both regions linear / default
grid = do_make_grid(rgs_geo , np_rgs , grid_types , geo_facts);  

%% sunstitutional interstitial
Sub_Int = [ 0, 1
          ; 0, 1 ]; 

%% Local equilibrium values from TC (mole fraction)
xLE = [ 7.91841E-04, 9.99208E-01 ...
      ; 2.83247E-02, 9.71675E-01 ];

uLE = do_calc_ufraction(Sub_Int, xLE)

%% time loop
for tstp = 1 : floor(t/dt);   
    tstp_keep(tstp)=tstp*dt;
    if tstp < 2
      % Initiate concentration
      [er,u1,u2,grid1,grid2, reg1,reg2, np_rgs, uBC1] = do_init();
      % update keepers
      grid1_keep(tstp,:) = grid1;
      grid2_kepp(tstp,:) = grid2;
      u2_keep(tstp,:,:) = u2;
    else 
      % solve PDE to calculate flux
      D21 = 2.5e-10 * u2(1,:) + 5.9e-13; % D reg=2 el=1
      [j21] = do_calc_flux(u2(1,:), grid2, D21);
    
      % solve flux equations, calculate intevace speed and displacement
      v = j21(1) / (uBC1(3) - uBC1(2));
      dx = v * dt;
      % update concentration values 
      for i = 2 : size(u2,2)
        u2(1,i) = u2(1,i) + j21(i-1); 
      end
      %hold on
      % update grid
      if reg1(2)+dx>reg2(2) 
          er = 'larger than doamin'
          tstp
          tstp*dt
          plot_t(1,end)= tstp*dt;
          break
      else
        % update keepers
        grid1(:) = linspace(reg1(1), reg1(2)+dx, np_rgs(1)); %% grid[tstp, region]
        grid2(:) = linspace(reg2(1)+dx, reg2(2), np_rgs(2)); %% grid[tstp, region]    
        
      end
      % update keepers
      v_keep(tstp) = v;
      grid1_keep(tstp,:) = grid1;
      grid2_keep(tstp,:) = grid2;
      u2_keep(tstp,:,:) = u2;
      u1_keep(tstp,:,:) = u1;
    end

end
% Printing
figure
hold on
  for i = 1: size(plot_t,2)
    t_plot =floor(plot_t(1,i)/dt)
    grid1_plot(:) = grid1_keep(t_plot,:);
    grid2_plot(:) = grid2_keep(t_plot,:);
    u1_plot(:,:) = u1_keep(t_plot,:,:);     
    u2_plot(:,:) = u2_keep(t_plot,:,:);     
    plot( [grid1_plot,grid2_plot], [u1_plot(1,:),u2_plot(1,:)] );
  end
figure
hold on
plot(tstp_kepp,v_keep);


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
function [u1, u2, grid1, grid2, reg1,reg2, np_regs, uBC1, er] = do_init(reg1, reg2, np_regs,  xLE1, xLE2)
    er = true;
     
   
    %% calculate u-fractions
    %[alpha/Beta: uC-uFe
    % Beta/alpha: uC-uFe ]
    uLE1 =  xLE1(1)/xLE1(2);%, xLE1(2)/xLE1(2)];
    uLE2 =  xLE2(1)/xLE2(2);%, xLE2(2)/xLE2(2)];
    
    %% Boundary condition
    xBC1 = [ 7.91841E-04, 9.99208E-01 ];  
    xBC2 = [ 7.91841E-04, 9.99208E-01 ];
    xBC3 = [ 2.83247E-02, 9.71675E-01 ];
    xBC4 = [ 9.23190E-03, 9.90770E-01 ];
    
    %% Boundary condition : u-fraction
    uBC1 =  [ xBC1(1)/xBC1(2), ...%, xBC1(2)/xBC1(2) ]; 
              xBC2(1)/xBC2(2), ...%, xBC2(2)/xBC2(2) ]; 
              xBC3(1)/xBC3(2), ...%, xBC3(2)/xBC3(2) ];
              xBC4(1)/xBC4(2)];%, xBC4(2)/xBC4(2) ]; 

    %% Initiate concentration profiles - all regions all elements
    %% x-u-[tstp, region, phase, element]
    x1(1,:) = grid1 * 0 + xBC1(1);%% x_reg_1[Phase_1:alpha,element_1:C]
    x1(2,:) = grid1 * 0 + xBC1(2);%% x_reg_1[Phase_1:alpha,element_1:Fe]
    x2(1,:) = grid2 * 0 + xBC4(1);%% x_reg_2[phase_2:alpha,elemet_1:C]
    x2(2,:) = grid2 * 0 + xBC4(2);%% x_reg_2[phase_2:alpha,elemet_1:Fe]
    u1(1,:) = x1(1,:) ./ x1(2,:);%% u_reg_1[Phase_1:alpha,element_1:C]
    u1(2,:) = x1(2,:) ./ x1(2,:);%% u_reg_1[Phase_1:alpha,element_1:Fe]
    u2(1,:) = x2(1,:) ./ x2(2,:);%% u_reg_2[phase_2:alpha,elemet_1:C]
    u2(2,:) = x2(2,:) ./ x2(2,:);%% u_reg_2[Phase_1:alpha,element_1:Fe]
    
    %% apply boundary condition
    [u1,u2]=do_apply_BC(uBC1,u1,u2);
end
%%
