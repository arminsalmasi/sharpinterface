clc; clear variables, close all

format long
%% time setting
t = 1000 %integration time
dt = 0.0001; %time step
    
%% Cunstract initial domain and grid 
  dom = [0, 1e-9, 1e-6]; % domain specification: [left, position of interface, right]
  np_reg = [5, 5]; % number of points in the domain [phase1, phase2]
  npts = sum(np_reg)
  z = [linspace(dom(1,1), dom(1,2), np_reg(1)), ... 
       linspace(dom(1,2), dom(1,3), np_reg(2))] ; % initial grid cunstruction - linear
    
%% Initial LE ad BC
%[xC-alpha,xFe-alpha
% xC_Beta ,xFe-Beta 
% xC_final ,xFe-final ] %ALPHABETICAL ORDER %Final last
xLE = [ 7.91841/10000, 9.99208/10 ...
      ; 2.83247/100  , 9.71675/10 ...
      ; 9.2319/1000  , 9.9077/10  ];

  uLE = [ xLE(1,1)/xLE(2,1), xLE(2,1)/xLE(2,1) ... 
        ; xLE(1,2)/xLE(2,2), xLE(2,2)/xLE(2,2) ]
%% Initial concentration profiles
x(1,:) = [z(1:np_reg(1,1))*0 + xLE(1,1), z(np_reg(1,1)+1: npts)*0 + xLE(3,1)]%% X(C)
x(2,:) = [z(1:np_reg(1,1))*0 + xLE(1,2), z(np_reg(1,1)+1: npts)*0 + xLE(3,2)]%% X(Fe)
u(1,:) = x(1,:)./x(2,:)%% U(C)
u(2,:) = x(2,:)./x(2,:)%% U(Fe)
%% apply boundary condition
u(1,1) = 10;
u(1,np_reg(1,1)) = 10;
u(1,np_reg(1,1)+1) = 10;
u(1,npts) = 10;

u(2,1) = 10;
u(2,np_reg(1,1)) = 10;
u(2,np_reg(1,1)+1) = 10;
u(2,npts) = 10;

%% timestep loop
for tstp = 1 : t/dt;
    %% update domain, grid amd position of interface
    if tstp < 2
    else
    %% calculate u frations / initial Local equilibrium
    %% C is 1 - alphabetical order
    %% Fe is 2 - alphabetical order
    %% phase sys is 0 
    %% phase alpha bcc is 1 - alphabetical order
    %% phase Beta fcc is 2 - alphabetical order
    % order is Var_element_phase
    
    
    for i = 2 : Npts-1 %% dU/dZ
      U_C_dotz(i) = ( U_C(i+1) -U_C(i-1) ) / (h(i+1)+h(i))  ;
    end

    DU_C_dotz = U_C_dotz(:) .* (2.5e-10 * U_C(:) + 5.9e-13); %% D.dU/dZ

    for i = 2 : Npts-1 %% d(D.dU/dZ)/dz
      DU_C_dotz_dotz(i) = ( DU_C_dotz(i+1) -DU_C_dotz(i-1) ) / (h(i+1)+h(i))  ;%% dU/dt = d(D.dU/dZ)/dz
    end


    %% solve flux equations

    j = sum(DU_C_dotz_dotz(:) .* h(:)); %/ ( u_C_1_init - u_C_inf_init)
    v =  j * dt / ( u_1_1_init - u_1_0_init);


    %% integrate concentration values 
    U_C(:) = U_C(:) + DU_C_dotz_dotz(:) * dt
    UKEEP(cnt,:) =  U_C(:);
    cnt
    U_C(1) = bc_C_1_0;
    U_Fe(1) = bc_Fe_1_0;
    U_C(npts(1)) = bc_C_1_1;
    U_Fe(npts(1)) = bc_Fe_1_1;
    U_C(npts(1)+1) = bc_C_2_0;
    U_Fe(npts(1)+1) = bc_Fe_2_0;
    U_C(Npts) = bc_C_2_1;
    U_Fe(Npts) = bc_Fe_2_1;
    %% finalize timestep

    %% check if mass concervation is violated (only if there is impingment or not!)

    %% check convergance?! or just simply jump to the next time step 
    end
end