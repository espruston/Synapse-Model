function par = model_parameters_det(par_init)

%%%%%%%%%%%%%%%%%%%%%
%%%%%%This script contains parameters. Not all are used, since there are
%%%%%%leftovers other models.
%%%%%%%%%%%%%%%%%%%%%

par = zeros(45,1);

num_sim = 200;

%Fixed assumptions
sigm = 76.5154; %Scale parameter in vesicles distribution

act_const_init = 1e6;

% Calcium
CaMax_rest = par_init(20); %M

z_dist = 0.01; %z distance of calcium reading
grid_points = 71;
hill = 1;

num_stim = 2; %Number of stimulations
height = 1; %Height of simulation volume
min_step = 1e-6; %Minimum time step in stochastic simulation
size_of_mini = 0.6e-9;


AZ_size = par_init(24); %0.6239936794; %Volume is adjusted to match the box volume with 0.553; %?m

%Fast sensor
n_max = 5; %Cooperativity of fast sensor
k_3 = 1.4e8; %(1/M*s)
k_min3 = 4000; %(1/s)
b_f = 0.5; %Cooperativity factor
l_0 = 3.50e-4; %Basal fusion rate
f_val = 27.9780; % = (6000/3.5e-4)^(1/5)


%Parameters not used anymore
Ca_rest = NaN; %Not used anymore
k_rep_constraint = 1000; %Not used anymore
num_rep = 1; %Not used anymore
pool_size = 180*num_rep; %Not used anymore



%Basic
CalC_geometry = par_init(1);
Q_max = par_init(2);
k_M_calc = par_init(3);
k_M_rest = par_init(19); %mM
num_ves = par_init(4);
k_rep = par_init(5);
%Model 2
m_max = par_init(6);
SS_PM = par_init(7);
k_d_second = par_init(8);
k_4 = par_init(9);
s_val = par_init(10);
b_s = (m_max ~= 0) * 0.5;
k_min4 = k_4 * k_d_second;
%Model 3
act_model_type = par_init(11);
k_M_act = (act_model_type ~= 0) * par_init(12)^act_model_type;
delay_rate = par_init(13);
invdelay_rate = par_init(14);
act_rate_const = (act_model_type ~=0) * (act_const_init)^act_model_type;
inact_rate_const = act_rate_const*k_M_act;
%Model 4
Ca_prim_type = par_init(15);
prim_kM = (Ca_prim_type ~= 0)*par_init(16)^Ca_prim_type;
prim_rate_const = par_init(17);
unprim_rate_const = par_init(18);
unprim_rate_const_0 = par_init(21);
%Model 5
u_val = par_init(22);

model_type = par_init(25);
unprim_onestate = par_init(26);

%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par(1) = k_3;
par(2) = k_min3;
par(3) = b_f;
par(4) = l_0;
par(5) = f_val;
par(6) = k_4;
par(7) = k_min4;
par(8) = b_s;
par(9) = s_val;
par(10) = k_rep;
par(11) = k_rep_constraint;
par(12) = pool_size;
par(13) = n_max;
par(14) = m_max;
par(15) = num_rep;
par(16) = AZ_size;
par(17) = Q_max;
par(18) = k_M_calc;
par(19) = size_of_mini;
par(20) = sigm;
par(21) = Ca_rest;
par(22) = num_ves;
par(23) = CaMax_rest;
par(24) = k_M_rest;
par(25) = k_M_act;
par(26) = act_rate_const;
par(27) = inact_rate_const;
par(28) = delay_rate;
par(29) = invdelay_rate;
par(30) = z_dist;
par(31) = prim_kM;
par(32) = prim_rate_const;
par(33) = unprim_rate_const;
par(34) = grid_points;
par(35) = num_stim;
par(36) = height;
par(37) = min_step;
par(38) = CalC_geometry;
par(39) = Ca_prim_type;
par(40) = SS_PM;
par(41) = hill;
par(42) = act_model_type;
par(43) = u_val;
par(44) = unprim_rate_const_0;
par(45) = num_sim;
par(46) = model_type;
par(47) = unprim_onestate;

end