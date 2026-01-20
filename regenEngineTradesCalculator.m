%% Regen Engine Trades
% Units start in SI, converted to Imperial later for front end (MATLAB app)
% -------------------------------------------------------------------------

clc
clear
% 1) Make sure python is installed on your system and MATLAB can use it
% 2) Install the rocketcea library/module. Ask internet or AI for help. May
%   need to do some testing to get it working first before trying this
%   script
% 3) Need this py.import... line in order to load the module every time.
%   Like you would for a normal python script.
py.importlib.import_module("rocketcea");
% For data export, to be graphed later in MATLAB app
filename_datastorage = append(pwd, "\regenEngineTradesData.mat");

%% Target Performance Parameters
% Derived from requirements
% c* efficiency
% Maximum allowable cost of consumed propellants
% Burntime

eta_cstar_min = .6;
cost_max = 200; % USD
% Normalized now to 1s because the 10s burntime requirement was removed
%   (01/10/2026)
burntime = 5; % s

%% Constants
% g0                        Accel. due to gravity
% Lstar                     L* for props: liquid eth and liquid nit
% CR                        Contraction ratio
% R_universal               Universal gas constant
% pa                        Ambient pressure
% DF*                       Design factor *(for multiple things)
% T_AlSi10Mg_melt           Chamber alloy melting temp
% cost_nit                  USD per kgm
% cost_eth                  USD per kgm
% yieldstress_alloy*        Linearly interpolated with separate function
%                           (dependent on hot wall temperature)
% E_alloy                   Modulus of elasticity of chamber alloy
% alpha_alloy               CTE
% lchannel                  Conservative estimate of coolant channel length
%                               (overestimate- leads to overestimating
%                               minimum required coolant pressure)

g0 = 9.81;
Lstar = .6; % estimate for MIT's value for their eth/nit biprop was .6
CR = 6;
R_universal = 8.314;
% Calculate pa at Rolla's altitude
% psi to pa
p_sea = 14.7 .* 6894.76;
% K
T_sea = 288.15;
% kg/mol
mlr_wgt_air = .02896968;
L = .00976;
% Thought I found this on wikipedia the other day but now I can't find it
fPAtAlt = @(alt) p_sea .* (1 + (L .* alt) ./ T_sea).^(-g0 .* mlr_wgt_air ./ (R_universal .* L));
alt = 342;
pa = fPAtAlt(alt);

% RDT technically has no DF for these
DFstress = 1.5;
DFtemp = 1.2;
DFcoolantpress = 1.2;
% g/cc to kg/m^3
densalloy = 2.665 .* 10.^-3;
% deg C to K
T_AlSi10Mg_melt = 570 + 273.15;
% USD/lbm to USD/kgm
cost_nit = 6 .* 2.20462;
cost_eth = .46 .* 2.20462;
% MPa to Pa (75 +- 10) - heat treated would be 80, 70 +- 10
E_alloy = 85E9;
% From data sheet: https://www.eos.info/var/assets/03_system-related-assets/material-related-contents/metal-materials-and-examples/metal-material-datasheet/aluminium/material_datasheet_eos_aluminium-alsi10mg_en_web.pdf
% Took most conservative/highest estimate
% K^-1
alpha_alloy = 27E-6;
% in to m
% l_channel = 5 .* .0254;
numelsnoz = 10;
numsectionsnoz = 5;

%% Independent, Variable Parameters
% pc
% OF                        O to F ratio
% expansion_ratio
% mdot
% d_channel
% num_channels
% T_coolant
% wallt                     hot-wall thickness
% k_wall                    wall material thermal conductivity

list_var_names = ["pc";
                  "OF";
                  "expansion_ratio";
                  "mdot";
                  "d_channel";
                  "num_channels";
                  "T_coolant";
                  "wallt";
                  "k_wall"];

%% Some precalcs
% Applying DF's to hot wall temp
T_w_max = T_AlSi10Mg_melt ./ DFtemp;

%% Make iteration ranges for independent vars
stdvarrange = [1];
highdefvarrange = linspace(.8, 1.2, 5);
stdvarrange = highdefvarrange;
% psi to Pa
pc_nom = 160 .* 6894.76;
OF_nom = 4.4;
expansion_ratio_nom = 1.9;
% lbm/s -> kgm/s
mdot_nom = 1 .* .45359237;
% mm to m
d_channel_nom = 1 .* 1E-3;
num_channels_nom = 60;
% convert deg C to K LATER
T_coolant_nom = 30;
% mm to m
wallt_nom = 1 ./ 1000;
% W/m-K - Largest range I found in the paper: On the thermal conductivity of AlSi10Mg and lattice structures made by laser powder bed fusion Richard R.J. Sélo, Sam Catchpole-Smith, Ian Maskery, Ian Ashcroft, Christopher Tuck
% k_wall range is from ~98 to ~183
k_wall_nom = (98 + 183) / 2;

pc_range = pc_nom .* stdvarrange;
% OF_range = OF_nom .* stdvarrange;
OF_range = linspace(1.5, 1.7, 5);
% expansion_ratio_range = expansion_ratio_nom .* stdvarrange;
expansion_ratio_range = linspace(1.8, 2.2, 5);
mdot_range = mdot_nom .* stdvarrange;
d_channel_range = d_channel_nom .* stdvarrange;
num_channels_range = num_channels_nom .* stdvarrange;
T_coolant_range = T_coolant_nom .* stdvarrange;
% deg C to K
T_coolant_range = T_coolant_range + 273.15;
wallt_range = wallt_nom .* stdvarrange;
k_wall_range = k_wall_nom .* stdvarrange;
% For sizing resultant nd arrays to be compatible in one data structure -
% allows for efficient export and reading later
ranges = {pc_range, OF_range, expansion_ratio_range, mdot_range, d_channel_range, num_channels_range, T_coolant_range, wallt_range, k_wall_range};
ndsvar = cellfun(@length, ranges);
[mdot_range, d_channel_range, num_channels_range, T_coolant_range, wallt_range, k_walls] = ndgrid(mdot_range, d_channel_range, num_channels_range, T_coolant_range, wallt_range, k_wall_range);
% Channels are pentagon shape, square with a 45-45-90 triangle on one side
fChannelArea = @(d) d.^2 + (d ./ sqrt(2)).^2 .* .5;
A_channel_range = fChannelArea(d_channel_range);
% + 3 comes from pc, OF, and expansion_ratio, which are iterated through
% via for loops - NOT vectorized like the other indep. vars. This is
% because CEA is not vectorizeable, and those are required CEA inputs
nds = ndims(mdot_range) + 3;
% Max possible mdot according to prop_cost and burntime limits
mdot_max = cost_max ./ burntime ./ (cost_nit .* (1 ./ (1 + OF_range) .* OF_range) + cost_eth .* (1 ./ (1 + OF_range)));

% For CEA later:
fuel_name = "ETHANOL";
ox_name = "N2O";

%% Things for film coefficient calcs later
% Coolant is ethanol
mlr_wgt_eth = .0460684; % kg/mol
% -------------------------------------------------------------------------
% Stuff for coolant density calcs
% For liquid ethanol cp(T) info - https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=187F
% From Green J.H.S 1961, found at temp of 298.15 K ->
% J/mol-K to J/kg-K
cp_coolant_avg = 111.96 ./ mlr_wgt_eth;
Tcrit_eth = 514; % K
Pcrit_eth = 6.14E6; % Pa
% liquid eth compressibility factor
Zc_eth = .241;
Vc_eth = Zc_eth .* R_universal .* Tcrit_eth ./ Pcrit_eth;
rho_crit_eth = 276;
fTr = @(T, Tcrit) T ./ Tcrit;
fRackett = @(rho_crit, Zc, T, Tcrit) rho_crit ./ Zc.^(1 - fTr(T, Tcrit));
% -------------------------------------------------------------------------
% For coolant thermal conductivity calcs later - all these functions are
% just meant to organize the summation expressions and other things from
% the papers more easily (read the papers)
% From Assael, M.J., et al., Reference Data for the Thermal Conductivity of Ethanol , J. Phys. Chem. Ref. Data, 2010
flambda0_eth = @(T) (-2.09575 + 19.9045 .* fTr(T, Tcrit_eth) - 53.964 .* fTr(T, Tcrit_eth).^2 + 82.1223 .* fTr(T, Tcrit_eth).^3 - 1.98864 .* fTr(T, Tcrit_eth).^4 - .495513 .* fTr(T, Tcrit_eth).^5) ./ (.17223 - .078273 .* fTr(T, Tcrit_eth) + fTr(T, Tcrit_eth).^2) .* 10.^-3; % W/m-K
B1i_eth = [2.67222E-2, 1.48279E-1, -1.30429E-1, 3.46232E-2, -2.44293E-3];
B2i_eth = [1.77166E-2, -8.93088E-2, 6.84664E-2, -1.45702E-2, 8.09189E-4];
is_deltalambda = 1:5;
fRshpSmmtnArrys = @(arry, nds) reshape(arry, [ones(1, nds) numel(arry)]);
B1i_eth = fRshpSmmtnArrys(B1i_eth, nds);
B2i_eth = fRshpSmmtnArrys(B2i_eth, nds);
is_deltalambda = fRshpSmmtnArrys(is_deltalambda, nds);
frho_r = @(rho, rho_crit) rho ./ rho_crit;
fdeltalambda_eth = @(rho, T) sum((B1i_eth + B2i_eth .* fTr(T, Tcrit_eth)) .* (frho_r(rho, rho_crit_eth)).^is_deltalambda, nds + 1);
fdeltalambdac_eth = @(rho, T) 1.7E-3 ./ (7E-2 + abs(fTr(T, Tcrit_eth) - 1)) .* exp(-(1.7 .* (frho_r(rho, rho_crit_eth) - 1)).^2);
fk_eth = @(rho, T) flambda0_eth(T) + fdeltalambda_eth(rho, T) + fdeltalambdac_eth(rho, T); % inputs: T in K, rho in kg/m^3; output in W/m-K
% k_coolant sanity check
T_exp1 = 300;
rho_exp1 = 850;
k_exp1_act = 209.68E-3;
k_exp1_calc = fk_eth(rho_exp1, T_exp1);
% -------------------------------------------------------------------------
% For coolant viscosity calcs - read the paper
% From Sotiriadou, S., Ntonti, E., Velliadou, D., Antoniadis, K. D., Assael, M. J., & Huber, M. L. (2023). Reference Correlation for the Viscosity of Ethanol from the Triple Point to 620 K and Pressures Up to 102 MPa
bi_eth = [.422373, -3.78868, 23.8708, -7.89204, 2.09783, -.247702];
ci_eth = [-.0281703, 1];
is_eta0b = 0:5;
is_eta0c = 0:1;
bi_eth = fRshpSmmtnArrys(bi_eth, nds);
ci_eth = fRshpSmmtnArrys(ci_eth, nds);
is_eta0b = fRshpSmmtnArrys(is_eta0b, nds);
is_eta0c = fRshpSmmtnArrys(is_eta0c, nds);
feta0_eth = @(T) sum(bi_eth .* fTr(T, Tcrit_eth).^is_eta0b, nds + 1) ./ sum(ci_eth .* fTr(T, Tcrit_eth).^is_eta0c, nds + 1);
epsilon_k_B_eth = 265;
fTstar = @(T, epsilon_k_B) T ./ epsilon_k_B;
di = [-1.9572881E1, 2.1973999E2, -1.0153226E3, 2.4710125E3, -3.3751717E3, 2.4916597E3, -7.8726086E2];
di = fRshpSmmtnArrys(di, nds);
d7 = 1.4085455E1;
d8 = -3.4664158E-1;
is_Bstar_eta = 0:6;
is_Bstar_eta = fRshpSmmtnArrys(is_Bstar_eta, nds);
Bstar_eta_eth = @(Tstar) sum(di .* Tstar.^(-.25 .* is_Bstar_eta), nds + 1) + d7 .* Tstar.^-2.5 + d8 .* Tstar.^-5.5;
NA = 6.02214076E23;
% nm to m
sigma_eth = .479E-9;
fB_eta_eth = @(T) Bstar_eta_eth(fTstar(T, epsilon_k_B_eth)) .* NA .* sigma_eth.^3 ./ mlr_wgt_eth;
feta1_eth = @(T) feta0_eth(T) .* fB_eta_eth(T);
fdeltaeta_eth = @(rho, T) frho_r(rho, rho_crit_eth).^(2./3) .* fTr(T, Tcrit_eth).^.5 .* (8.32575272 .* frho_r(rho, rho_crit_eth) + 9.66535242E-2 .* (frho_r(rho, rho_crit_eth).^8 ./ (fTr(T, Tcrit_eth).^4 .* (1 + frho_r(rho, rho_crit_eth).^2) - fTr(T, Tcrit_eth).^2)));
fmu_eth = @(T, rho) (feta0_eth(T) + feta1_eth(T) .* rho + fdeltaeta_eth(rho, T)) .* 10.^-6; % inputs: T in K, rho in kg/m^3; output in Pa-s
% mu_coolant sanity check
T_exp1 = 300;
rho_exp1 = 10;
mu_exp1_act = 8.9382E-6;
mu_exp1_calc = fmu_eth(T_exp1, rho_exp1);

% Results of interest, things to graph in the app
thermstress = zeros(ndsvar);
T_coolant_f = zeros(ndsvar);
Twg = zeros(ndsvar);
yieldstress_alloy = zeros(ndsvar);
Twl = zeros(ndsvar);
q = zeros(ndsvar);
Isp = zeros(ndsvar);
prop_cost = zeros(ndsvar);
At = zeros(ndsvar);
cstar = zeros(ndsvar);
cstar_theo = zeros(ndsvar);
thrust = zeros(ndsvar);
Ae = zeros(ndsvar);
eta_cstar = zeros(ndsvar);
throat_flow_temp = zeros(ndsvar);
Vc = zeros(ndsvar);
flow_sonic = zeros(ndsvar);
vol_engine = zeros(ndsvar);
r_engine = zeros([ndsvar, numelsnoz .* numsectionsnoz]);
z_engine = zeros([ndsvar, numelsnoz .* numsectionsnoz]);
L_nozzle_parabolic = zeros(ndsvar);
Re = zeros(ndsvar);
thetaN = zeros(ndsvar);
xN = zeros(ndsvar);
R1 = zeros(ndsvar);
R1p = zeros(ndsvar);
alpha = zeros(ndsvar);
L_chamber_circular_narrow = zeros(ndsvar);
Rc = zeros(ndsvar);
L_chamber_linear = zeros(ndsvar);
TWR = zeros(ndsvar);
hoop_stress_doghouse = zeros(ndsvar);
hoop_stress_fins = zeros(ndsvar);


%% Begin calcs
% Workflow - hybrid iterative vectorized method
    % Vary pc
    %   Vary OF
    %       Vary over expansion_ratio
    %           Get cstar_theo (do CEA)
    %           Begin vectorized equations
    %               Get post-throat gamma and other things
    %               Get exit press
    %               Get CF
    %               Get throat temp from CEA
    %               Get required At
    %               Get required Ae
    %               Get Thrust
    %               Already know burntime --> get prop_cost
    %               Get Isp
    %               Get cstar
    %               Get eta_cstar
    %               Get throat wall temp
    % Get thermal stress at throat on hot wall for all calc'd engines
    % Get final coolant pressure
    % Get required coolant pressure to keep it liquid
    % Record all data in .mat file and load in app for user
    %   friendly interaction
    % -> to app

% Start runtime tracking for engine iterations
runtime = tic;
% Condensed engine iteration structure brancher
% CEA syntax:
% py.rocketcea.cea_obj.CEA_Obj(fuelName, oxName) - see input_cards.py
% CEA_Obj class functions: (Pc, MR, eps); Pc in psi, MR is OF, eps is expansion_ratio
cea_out = py.rocketcea.cea_obj.CEA_Obj(fuelName = fuel_name, oxName = ox_name, fac_CR = CR);
%% Vary pc
for i_pc = 1:length(pc_range)
    pc = pc_range(i_pc);

    %% Vary OF
    for i_OF = 1:length(OF_range)
        OF = OF_range(i_OF);

        %% Vary expansion ratio
        for i_eps = 1:length(expansion_ratio_range)
            expansion_ratio = expansion_ratio_range(i_eps);

            %% CEA Calcs
            % Needed for basic results
            % Remember to convert Pa to psi (for inputs), and ft/s to m/s
            % (for outputs)
            % Get specific heat ratio and flow molar weight at throat - for
            % later use in thrust coeff and sonic throat flow velocity calc
            throat_flow_props = double(cea_out.get_Throat_MolWt_gamma(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio));
            throat_flow_mlr_wgt = throat_flow_props(1) .* 10.^-3; % lbm/lbm-mol -> kg/mol
            throat_flow_gamma = throat_flow_props(2);
            exit_flow_props = double(cea_out.get_exit_MolWt_gamma(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio));
            exit_flow_gamma = exit_flow_props(2);
            % Needed for thrust calcs later
            pex = pc ./ cea_out.get_PcOvPe(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio);
            Cf = sqrt(2 .* exit_flow_gamma.^2 ./ (exit_flow_gamma - 1) .* (2 ./ (exit_flow_gamma + 1)).^((exit_flow_gamma + 1) ./ (exit_flow_gamma - 1)).*(1 - (pex ./ pc).^((exit_flow_gamma - 1) ./ exit_flow_gamma))) + expansion_ratio .* ((pex - pa) ./ pc);
            % Thought theoretical Isp might be interesting to compare to
            % predicted-actual, not really used atm though (ca 11/04/2025)
            Isp_theo = cea_out.get_Isp(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio);
            % Gets throat temp needed for several calculations later
            % deg R to K
            CEA_temps = double(cea_out.get_Temperatures(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio)) .* 5 ./ 9;
            throat_flow_temp(i_pc, i_OF, i_eps, :, :, :, :, :, :) = CEA_temps(2);
            chamber_flow_temp = CEA_temps(1);
            exit_flow_temp = CEA_temps(3);
            % Fix units from CEA output lbm/ft^3 -> kg/m^3
            CEA_densities = double(cea_out.get_Densities(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio)) .* 16.018463;
            throat_flow_density = CEA_densities(2);
            % Gets things needed for gas film coeff calc later
            throat_trans = double(cea_out.get_Throat_Transport(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio));
            % milipoise to Pa-s
            mu_flow = throat_trans(2) .* 10.^-4;
            % millicalories / (cm-K-sec) to W/m-K
            k_flow = throat_trans(3) .* .4184;
            Pr_flow = throat_trans(4);
            % ft/s -> m/s
            cstar_theo(i_pc, i_OF, i_eps, :, :, :, :, :, :) = cea_out.get_Cstar(Pc = pc * 0.000145038, MR = OF) .* .3048;

            % Get some KPPs and other sizing info for trade alternatives
            flow_sonic(i_pc, i_OF, i_eps, :, :, :, :, :, :) = sqrt(throat_flow_gamma * R_universal / throat_flow_mlr_wgt * squeeze(throat_flow_temp(i_pc, i_OF, i_eps, :, :, :, :, :, :)));
            % Check sonic flow velocity
            CEA_velocity_sonic_flow = cea_out.get_Chamber_SonicVel(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio);
            At(i_pc, i_OF, i_eps, :, :, :, :, :, :) = mdot_range ./ squeeze(flow_sonic(i_pc, i_OF, i_eps, :, :, :, :, :, :)) ./ throat_flow_density;
            Vc(i_pc, i_OF, i_eps, :, :, :, :, :, :) = squeeze(At(i_pc, i_OF, i_eps, :, :, :, :, :, :)) .* Lstar;
            Ae(i_pc, i_OF, i_eps, :, :, :, :, :, :) = squeeze(At(i_pc, i_OF, i_eps, :, :, :, :, :, :)) .* expansion_ratio;
            thrust(i_pc, i_OF, i_eps, :, :, :, :, :, :) = Cf .* squeeze(At(i_pc, i_OF, i_eps, :, :, :, :, :, :)) .* pc;
            prop_cost(i_pc, i_OF, i_eps, :, :, :, :, :, :) = burntime .* mdot_range .* (cost_nit .* (1 ./ (1 + OF) .* OF) + cost_eth .* (1 ./ (1 + OF)));
            Isp(i_pc, i_OF, i_eps, :, :, :, :, :, :) = squeeze(thrust(i_pc, i_OF, i_eps, :, :, :, :, :, :)) ./ mdot_range ./ g0;
            cstar(i_pc, i_OF, i_eps, :, :, :, :, :, :) = squeeze(Isp(i_pc, i_OF, i_eps, :, :, :, :, :, :)) .* g0 ./ Cf;
            eta_cstar(i_pc, i_OF, i_eps, :, :, :, :, :, :) = squeeze(cstar(i_pc, i_OF, i_eps, :, :, :, :, :, :)) ./ squeeze(cstar_theo(i_pc, i_OF, i_eps, :, :, :, :, :, :));

            % -------------------------------------------------------------
            % Needed to find film coeffs
            % flow velocity in chamber (just throat flow vel or sonic flow
            % vel bc we're at throat)
            flow_vel = flow_sonic;
            % hydraulic radius of channel is 4 * area / perimeter
            dH_channels = 4 .* A_channel_range ./ (3 .* d_channel_range + 2 .* d_channel_range ./ sqrt(2));
            % Remember it's just fuel in channels
            mdot_channel = mdot_range .* (1 ./ (1 + OF)) ./ num_channels_range;
            % rackett equation used to calc coolant density
            rho_coolant = fRackett(rho_crit_eth, Zc_eth, T_coolant_range, Tcrit_eth);
            k_coolant = fk_eth(rho_coolant, T_coolant_range);
            mu_coolant = fmu_eth(T_coolant_range, rho_coolant);
            fluid_vel = mdot_channel ./ A_channel_range ./ rho_coolant;
            % -------------------------------------------------------------

            % liquid and gas film coeffs - equas 8-21 and 8-24 in RPE by
            % Sutton
            % Check throat Re is high enough (> 10k) to use Sutton hl
            % equation
            % hydraulic diam is 4 * area / perimeter
            dt = sqrt(squeeze(At(i_pc, i_OF, i_eps, :, :, :, :, :, :)) .* 4 ./ pi);
            Ret = throat_flow_density .* squeeze(flow_vel(i_pc, i_OF, i_eps, :, :, :, :, :, :)) .* dt ./ mu_flow;
            Re_coolant = dH_channels .* fluid_vel .* rho_coolant ./ mu_coolant;
            Pr_coolant = cp_coolant_avg .* mu_coolant ./ k_coolant;
            % Exponent on Pr_coolant is .4 bc fluid is being heated - https://www.sciencedirect.com/topics/engineering/dittus-boelter-correlation?
            NU_coolant = .023 .* Re_coolant.^.8 .* Pr_coolant.^.4;
            hl_laminar_mask = Re_coolant <= 2300;
            hl_turbulent_mask = ~hl_laminar_mask;
            hl = zeros(size(mdot_channel));
            % only valid for laminar
            hl(hl_laminar_mask) = .023 .* cp_coolant_avg .* mdot_channel(hl_laminar_mask) ./ A_channel_range(hl_laminar_mask) .* Re_coolant(hl_laminar_mask).^-.2 .* (mu_coolant(hl_laminar_mask) .* cp_coolant_avg ./ k_coolant(hl_laminar_mask)).^(-2./3);
            % Dittus Boelter correlation, used with caution from Re = 2300 to Re = 10k, above 10k it's pretty good
            hl(hl_turbulent_mask) = NU_coolant(hl_turbulent_mask) .* k_coolant(hl_turbulent_mask) ./ dH_channels(hl_turbulent_mask);
            hg = .023 .* (throat_flow_density .* squeeze(flow_vel(i_pc, i_OF, i_eps, :, :, :, :, :, :))).^.8 ./ dt.^.2 .* Pr_flow.^.4 .* k_flow ./ mu_flow.^.8;

            % Recommended by Claude after I asked about needing Re, Pr,
            % Nusselt, and Biot numbers, and fric coeff for Brennen Kohlman
            % (for ANSYS) - PENDING REVIEW
            % Re_gas = throat_flow_density .* flow_vel .* dH_channels ./ mu_flow;
            % Re_coolant = dH_channels .* fluid_vel .* rho_coolant ./ mu_coolant;
            % % Coolant Prandtl (you have gas Pr from CEA)
            % Pr_coolant = mu_coolant .* cp_coolant_avg ./ k_coolant;
            % % Nusselt numbers (dimensionless heat transfer coefficient)
            % Nu_gas = hg .* dH_channels ./ k_flow;
            % Nu_coolant = hl .* dH_channels ./ k_coolant;
            % % Biot numbers (ratio of internal thermal resistance to surface resistance)
            % Bi_gas = hg .* wallts ./ k_walls;
            % Bi_coolant = hl .* wallts ./ k_coolant;
            % % Friction factors (Blasius for turbulent smooth pipes)
            % f_gas = 0.079 ./ Re_gas.^0.25;
            % f_coolant = 0.079 ./ Re_coolant.^0.25;

            % Get radial heat transfer through wall and hot wall temps
            % Accounts for hl coefficient being for heat transfer over
            %   larger area than hg coefficient (hot wall area facing
            %   channel per-cross-section is larger inside the channel than
            %   in the chamber)
            q(i_pc, i_OF, i_eps, :, :, :, :, :, :) = (squeeze(throat_flow_temp(i_pc, i_OF, i_eps, :, :, :, :, :, :)) - T_coolant_range) ./ (1 ./ hg + wallt_range ./ k_walls + 1 ./ hl ./ ((dt + wallt_range) ./ dt));
            Twg(i_pc, i_OF, i_eps, :, :, :, :, :, :) = squeeze(throat_flow_temp(i_pc, i_OF, i_eps, :, :, :, :, :, :)) - squeeze(q(i_pc, i_OF, i_eps, :, :, :, :, :, :)) ./ hg;
            Twl(i_pc, i_OF, i_eps, :, :, :, :, :, :) = T_coolant_range + squeeze(q(i_pc, i_OF, i_eps, :, :, :, :, :, :)) ./ hl;

            yieldstress_alloy(i_pc, i_OF, i_eps, :, :, :, :, :, :) = getYieldStress(squeeze(Twg(i_pc, i_OF, i_eps, :, :, :, :, :, :)));

            %% Calculate nozzle contour geometry

            % Doing Rao approximation bell contour nozzle
            % GT slides:
            % From: https://www.seitzman.gatech.edu/classes/ae6450/nozzle_geometries.pdf
            
            % Assemble bottom up, starting with nozzle exit, but not
            % necessarily calculated in that order:

            % Sections in order by calculation:
            % ----------------
            % Nozzle circular
            % Chamber circular widen
            % Nozzle parabolic
            % Chamber circular narrow
            % Chamber linear

            %% Nozzle circular
            % Give the contour arrays a few datapoints
            % Keep using numelsnoz
            % Makes an arc around 0,0 and shifts later for simplicity
            Rt = dt ./ 2;
            R1(i_pc, i_OF, i_eps, :, :, :, :, :, :) = .382 .* Rt;
            % GT slides say 15 degrees, but using 20 because for
            % 100% fpct and ~5 expansion ratio, 20's a good thetaN and it's
            % close enough to 15 degrees
            % Convert to radians
            thetaN(i_pc, i_OF, i_eps, :, :, :, :, :, :) = 20 ./ 180 .* pi;
            xN(i_pc, i_OF, i_eps, :, :, :, :, :, :) = squeeze(R1(i_pc, i_OF, i_eps, :, :, :, :, :, :)) .* sin(squeeze(thetaN(i_pc, i_OF, i_eps, :, :, :, :, :, :)));
            L_nozzle_circular = squeeze(xN(i_pc, i_OF, i_eps, :, :, :, :, :, :));
            
            % Basically makes a linspace from 0 to L_nozzle_circular's
            % value at every index in L_nozzle_circular
            z_nozzle_circular = fRshpSmmtnArrys(linspace(0, 1, numelsnoz), ndims(L_nozzle_circular));
            z_nozzle_circular = z_nozzle_circular .* L_nozzle_circular;

            rc = Rt + squeeze(R1(i_pc, i_OF, i_eps, :, :, :, :, :, :));

            r_nozzle_circular = rc - (squeeze(R1(i_pc, i_OF, i_eps, :, :, :, :, :, :)).^2 - z_nozzle_circular.^2).^.5;

            % Flip upside down (flip last dimension)
            r_nozzle_circular = flip(r_nozzle_circular, ndims(r_nozzle_circular));
            % Shift in z
            fpct = 1;
            L_nozzle = fpct .* Rt ./ tan(squeeze(thetaN(i_pc, i_OF, i_eps, :, :, :, :, :, :))) .* (sqrt(expansion_ratio) - 1 + 1.5 .* (1 ./ cos(squeeze(thetaN(i_pc, i_OF, i_eps, :, :, :, :, :, :))) - 1));
            z_nozzle_circular = z_nozzle_circular + L_nozzle - L_nozzle_circular;

            %% Chamber circular widen
            R1p(i_pc, i_OF, i_eps, :, :, :, :, :, :) = Rt .* 1.5;
            Rc(i_pc, i_OF, i_eps, :, :, :, :, :, :) = CR.^.5 .* Rt;
            % Side note: can calculate max hoop stress now
            % Equa: hoop_stress_doghouse = pd/2t
            % This hoop stress calc is calculating the "doghouse effect"
            %   stress
            hoop_stress_doghouse(i_pc, i_OF, i_eps, :, :, :, :, :, :) = pc .* 2 .* squeeze(Rc(i_pc, i_OF, i_eps, :, :, :, :, :, :)) ./ 2 ./ wallt_range;
            % This hoop stress calc is calculating the stress on the fins
            %   between the channels, starting from one hot-wall-thickness
            %   into the hot-wall and ending where the fin contacts the
            %   outer jacket wall
            hoop_stress_fins(i_pc, i_OF, i_eps, :, :, :, :, :, :) = pc .* 2 .* (squeeze(Rc(i_pc, i_OF, i_eps, :, :, :, :, :, :)) + wallt_range) ./ 2 ./ d_channel_range;
            % Arbitrary chamber converging section geometry choice. Made it
            % a circle for simplicity
            % Try deriving alpha from scale of how far in between throat
            % radius and chamber radius I want arc to stop
            % Using a while loop to avoid impossible geometric constraints
            etaRc = .2;
            while true
                alpha(i_pc, i_OF, i_eps, :, :, :, :, :, :) = acos(1 - (etaRc .* (squeeze(Rc(i_pc, i_OF, i_eps, :, :, :, :, :, :)) - Rt)) ./ squeeze(R1p(i_pc, i_OF, i_eps, :, :, :, :, :, :)));
    
                L_chamber_circular_widen = squeeze(R1p(i_pc, i_OF, i_eps, :, :, :, :, :, :)) .* sin(squeeze(alpha(i_pc, i_OF, i_eps, :, :, :, :, :, :)));
    
                z_chamber_circular_widen = fRshpSmmtnArrys(linspace(0, 1, numelsnoz), ndims(L_chamber_circular_widen));
                z_chamber_circular_widen = z_chamber_circular_widen .* L_chamber_circular_widen;
    
                rc = Rt + squeeze(R1p(i_pc, i_OF, i_eps, :, :, :, :, :, :));
    
                % No need to flip upside down
                r_chamber_circular_widen = rc - (squeeze(R1p(i_pc, i_OF, i_eps, :, :, :, :, :, :)).^2 - z_chamber_circular_widen.^2).^.5;
    
                % Shift in z
                z_chamber_circular_widen = z_chamber_circular_widen + z_nozzle_circular(:, :, :, :, :, :, end);
    
                %% Nozzle parabolic
                % Start with coord system the same way as done with the nozzle
                % circular section, then shift later
                % Based on the low expansion ratios we'll be using and high
                % fpct for low thetae (in context of Rao):
                % Use thetaN, the nozzle-circular end-part's angle.
                % thetaN stays the same
                Re(i_pc, i_OF, i_eps, :, :, :, :, :, :) = expansion_ratio.^.5 .* Rt; 
                L_nozzle_parabolic(i_pc, i_OF, i_eps, :, :, :, :, :, :) = L_nozzle - L_nozzle_circular;
    
                z_nozzle_parabolic = fRshpSmmtnArrys(linspace(0, 1, numelsnoz), ndims(squeeze(L_nozzle_parabolic(i_pc, i_OF, i_eps, :, :, :, :, :, :))));
                z_nozzle_parabolic = z_nozzle_parabolic .* squeeze(L_nozzle_parabolic(i_pc, i_OF, i_eps, :, :, :, :, :, :));
    
                % Nonlinear three variable system - solved by hand
                k = squeeze(Re(i_pc, i_OF, i_eps, :, :, :, :, :, :)).^2 ./ (2 .* tan(squeeze(thetaN(i_pc, i_OF, i_eps, :, :, :, :, :, :))) .* squeeze(L_nozzle_parabolic(i_pc, i_OF, i_eps, :, :, :, :, :, :)) + 2 .* squeeze(Re(i_pc, i_OF, i_eps, :, :, :, :, :, :)));
                b = 2 .* tan(squeeze(thetaN(i_pc, i_OF, i_eps, :, :, :, :, :, :))) .* k;
                h = -b ./ (4 .* tan(squeeze(thetaN(i_pc, i_OF, i_eps, :, :, :, :, :, :))).^2);            
    
                r_nozzle_parabolic = (b .* (z_nozzle_parabolic - h)).^.5 + k;
                
                % r not controlled for thetaN at z=0 constraint, so have to
                % account for that with a shift
                rShift = r_nozzle_circular(:, :, :, :, :, :, 1) - r_nozzle_parabolic(:, :, :, :, :, :, 1);
                r_nozzle_parabolic = r_nozzle_parabolic + rShift;
            
                % Flip z
                z_nozzle_parabolic = flip(z_nozzle_parabolic, ndims(z_nozzle_parabolic));
    
                %% Chamber circular narrow
                % Yes I realize it isn't circular, it used to be and I
                % never changed the name
                % Angle stays same as where it meets chamber circular widen
                % section.
                % alpha stays the same
    
                r1 = squeeze(Rc(i_pc, i_OF, i_eps, :, :, :, :, :, :)) - r_chamber_circular_widen(:, :, :, :, :, :, end);
    
                % Parabolic curve, solved by hand for a and b
                b = tan(squeeze(alpha(i_pc, i_OF, i_eps, :, :, :, :, :, :)));
                a = -b.^2 ./ (4 .* r1);
    
                z1 = (-b + (b.^2 + 4 .* a .* r1).^.5) ./ (2 .* a);
                if isreal(z1)
                   break
                else
                    etaRc = etaRc + .1 .* (1 - etaRc);
                end
            end
            L_chamber_circular_narrow(i_pc, i_OF, i_eps, :, :, :, :, :, :) = z1;

            z_chamber_circular_narrow = fRshpSmmtnArrys(linspace(0, 1, numelsnoz), ndims(squeeze(L_chamber_circular_narrow(i_pc, i_OF, i_eps, :, :, :, :, :, :))));
            z_chamber_circular_narrow = z_chamber_circular_narrow .* squeeze(L_chamber_circular_narrow(i_pc, i_OF, i_eps, :, :, :, :, :, :));

            r_chamber_circular_narrow = a .* z_chamber_circular_narrow.^2 + b .* z_chamber_circular_narrow;

            % r and z should already be in the proper orientation, so just
            % shift r out and z out
            z_chamber_circular_narrow = z_chamber_circular_narrow + z_chamber_circular_widen(:, :, :, :, :, :, end);
            r_chamber_circular_narrow = r_chamber_circular_narrow + r_chamber_circular_widen(:, :, :, :, :, :, end);

            %% Chamber linear

            % Get volume of chamber circular narrow and chamber circular
            % widen sections to calculate what remains of Vc, derive
            % chamber linear dimensions from that remainder of Vc
            % Account for widen section volume first
            dz_chamber_circular_widen = diff(z_chamber_circular_widen, [], ndims(z_chamber_circular_widen));
            Vcl = squeeze(Vc(i_pc, i_OF, i_eps, :, :, :, :, :, :)) - pi .* sum(r_chamber_circular_widen(:, :, :, :, :, :, 1:end-1).^2 .* dz_chamber_circular_widen, ndims(r_chamber_circular_widen));
            % Then for narrow section volume
            dz_chamber_circular_narrow = diff(z_chamber_circular_narrow, [], ndims(z_chamber_circular_narrow));
            Vcl = Vcl - pi .* sum(r_chamber_circular_narrow(:, :, :, :, :, :, 1:end-1).^2 .* dz_chamber_circular_narrow, ndims(r_chamber_circular_narrow));

            L_chamber_linear(i_pc, i_OF, i_eps, :, :, :, :, :, :) = Vcl ./ (pi .* squeeze(Rc(i_pc, i_OF, i_eps, :, :, :, :, :, :)).^2);

            % r_chamber_linear = Rc .* ones(1, 50);
            r_chamber_linear = fRshpSmmtnArrys(ones(1, numelsnoz), ndims(squeeze(L_chamber_circular_narrow(i_pc, i_OF, i_eps, :, :, :, :, :, :))));
            r_chamber_linear = r_chamber_linear .* squeeze(Rc(i_pc, i_OF, i_eps, :, :, :, :, :, :));

            z_chamber_linear = fRshpSmmtnArrys(linspace(0, 1, numelsnoz), ndims(squeeze(L_chamber_linear(i_pc, i_OF, i_eps, :, :, :, :, :, :))));
            z_chamber_linear = z_chamber_linear .* squeeze(L_chamber_linear(i_pc, i_OF, i_eps, :, :, :, :, :, :));

            % Shift z
            z_chamber_linear = z_chamber_linear + max(z_chamber_circular_narrow, [], 7);

            %% Compile all sections to export to app
            r_engine(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = cat(7, r_nozzle_parabolic, r_nozzle_circular, r_chamber_circular_widen, r_chamber_circular_narrow, r_chamber_linear);
            z_engine(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = cat(7, z_nozzle_parabolic, z_nozzle_circular, z_chamber_circular_widen, z_chamber_circular_narrow, z_chamber_linear);
            vol_engine(i_pc, i_OF, i_eps, :, :, :, :, :, :) = pi .* sum((squeeze(r_engine(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)) + 3 .* wallt_range).^2 - squeeze(r_engine(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)).^2 - num_channels_range .* d_channel_range.^2, 7);
            
            % Axial temperature
            area_ratios = squeeze(r_engine(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)).^2 ./ Rt.^2;
            % area_ratios is the same for every for loop iteration varying
            %   across expansion_ratio_range
            area_ratios = squeeze(area_ratios(1, 1, 1, 1, 1, 1, :));
            % Hardcoded 7 because ndims(area_ratios) might not return the
            %   right number if the final dimension is singleton
            idx_throat = 20; % Where the end of r_chamber_circular_narrow is in r_engine, because r_chamber_circular_narrow actually progresses correctly- forwards
            chamber_flow_mach = cea_out.get_Chamber_MachNumber(Pc = pc * 0.000145038, MR = OF, fac_CR = CR);
            chamber_flow_props = double(cea_out.get_Chamber_MolWt_gamma(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio));
            chamber_flow_gamma = chamber_flow_props(2);
            % Rearranged stagnation temp equation with knowns at chamber
            stag_flow_temp = chamber_flow_temp ./ (1 + (chamber_flow_gamma - 1) ./ 2 .* chamber_flow_mach.^2).^-1;
            % Can get some axial temps from get_Temperatures by CEA
            %   (post-throat/supersonic)
            % size_super = )
            % area_ratios is 1D array
            numels_area_ratio = numel(area_ratios);
            size_array_temp = idx_throat;
            axial_temps_super = zeros(1, size_array_temp);
            for i_station = 1:(idx_throat - 1)
                area_ratio = area_ratios(i_station);
                placeholder_temp_array = double(cea_out.get_Temperatures(Pc = pc * 0.000145038, MR = OF, eps = area_ratio));
                axial_temps_super(i_station) = placeholder_temp_array(3);
            end
            % throat_flow_temp is basically a scalar past the first three
            %   iterators
            axial_temps_super(idx_throat) = squeeze(throat_flow_temp(i_pc, i_OF, i_eps, 1));
            % Get pre-throat/subsonic station temperatures
            size_array_temp = numels_area_ratio - idx_throat;
            axial_temps_sub = zeros(1, size_array_temp);
            for i_station = 1:numels_area_ratio
                % Use get_full_cea_output, parse large output string
                % axial_temps_sub = cea_out.get_full_cea_output()
            end
            axial_temps = cat(7, axial_temps_super, axial_temps_sub);

            % For now, this is just approximated using the throat heat
            %   transfer
            % delta T of coolant just at throat cross section
            deltaT_coolantslice = squeeze(q(i_pc, i_OF, i_eps, :, :, :, :, :, :)) .* d_channel_range ./ mdot_channel ./ cp_coolant_avg;
            % assume at worst, it's that much at EVERY cross section
            % per-channel delta T
            l_channel = 2 .* max(squeeze(z_engine(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)), [], 7);
            deltaT_coolant = deltaT_coolantslice .* l_channel;
            T_coolant_f(i_pc, i_OF, i_eps, :, :, :, :, :, :) = T_coolant_range + deltaT_coolant;

            % Very crude estimate, basically assumes solid engine. no
            %   cavities inside aka no chamber
            % Honestly meaningless atm, was used earlier but now not so
            %   much. May implement more in the future
            TWR(i_pc, i_OF, i_eps, :, :, :, :, :, :) = squeeze(thrust(i_pc, i_OF, i_eps, :, :, :, :, :, :)) ./ (g0 .* densalloy .* squeeze(vol_engine(i_pc, i_OF, i_eps, :, :, :, :, :, :)));
            if ~isreal(r_engine)
                disp('non real physical engine design');
            end
        end
    end
end
% min coolant press to prevent it boiling in the channels (w/in DF)
P_coolant_min = fpvapeth(T_coolant_f) ./ DFcoolantpress;
thermstress(:, :, :, :, :, :, :, :, :) = E_alloy .* alpha_alloy .* (Twg - Twl);

% Get runtime diagnostics
runtime = toc(runtime)
num_engines = prod(ndsvar)
timeperengine = runtime / num_engines

% Haven't finished doing all the unit conversions yet for imperial system
% outputs
% constants
% dimensionless to %
eta_cstar_min = eta_cstar_min .* 100;
% K to deg C
T_AlSi10Mg_melt = T_AlSi10Mg_melt - 273.15;
T_w_max = T_w_max - 273.15;
% Already in MPa
yieldstress_max = yieldstress_alloy ./ DFstress;
% independent vars
% update this compact cell array too
% ranges = {pc_range, OF_range, expansion_ratio_range, mdot_range, d_channel_range, num_channels_range, T_coolant_range, wallt_range, k_wall_range};
% pc_range -> Pa to psi
ranges{1} = ranges{1} .* .000145038;
% mdot_range -> kgm/s to lbm/s
ranges{4} = ranges{4} .* 2.20462;
% d_channel_range -> m to mm
ranges{5} = ranges{5} .* 10.^3;
% T_coolant_range -> K to deg C
ranges{7} = ranges{7} - 273.15;
% wallt_range -> m to mm
ranges{8} = ranges{8} .* 10.^3;
% dependent vars
% K to deg C
T_coolant_f = T_coolant_f - 273.15;
Twg = Twg - 273.15;
Twl = Twl - 273.15;
% m to in
dt = sqrt(At ./ pi .* 4) .* 39.3701;
de = sqrt(Ae ./ pi .* 4) .* 39.3701;
% m^3 to mm^3
Vc = Vc .* 1000.^3;
% m/s to ft/s
cstar = cstar .* 3.28084;
cstar_theo = cstar_theo .* 3.28084;
% N to lbf
thrust = thrust .* 0.224809;
% Pa to psi
P_coolant_min = P_coolant_min .* .000145038;
% dimensionless to %
eta_cstar = eta_cstar .* 100;
% Pa to MPa
thermstress = thermstress ./ 10.^6;
hoop_stress_doghouse = hoop_stress_doghouse ./ 10.^6;
hoop_stress_fins = hoop_stress_fins ./ 10.^6;
total_stress_doghouse = thermstress + hoop_stress_doghouse;
total_stress_fins = thermstress + hoop_stress_fins;
% Engine volume and contour info
% m to mm
r_engine = r_engine .* 10.^3;
z_engine = z_engine .* 10.^3;
% m^3 to mm^3
vol_engine = vol_engine .* 10.^9;

contourvalnames = ["vol_engine", "L_nozzle_parabolic", "Re", "thetaN", "xN", "R1", "R1p", "alpha", "L_chamber_circular_narrow", "Rc", "L_chamber_linear", "dt"];

% Export to mat file for use in app
save(filename_datastorage, "T_AlSi10Mg_melt", "T_w_max", "yieldstress_alloy", "yieldstress_max", "eta_cstar_min", "list_var_names", "ranges", "pc_range", "OF_range", "expansion_ratio_range", "mdot_range", "d_channel_range", "num_channels_range", "T_coolant_range", "wallt_range", "k_walls", "thermstress", "total_stress_doghouse", "total_stress_fins", "hoop_stress_doghouse", "hoop_stress_fins", "T_coolant_f", "P_coolant_min", "Twg", "Twl", "q", "Isp", "prop_cost", "dt", "cstar", "cstar_theo", "thrust", "de", "eta_cstar", "Vc", "CR", "flow_sonic", "r_engine", "z_engine", "vol_engine", "L_nozzle_parabolic", "Re", "thetaN", "xN", "R1", "R1p", "alpha", "L_chamber_circular_narrow", "Rc", "L_chamber_linear", "contourvalnames", "TWR");


%% More complex functions
% For vaporization pressure calclulation to check coolant doesn't boil in
% coolant channels
% antoine equation (T in K) - source: NIST, pvapeth in Pa
function pvapeth = fpvapeth(T)
    Tcrit_eth = 513.91;

    idxslw = T < 273.15;
    idxs1 = 273.15 <= T & T < 292.77;
    idxs2 = 292.77 <= T & T < 364.8;
    idxs3 = 364.8 <= T & T < Tcrit_eth;
    idxshgh = Tcrit_eth <= T;
    idxsvalid = idxs1 | idxs2 | idxs3;

    pvapeth = nan(size(T));
    A = pvapeth;
    B = pvapeth;
    C = pvapeth;

    A(idxs1) = 5.37229;
    B(idxs1) = 1670.409;
    C(idxs1) = -40.191;

    A(idxs2) = 5.24677;
    B(idxs2) = 1598.673;
    C(idxs2) = -46.424;

    A(idxs3) = 4.92531;
    B(idxs3) = 1432.526;
    C(idxs3) = -61.819;

    if any(idxslw)
        warning("some input T's too low in fpvapeth function");
    end
    if any(idxshgh)
        warning("some input T's above T crit eth");
    end

    % bar to Pa
    pvapeth(idxsvalid) = 10.^(A(idxsvalid) - B(idxsvalid) ./ (T(idxsvalid) + C(idxsvalid))) .* 10.^5;
end

% Assuming material is AlSi10Mg
% Vectorized
% temp is Twg, hot wall temperature, in Kelvin
% Output is in MPa
function yieldStress = getYieldStress(temp)
    % Data sheet on material properties for AlSi10Mg - https://fathommfg.com/wp-content/uploads/2020/11/EOS_Aluminium_AlSi10Mg_en.pdf
    % MPa to Pa (270 +- 10) - heat treated would be 245, 230 +- 15
    % Also got these tabulated value from "The variation of the yield stress vs
    %   temperature for experimented materials"
    % K to deg C
    temp = temp - 273.15;
    % deg C, MPa
    % Also adding really high value just to avoid having to fix the NaN
    %   error returned by interp1 when you feed it a value outside of the x
    %   range
    temps = [0;
             150;
             250;
             350;
             450;
             10000];
    stresses = [280;
                280;
                210;
                135;
                65;
                65];
    yieldStress = interp1(temps, stresses, temp);
end