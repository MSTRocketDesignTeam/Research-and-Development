%% Engine Performance Trades
% Units start in SI, converted to Imperial on front end
% -------------------------------------------------------------------------

clc 
clear
% Needed to get NASA CEA to work
py.importlib.import_module("rocketcea");
% For data export later to be graphed in MATLAB app
filename_datastorage = append(pwd, "\engines.mat");

%% Target Performance Parameters
% Derived from requirements
% c* efficiency
% total cost of consumed propellants
% burntime

eta_cstar_min = .6;
cost_max = 200; % USD
burntime = 10; % s

%% Constants
% g0                        accel. due to gravity
% R_universal               universal gas constant
% pa                        ambient pressure
% MoS                       (for some things)
% T_AlSi10Mg_melt           chamber alloy melting temp
% cost_nit                  $ per kgm
% cost_eth                  $ per kgm
% yieldstress_alloy
% E_alloy
% alpha_alloy               % CTE
% lchannel                  worst-case scenario length of a single coolant
                            %   channel

g0 = 9.81;
R_universal = 8.314;

% Calculate pa at Rolla altitude - Starting these calcs in Imperial
p_sea = 14.7 .* 6894.76; % psi to pa
T_sea = 288.15; % K
mlr_wgt_air = .02896968; % kg/mol
L = .00976;
fPAtAlt = @(alt) p_sea .* (1 + (L .* alt) ./ T_sea).^(-g0 .* mlr_wgt_air ./ (R_universal .* L));
alt = 342;
pa = fPAtAlt(alt);

% RDT technically has no MoS for these
MoSstress = 2;
MoStemp = 1.15;
MoScoolantpress = 1.2;
% deg C to K
T_AlSi10Mg_melt = 570 + 273.15; % K
% $/lbm to $/kgm
cost_nit = 6 .* 2.20462;
cost_eth = .46 .* 2.20462;
% Data sheet on material properties for AlSi10Mg - https://fathommfg.com/wp-content/uploads/2020/11/EOS_Aluminium_AlSi10Mg_en.pdf
yieldstress_alloy = 280E6; % MPa > Pa, 270 +- 10. heat treated would be 245, 230 +- 15
E_alloy = 85E6; % MPa > Pa, 75 +- 10. heat treated would be 80, 70 +- 10
% From data sheet: https://www.eos.info/var/assets/03_system-related-assets/material-related-contents/metal-materials-and-examples/metal-material-datasheet/aluminium/material_datasheet_eos_aluminium-alsi10mg_en_web.pdf
% thermal expansion coeff
alpha_alloy = 27E-6; % K^-1
% in to m
l_channel = 5 .* .0254;

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

% psi to Pa
pc_nom = 350 .* 6894.76;
expansion_ratio_nom = 4;
% lbm/s -> kgm/s
mdot_nom = 3 .* .45359237;
% mm to m
d_channel_nom = 1 .* 1E-3;
% convert deg C to K LATER
T_coolant_nom = 30;
% mm to m
wallt_nom = 1 ./ 1000;
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
% Applying MoS'es to hot wall temp and hot wall thermal stress
T_w_max = T_AlSi10Mg_melt ./ MoStemp;
yieldstress_max = yieldstress_alloy ./ MoSstress;

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

%% Make iteration ranges for independent vars
stdvarrange = [.8, .9, 1, 1.1];
highdefvarrange = linspace(.7, 1.3, 30);
pc_range = pc_nom .* stdvarrange;
OF_range = linspace(2, 7, 4);
expansion_ratio_range = expansion_ratio_nom .* stdvarrange;
mdot_range = mdot_nom .* stdvarrange;
d_channel_range = d_channel_nom .* stdvarrange;
num_channels_nom = 100;
num_channels_range = num_channels_nom .* stdvarrange;
T_coolant_range = T_coolant_nom .* stdvarrange;
% deg C to K
T_coolant_range = T_coolant_range + 273.15;
wallt_range = wallt_nom .* stdvarrange;
% W/m-K - Largest range I found in the paper: On the thermal conductivity of AlSi10Mg and lattice structures made by laser powder bed fusion Richard R.J. Sélo, Sam Catchpole-Smith, Ian Maskery, Ian Ashcroft, Christopher Tuck
% k_wall range is from ~98 to ~183
k_wall_nom = (98 + 183) / 2;
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

% Start runtime tracking for engine iterations
runtime = tic;
% Condensed engine iteration structure brancher
% CEA syntax:
% py.rocketcea.cea_obj.CEA_Obj(fuelName, oxName) - see input_cards.py
% CEA_Obj class functions: (Pc, MR, eps); Pc in psi, MR is OF, eps is expansion_ratio
cea_out = py.rocketcea.cea_obj.CEA_Obj(fuelName = fuel_name, oxName = ox_name);
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
            gamma_es = double(cea_out.get_Throat_MolWt_gamma(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio));
            mlr_wgt_flow_t = gamma_es(1) .* 10.^-3; % lbm/lbm-mol -> kg/mol
            gamma_e = gamma_es(2);
            % Needed for thrust calcs later
            pex = pc ./ cea_out.get_PcOvPe(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio);
            Cf = sqrt(2 .* gamma_e.^2 ./ (gamma_e - 1) .* (2 ./ (gamma_e + 1)).^((gamma_e + 1) ./ (gamma_e - 1)).*(1 - (pex ./ pc).^((gamma_e - 1) ./ gamma_e))) + expansion_ratio .* ((pex - pa) ./ pc);
            % Thought theoretical Isp might be interesting to compare to
            % predicted-actual, not really used atm though (ca 11/04/2025)
            Isp_theo = cea_out.get_Isp(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio);
            % Gets throat temp needed for several calculations later
            % deg R to K
            Tts = double(cea_out.get_Temperatures(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio)) .* 5 ./ 9;
            Tt = Tts(2);
            % Fix units from CEA output lbm/ft^3 -> kg/m^3
            rho_ts = double(cea_out.get_Densities(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio)) .* 16.018463;
            rho_t = rho_ts(2);
            % Gets things needed for gas film coeff calc later
            throat_trans = double(cea_out.get_Throat_Transport(Pc = pc * 0.000145038, MR = OF, eps = expansion_ratio));
            mu_flow = throat_trans(2) .* 10.^-3;
            k_flow = throat_trans(3) .* .4184;
            Pr_flow = throat_trans(4);
            % ft/s -> m/s
            cstar_theo(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = cea_out.get_Cstar(Pc = pc * 0.000145038, MR = OF) .* .3048;

            % Get some KPPs and other sizing info for trade alternatives
            flow_sonic = sqrt(gamma_e * R_universal / mlr_wgt_flow_t * Tt);
            At(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = mdot_range ./ flow_sonic ./ rho_t;
            Ae(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = squeeze(At(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)) .* expansion_ratio;
            thrust(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = Cf .* squeeze(At(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)) .* pc;
            prop_cost(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = burntime .* mdot_range .* (cost_nit .* (1 ./ (1 + OF) .* OF) + cost_eth .* (1 ./ (1 + OF)));
            Isp(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = squeeze(thrust(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)) ./ mdot_range ./ g0;
            cstar(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = squeeze(Isp(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)) .* g0 ./ Cf;
            eta_cstar(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = squeeze(cstar(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)) ./ squeeze(cstar_theo(i_pc, i_OF, i_eps, :, :, :, :, :, :, :));

            % -------------------------------------------------------------
            % Needed to find film coeffs
            % flow velocity in chamber (just throat flow vel or sonic flow
            % vel bc we're at throat)
            flow_vel = flow_sonic;
            % hydraulic radius of channel
            dH_channels = 4 .* A_channel_range ./ (3 .* d_channel_range + 2 .* d_channel_range ./ sqrt(2));
            % Remember it's just fuel in channels
            mdot_channel = mdot_range .* (1 / (1 + OF)) ./ num_channels_range;
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
            % hydraulic radius is perimeter / circumference here
            Dht = sqrt(squeeze(At(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)) .* 4 ./ pi) .* pi;
            Ret = rho_t .* flow_vel .* squeeze(Dht(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)) ./ mu_flow;
            hl = .023 .* cp_coolant_avg .* mdot_channel ./ A_channel_range .* (dH_channels .* fluid_vel .* rho_coolant ./ mu_coolant).^-.2 .* (mu_coolant .* cp_coolant_avg ./ k_coolant).^(-2./3);
            hg = .023 .* (rho_t .* flow_vel).^.8 ./ dH_channels.^.2 .* Pr_flow.^.4 .* k_flow ./ mu_flow.^.8;
         
            % Recommended by Claude after I asked about needing Re, Pr,
            % Nusselt, and Biot numbers, and fric coeff for Brennen Kohlman
            % (for ANSYS) - PENDING REVIEW
            % Re_gas = rho_t .* flow_vel .* dH_channels ./ mu_flow;
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
            q(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = (Tt - T_coolant_range) ./ (1 ./ hg + wallt_range ./ k_walls + 1 ./ hl);
            Twg(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = Tt - squeeze(q(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)) ./ hg;
            Twl(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = T_coolant_range + squeeze(q(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)) ./ hl;

            % delta T of coolant just at throat cross section
            deltaT_coolantslice = squeeze(q(i_pc, i_OF, i_eps, :, :, :, :, :, :, :)) .* d_channel_range ./ mdot_channel ./ cp_coolant_avg;
            % assume at worst, it's that much at EVERY cross section
            % per-channel delta T
            deltaT_coolant = deltaT_coolantslice .* l_channel;
            T_coolant_f(i_pc, i_OF, i_eps, :, :, :, :, :, :, :) = T_coolant_range + deltaT_coolant;
        end
    end
end
% min coolant press to prevent it boiling in the channels (w/in MoS)
P_coolant_min = fpvapeth(T_coolant_f) ./ MoScoolantpress;
thermstress(:, :, :, :, :, :, :, :, :, :) = E_alloy .* alpha_alloy .* (Twg - Twl);

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
% Pa to psi
yieldstress_alloy = yieldstress_alloy .* .000145038;
yieldstress_max = yieldstress_max .* .000145038;
E_alloy = E_alloy .* .000145038;
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
% m/s to ft/s
cstar = cstar .* 3.28084;
cstar_theo = cstar_theo .* 3.28084;
% N to lbf
thrust = thrust .* 0.224809;
% Pa to psi
P_coolant_min = P_coolant_min .* .000145038;
thermstress = thermstress .* .000145038;
% dimensionless to %
eta_cstar = eta_cstar .* 100;

% Export to mat file for use in app
save(filename_datastorage, "T_AlSi10Mg_melt", "T_w_max", "yieldstress_alloy", "yieldstress_max", "eta_cstar_min", "list_var_names", "ranges", "pc_range", "OF_range", "expansion_ratio_range", "mdot_range", "d_channel_range", "num_channels_range", "T_coolant_range", "wallt_range", "k_walls", "thermstress", "T_coolant_f", "P_coolant_min", "Twg", "q", "Isp", "prop_cost", "dt", "cstar", "cstar_theo", "thrust", "de", "eta_cstar");


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


%% Next
%   Claude, AI agent, gave some code to calc Re, Pr, Biot #, Nusselt #, and
%       friction coeff if wanted, but will leave commented out for now