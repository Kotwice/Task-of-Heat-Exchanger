
% PARAMETERS LIST

parameters.d = 10;
parameters.l_in = 10 * parameters.d;
parameters.l_out = 10 * parameters.d;

parameters.h_r = 0.5;
parameters.s_r = 5;
parameters.w_r = 1;
parameters.n = 10;

parameters.Re = 1e4;

parameters.p_0 = 0;
parameters.p_out = 0;

parameters.T_w = 1e3;
parameters.T_0 = 298.15;
parameters.T_in = 300;

% EXECUTION CODE

model = create_model(parameters);
model = build_component_smooth(model);
model = build_component_rough(model);
% model = set_studies(model);

mphsave(model, strcat(pwd, '\build_project_rectangle_3d.mph'));

% INITIALIZATION BLOCK

function model = create_model(parameters)

    import com.comsol.model.*
    import com.comsol.model.util.*

    model = ModelUtil.create('Model');

        model.param.set('d', strcat(num2str(parameters.d), '[mm]'), 'Diameter of Tube');
        model.param.set('h_r', strcat(num2str(parameters.h_r), '[mm]'), 'High of Roughness');
        model.param.set('s_r', strcat(num2str(parameters.s_r), '[mm]'), 'Step Length of Roughness');
        model.param.set('w_r', strcat(num2str(parameters.w_r), '[mm]'), 'Width of Roughness');

        model.param.set('rho_air', '1.1839 [kg/m^3]', 'Density of Air at 25 Cdeg');
        model.param.set('mu_air', '18.35 [uPa*s]', 'Dynamic Viscosity of Air at 25 Cdeg');
        model.param.set('Re', strcat(num2str(parameters.Re)), 'Reynolds Number');
        model.param.set('v_in', 'mu_air * Re / (d * rho_air)', 'Inlet Velocity');

        model.param.set('p_0', strcat(num2str(parameters.p_0), '[Pa]'), 'Initial Pressure');
        model.param.set('p_out', strcat(num2str(parameters.p_out), '[Pa]'), 'Output Pressure');
        model.param.set('T_w', strcat(num2str(parameters.T_w), '[K]'), 'Temperature of Wall');
        model.param.set('T_0', strcat(num2str(parameters.T_0), '[K]'), 'Initial Temperature');
        model.param.set('T_in', strcat(num2str(parameters.T_in ), '[K]'), 'Inlet Temperature');

        model.param.set('l_in', strcat(num2str(parameters.l_in), ' [mm]'), 'Length of Input Stabilization Domain');
        model.param.set('n', strcat(num2str(parameters.n)), 'Number of Ribs'); 
        model.param.set('l_out', strcat(num2str(parameters.l_out), ' [mm]'), 'Length of Output Stabilization Domain');

        model.param.set('h_r_min', '0.05 * d', 'Minimal Value of High of Roughness');
        model.param.set('h_r_step', '1 [mm]', 'Step Value of High of Roughness');
        model.param.set('h_r_max', '0.3 * d', 'Maximal Value of High of Roughness');

        model.param.set('s_r_min', '0.25 * d', 'Minimal Value of Step Length of Roughness');
        model.param.set('s_r_step', '1 [mm]', 'Step Value of Step Length of Roughness');
        model.param.set('s_r_max', 'd', 'Maximal Step Length of High of Roughness');   

        model.param.set('Re_min', '1e4', 'Minimal Reynolds Number');
        model.param.set('Re_step', '5 * 1e3', 'Step Reynolds Number');
        model.param.set('Re_max', '1e5', 'Maximal Reynolds Number');  

        model.param.set('F_ref_s', 'pi*d*n*(w_r+s_r)', 'Square of Exchange Surface at Smoothness Wall');
        model.param.set('P_ref_s', '2*(pi*d + n*(w_r+s_r))', 'Perimeter of Exchange Surface at Smoothness Wall');
        model.param.set('l_ref_s', '4*F_ref_s/P_ref_s', 'Reference Length at Smoothness Wall');

        model.param.set('F_ref_r', 'n*pi*(2*(d-2*h_r)*w_r+((d/2)^2-(d/2-h_r)^2))-pi*(d-2*h_r)*w_r', 'Square of Exchange Surface at Roughness Wall');
        model.param.set('P_ref_r', '2*(pi*d + n*(w_r+s_r))', 'Perimeter of Exchange Surface at Roughness Wall');
        model.param.set('l_ref_r', '4*F_ref_r/P_ref_r', 'Reference Length at Roughness Wall');

end

function model = build_component_smooth (model)

% GEOMETRY

    model.component.create('comp_s', true);

    model.component('comp_s').label('Smooth Wall');

    model.component('comp_s').geom.create('geom_s', 3);

    model.component('comp_s').geom('geom_s').lengthUnit('mm');

    model.component('comp_s').geom('geom_s').create('cyl1', 'Cylinder');
        model.component('comp_s').geom('geom_s').feature('cyl1').set('r', 'd / 2');
        model.component('comp_s').geom('geom_s').feature('cyl1').set('h', 'l_in');

    model.component('comp_s').geom('geom_s').create('cyl2', 'Cylinder');
        model.component('comp_s').geom('geom_s').feature('cyl2').set('r', 'd / 2');
        model.component('comp_s').geom('geom_s').feature('cyl2').set('h', 'n * (w_r + s_r)');
        model.component('comp_s').geom('geom_s').feature('cyl2').set('pos', {'0' '0' 'l_in'});

    model.component('comp_s').geom('geom_s').create('cyl3', 'Cylinder');
        model.component('comp_s').geom('geom_s').feature('cyl3').set('r', 'd / 2');
        model.component('comp_s').geom('geom_s').feature('cyl3').set('h', 'l_out');
        model.component('comp_s').geom('geom_s').feature('cyl3').set('pos', {'0' '0' 'l_in + n * (w_r + s_r)'});

    model.component('comp_s').geom('geom_s').run('fin');

% SELECTIONS

    model.component('comp_s').selection.create('sel_task_s', 'Cylinder');
        model.component('comp_s').selection('sel_task_s').label('Task');
        model.component('comp_s').selection('sel_task_s').set('r', 'd / 2');

    model.component('comp_s').selection.create('sel_fluid_d_s', 'Cylinder');
        model.component('comp_s').selection('sel_fluid_d_s').label('Fluid');

    model.component('comp_s').selection.create('sel_fluid_b_s', 'Adjacent');
        model.component('comp_s').selection('sel_fluid_b_s').label('Boundaries of Fluid');
        model.component('comp_s').selection('sel_fluid_b_s').set('input', {'sel_fluid_d_s'});

    model.component('comp_s').selection.create('sel_in_d_s', 'Cylinder');
        model.component('comp_s').selection('sel_in_d_s').set('top', '0');
        model.component('comp_s').selection('sel_in_d_s').set('bottom', '0');
        model.component('comp_s').selection('sel_in_d_s').label('Domain of Input');

    model.component('comp_s').selection.create('sel_out_d_s', 'Cylinder');
        model.component('comp_s').selection('sel_out_d_s').set('top', 'l_in + n * (w_r + h_r) + l_out');
        model.component('comp_s').selection('sel_out_d_s').set('bottom', 'l_in + n * (w_r + h_r) + l_out');
        model.component('comp_s').selection('sel_out_d_s').label('Domain of Output');

    model.component('comp_s').selection.create('sel_c_d_s', 'Difference');
        model.component('comp_s').selection('sel_c_d_s').set('add', {'sel_fluid_d_s'});
        model.component('comp_s').selection('sel_c_d_s').set('subtract', {'sel_in_d_s', 'sel_out_d_s'});
        model.component('comp_s').selection('sel_c_d_s').label('Domain of Center');

    model.component('comp_s').selection.create('sel_in_b_s', 'Adjacent');
        model.component('comp_s').selection('sel_in_b_s').label('Boundaries of Input');
        model.component('comp_s').selection('sel_in_b_s').set('input', {'sel_in_d_s'});

    model.component('comp_s').selection.create('sel_c_b_s', 'Adjacent');
        model.component('comp_s').selection('sel_c_b_s').label('Boundaries of Center');
        model.component('comp_s').selection('sel_c_b_s').set('input', {'sel_c_d_s'});

    model.component('comp_s').selection.create('sel_out_b_s', 'Adjacent');
        model.component('comp_s').selection('sel_out_b_s').label('Boundaries of Output');
        model.component('comp_s').selection('sel_out_b_s').set('input', {'sel_out_d_s'});

    model.component('comp_s').selection.create('sel_inlet_s', 'Cylinder');
        model.component('comp_s').selection('sel_inlet_s').set('entitydim', 2);
        model.component('comp_s').selection('sel_inlet_s').set('top', 0);
        model.component('comp_s').selection('sel_inlet_s').set('bottom', 0);
        model.component('comp_s').selection('sel_inlet_s').label('Inlet');

    model.component('comp_s').selection.create('sel_outlet_s', 'Cylinder');
        model.component('comp_s').selection('sel_outlet_s').set('entitydim', 2);
        model.component('comp_s').selection('sel_outlet_s').set('top', 'l_in + n * (w_r + s_r) + l_out');
        model.component('comp_s').selection('sel_outlet_s').set('bottom', 'l_in + n * (w_r + s_r) + l_out');
        model.component('comp_s').selection('sel_outlet_s').label('Outlet');

    model.component('comp_s').selection.create('sel_inlet_c_s', 'Cylinder');
        model.component('comp_s').selection('sel_inlet_c_s').set('entitydim', 2);
        model.component('comp_s').selection('sel_inlet_c_s').set('top', 'l_in');
        model.component('comp_s').selection('sel_inlet_c_s').set('bottom', 'l_in');
        model.component('comp_s').selection('sel_inlet_c_s').label('Inlet of Center');

    model.component('comp_s').selection.create('sel_outlet_c_s', 'Cylinder');
        model.component('comp_s').selection('sel_outlet_c_s').set('entitydim', 2);
        model.component('comp_s').selection('sel_outlet_c_s').set('top', 'l_in + n * (w_r + s_r)');
        model.component('comp_s').selection('sel_outlet_c_s').set('bottom', 'l_in + n * (w_r + s_r)');      
        model.component('comp_s').selection('sel_outlet_c_s').label('Outlet of Center');
    
    model.component('comp_s').selection.create('sel_c_wall_s', 'Difference');
        model.component('comp_s').selection('sel_c_wall_s').label('Wall of Center');
        model.component('comp_s').selection('sel_c_wall_s').set('entitydim', 2);
        model.component('comp_s').selection('sel_c_wall_s').set('add', {'sel_c_b_s'});
        model.component('comp_s').selection('sel_c_wall_s').set('subtract', {'sel_inlet_c_s' 'sel_outlet_c_s'});

    model.component('comp_s').selection.create('sel_in_wall_s', 'Difference');
        model.component('comp_s').selection('sel_in_wall_s').set('entitydim', 2);
        model.component('comp_s').selection('sel_in_wall_s').set('add', {'sel_in_b_s'});
        model.component('comp_s').selection('sel_in_wall_s').set('subtract', {'sel_inlet_s' 'sel_inlet_c_s'});
        model.component('comp_s').selection('sel_in_wall_s').label('Wall of Input');

    model.component('comp_s').selection.create('sel_out_wall_s', 'Difference');
        model.component('comp_s').selection('sel_out_wall_s').set('entitydim', 2);
        model.component('comp_s').selection('sel_out_wall_s').set('add', {'sel_out_b_s'});
        model.component('comp_s').selection('sel_out_wall_s').set('subtract', {'sel_outlet_s' 'sel_outlet_c_s'});
        model.component('comp_s').selection('sel_out_wall_s').label('Wall of Output');

    model.component('comp_s').selection.create('sel_wall_s', 'Difference');
        model.component('comp_s').selection('sel_wall_s').set('entitydim', 2);
        model.component('comp_s').selection('sel_wall_s').set('add', {'sel_fluid_b_s'});
        model.component('comp_s').selection('sel_wall_s').set('subtract', {'sel_inlet_s' 'sel_outlet_s'});
        model.component('comp_s').selection('sel_wall_s').label('Wall');

    model.component('comp_s').selection.create('sel_bc_s', 'Union');
        model.component('comp_s').selection('sel_bc_s').set('entitydim', 2);
        model.component('comp_s').selection('sel_bc_s').set('input', {'sel_c_wall_s' 'sel_in_wall_s' 'sel_out_wall_s'});
        model.component('comp_s').selection('sel_bc_s').label('Boundaries of Condition');
        
% PHYSICS

    model.component('comp_s').physics.create('spf_s', 'LESRBVM', 'geom_s');

        model.component('comp_s').physics('spf_s').selection.named('sel_fluid_d_s');

        model.component('comp_s').physics('spf_s').feature('init1').set('u_init', {'0' '0' 'v_in'});
        model.component('comp_s').physics('spf_s').feature('init1').set('p_init', 'p_0');

        model.component('comp_s').physics('spf_s').create('inl1', 'InletBoundary', 2);
        model.component('comp_s').physics('spf_s').feature('inl1').set('U0in', 'v_in');
            model.component('comp_s').physics('spf_s').feature('inl1').selection.named('sel_inlet_s');

        model.component('comp_s').physics('spf_s').create('out1', 'OutletBoundary', 2);
            model.component('comp_s').physics('spf_s').feature('out1').set('p0', 'p_out');
            model.component('comp_s').physics('spf_s').feature('out1').selection.named('sel_outlet_s');

        model.component('comp_s').physics('spf_s').tag('spf_s');

    model.component('comp_s').physics.create('ht_s', 'HeatTransferInSolidsAndFluids', 'geom_s');

        model.component('comp_s').physics('ht_s').selection.named('sel_task_s');

        model.component('comp_s').physics('ht_s').feature('init1').set('Tinit', 'T_0');

        model.component('comp_s').physics('ht_s').feature('fluid1').selection.named('sel_fluid_d_s');

        model.component('comp_s').physics('ht_s').create('ifl1', 'Inflow', 2);
            model.component('comp_s').physics('ht_s').feature('ifl1').set('Tustr', 'T_in');
            model.component('comp_s').physics('ht_s').feature('ifl1').selection.named('sel_inlet_s');

        model.component('comp_s').physics('ht_s').create('ofl1', 'ConvectiveOutflow', 2);
            model.component('comp_s').physics('ht_s').feature('ofl1').selection.named('sel_outlet_s');

        model.component('comp_s').physics('ht_s').create('temp1', 'TemperatureBoundary', 2);
            model.component('comp_s').physics('ht_s').feature('temp1').set('T0', 'T_w');
            model.component('comp_s').physics('ht_s').feature('temp1').selection.named('sel_bc_s');

        model.component('comp_s').physics('ht_s').tag('ht_s');

    model.component('comp_s').multiphysics.create('nitf_s', 'NonIsothermalFlow', 3);

% PROBES

    model.component('comp_s').probe.create('dom_int_1_s', 'Domain');
        model.component('comp_s').probe('dom_int_1_s').selection.named('sel_c_d_s');
        model.component('comp_s').probe('dom_int_1_s').set('type', 'integral');
        model.component('comp_s').probe('dom_int_1_s').set('expr', 'spf_s.rho * spf_s.U');
        model.component('comp_s').probe('dom_int_1_s').set('probename', 'int_rho_v');
        model.component('comp_s').probe('dom_int_1_s').set('unit', 'kg * m / s');
        model.component('comp_s').probe('dom_int_1_s').label('Domain Probe of Density Mass Rate');

    model.component('comp_s').probe.create('dom_int_2_s', 'Domain');
        model.component('comp_s').probe('dom_int_2_s').selection.named('sel_c_d_s');
        model.component('comp_s').probe('dom_int_2_s').set('type', 'integral');
        model.component('comp_s').probe('dom_int_2_s').set('expr', 'spf_s.rho * spf_s.U * ht_s.T');
        model.component('comp_s').probe('dom_int_2_s').set('probename', 'int_rho_v_T');
        model.component('comp_s').probe('dom_int_2_s').set('unit', 'kg * m * K / s');
        model.component('comp_s').probe('dom_int_2_s').label('Domain Probe of Density rho*v*T');

    model.component('comp_s').probe.create('var_T_m_s', 'GlobalVariable');
        model.component('comp_s').probe('var_T_m_s').set('probename', 'T_r_m');
        model.component('comp_s').probe('var_T_m_s').set('expr', 'int_rho_v_T / int_rho_v');
        model.component('comp_s').probe('var_T_m_s').set('unit', 'K');
        model.component('comp_s').probe('var_T_m_s').set('descractive', true);
        model.component('comp_s').probe('var_T_m_s').set('descr', 'Mass Average Temperature of Smooth Wall');
        model.component('comp_s').probe('var_T_m_s').label('Domain Probe of Mass Average Temperature');

    model.component('comp_s').probe.create('bnd_alpha_s', 'Boundary');
        model.component('comp_s').probe('bnd_alpha_s').selection.named('sel_c_wall_s');
        model.component('comp_s').probe('bnd_alpha_s').set('probename', 'alpha_s');
        model.component('comp_s').probe('bnd_alpha_s').set('expr', 'abs(ht_s.ntflux/(T_r_m - T_w))');
        model.component('comp_s').probe('bnd_alpha_s').set('descractive', true);
        model.component('comp_s').probe('bnd_alpha_s').set('descr', 'Average Alpha of Smooth Wall');
        model.component('comp_s').probe('bnd_alpha_s').label('Boundary Probe of Average Alpha');

    model.component('comp_s').probe.create('bnd_lambda_s', 'Boundary');
        model.component('comp_s').probe('bnd_lambda_s').selection.named('sel_c_wall_s');
        model.component('comp_s').probe('bnd_lambda_s').set('expr', 'ht_s.k_iso');
        model.component('comp_s').probe('bnd_lambda_s').set('probename', 'lambda_r');
        model.component('comp_s').probe('bnd_lambda_s').set('descractive', true);
        model.component('comp_s').probe('bnd_lambda_s').set('descr', 'Avarage Thermal Conductivity of Smooth Wall');
        model.component('comp_s').probe('bnd_lambda_s').label('Boundary Probe of Avarege Lambda');

    model.component('comp_s').probe.create('var_Nu_s', 'GlobalVariable');
        model.component('comp_s').probe('var_Nu_s').set('probename', 'Nu_s');
        model.component('comp_s').probe('var_Nu_s').set('expr', 'alpha_s * l_ref_s / lambda_r');
        model.component('comp_s').probe('var_Nu_s').set('descractive', true);
        model.component('comp_s').probe('var_Nu_s').set('descr', 'Nu of Smooth Wall');
        model.component('comp_s').probe('var_Nu_s').label('Global Variable Probe of Nu Number');

    model.component('comp_s').probe.create('bnd_p_in_s', 'Boundary');
        model.component('comp_s').probe('bnd_p_in_s').set('probename', 'p_in_s');
        model.component('comp_s').probe('bnd_p_in_s').selection.named('sel_inlet_c_s');
        model.component('comp_s').probe('bnd_p_in_s').set('expr', 'spf_s.p');
        model.component('comp_s').probe('bnd_p_in_s').set('unit', 'Pa');
        model.component('comp_s').probe('bnd_p_in_s').set('descractive', true);
        model.component('comp_s').probe('bnd_p_in_s').set('descr', 'Avarage Inlet Pressure of Smooth Wall');
        model.component('comp_s').probe('bnd_p_in_s').label('Boundary Probe of Avarage Inlet Pressure');

    model.component('comp_s').probe.create('bnd_p_out_s', 'Boundary');
        model.component('comp_s').probe('bnd_p_out_s').set('probename', 'p_out_s');
        model.component('comp_s').probe('bnd_p_out_s').selection.named('sel_outlet_c_s');
        model.component('comp_s').probe('bnd_p_out_s').set('expr', 'spf_s.p');
        model.component('comp_s').probe('bnd_p_out_s').set('unit', 'Pa');
        model.component('comp_s').probe('bnd_p_out_s').set('descractive', true);
        model.component('comp_s').probe('bnd_p_out_s').set('descr', 'Outlet Pressure of Smooth Wall');
        model.component('comp_s').probe('bnd_p_out_s').label('Boundary Probe of Avarage Outlet Pressure');

    model.component('comp_s').probe.create('dom_p_dyn_s', 'Domain');
        model.component('comp_s').probe('dom_p_dyn_s').selection.named('sel_c_d_s');
        model.component('comp_s').probe('dom_p_dyn_s').set('probename', 'p_dyn_s');
        model.component('comp_s').probe('dom_p_dyn_s').set('expr', '0.5 * spf_s.rho * spf_s.U');
        model.component('comp_s').probe('dom_p_dyn_s').set('unit', 'Pa');
        model.component('comp_s').probe('dom_p_dyn_s').set('descractive', true);
        model.component('comp_s').probe('dom_p_dyn_s').set('descr', 'Avarage Dynamic Pressure of Smooth Wall');
        model.component('comp_s').probe('dom_p_dyn_s').label('Boundary Probe of Avarage Dynamic Pressure');

    model.component('comp_s').probe.create('var_zeta_s', 'GlobalVariable');
        model.component('comp_s').probe('var_zeta_s').set('probename', 'zeta_s');
        model.component('comp_s').probe('var_zeta_s').set('expr', 'abs(p_out_s - p_in_s) / p_dyn_s');
        model.component('comp_s').probe('var_zeta_s').set('unit', '1');
        model.component('comp_s').probe('var_zeta_s').set('descractive', true);
        model.component('comp_s').probe('var_zeta_s').set('descr', 'Zeta of Smooth Wall');
        model.component('comp_s').probe('var_zeta_s').label('Global Variable Probe of Zeta');
    
end

function model = build_component_rough (model)

% GEOMETRY

model.component.create('comp_r', true);

    model.component('comp_r').label('Rough Wall');

    model.component('comp_r').geom.create('geom_r', 3);

    model.component('comp_r').geom('geom_r').lengthUnit('mm');

    model.component('comp_r').geom('geom_r').create('cyl1', 'Cylinder');
        model.component('comp_r').geom('geom_r').feature('cyl1').set('r', 'd / 2');
        model.component('comp_r').geom('geom_r').feature('cyl1').set('h', 'l_in');

    model.component('comp_r').geom('geom_r').create('cyl2', 'Cylinder');
        model.component('comp_r').geom('geom_r').feature('cyl2').set('r', 'd / 2');
        model.component('comp_r').geom('geom_r').feature('cyl2').set('h', 'w_r');
        model.component('comp_r').geom('geom_r').feature('cyl2').set('pos', {'0' '0' 'l_in'});

    model.component('comp_r').geom('geom_r').create('cyl3', 'Cylinder');
        model.component('comp_r').geom('geom_r').feature('cyl3').set('r', 'd / 2 - h_r');
        model.component('comp_r').geom('geom_r').feature('cyl3').set('h', 'w_r');
        model.component('comp_r').geom('geom_r').feature('cyl3').set('pos', {'0' '0' 'l_in'});

    model.component('comp_r').geom('geom_r').create('dif1', 'Difference');
        model.component('comp_r').geom('geom_r').feature('dif1').selection('input').set({'cyl2'});
        model.component('comp_r').geom('geom_r').feature('dif1').selection('input2').set({'cyl3'});
        model.component('comp_r').geom('geom_r').feature('dif1').set('keep', true);

    model.component('comp_r').geom('geom_r').feature.create('del1', 'Delete');
        model.component('comp_r').geom('geom_r').feature('del1').selection('input').init(3);
        model.component('comp_r').geom('geom_r').feature('del1').selection('input').set('cyl2', 1);

    model.component('comp_r').geom('geom_r').create('cyl4', 'Cylinder');
        model.component('comp_r').geom('geom_r').feature('cyl4').set('r', 'd / 2');
        model.component('comp_r').geom('geom_r').feature('cyl4').set('h', 's_r');
        model.component('comp_r').geom('geom_r').feature('cyl4').set('pos', {'0' '0' 'l_in + w_r'});

    model.component('comp_r').geom('geom_r').create('arr1', 'Array');
        model.component('comp_r').geom('geom_r').feature('arr1').selection('input').set({'cyl3' 'cyl4' 'dif1'});
        model.component('comp_r').geom('geom_r').feature('arr1').set('fullsize', {'1' '1' 'n'});
        model.component('comp_r').geom('geom_r').feature('arr1').set('displ', {'0' '0' 'w_r + s_r'});

    model.component('comp_r').geom('geom_r').create('cyl5', 'Cylinder');
        model.component('comp_r').geom('geom_r').feature('cyl5').set('r', 'd / 2');
        model.component('comp_r').geom('geom_r').feature('cyl5').set('h', 'l_out');
        model.component('comp_r').geom('geom_r').feature('cyl5').set('pos', {'0' '0' 'l_in + n * (w_r + s_r)'});

    model.component('comp_r').geom('geom_r').run('fin');

% SELECTIONS

    model.component('comp_r').selection.create('sel_task_r', 'Cylinder');
        model.component('comp_r').selection('sel_task_r').label('Task');
        model.component('comp_r').selection('sel_task_r').set('r', 'd / 2');

    model.component('comp_r').selection.create('sel_fluid_d_r', 'Cylinder');
        model.component('comp_r').selection('sel_fluid_d_r').label('Fluid');

    model.component('comp_r').selection.create('sel_ribs_d', 'Difference');
        model.component('comp_r').selection('sel_ribs_d').label('Ribs');        
        model.component('comp_r').selection('sel_ribs_d').set('add', {'sel_task_r'});
        model.component('comp_r').selection('sel_ribs_d').set('subtract', {'sel_fluid_d_r'});

    model.component('comp_r').selection.create('sel_fluid_b_r', 'Adjacent');
        model.component('comp_r').selection('sel_fluid_b_r').label('Boundaries of Fluid');
        model.component('comp_r').selection('sel_fluid_b_r').set('input', {'sel_fluid_d_r'});

    model.component('comp_r').selection.create('sel_ribs_d_b_r', 'Adjacent');
        model.component('comp_r').selection('sel_ribs_d_b_r').label('Boundaries of Ribs');
        model.component('comp_r').selection('sel_ribs_d_b_r').set('input', {'sel_ribs_d'});

    model.component('comp_r').selection.create('sel_in_d_r', 'Cylinder');
        model.component('comp_r').selection('sel_in_d_r').set('top', '0');
        model.component('comp_r').selection('sel_in_d_r').set('bottom', '0');
        model.component('comp_r').selection('sel_in_d_r').label('Domain of Input');

    model.component('comp_r').selection.create('sel_out_d_r', 'Cylinder');
        model.component('comp_r').selection('sel_out_d_r').set('top', 'l_in + n * (w_r + h_r) + l_out');
        model.component('comp_r').selection('sel_out_d_r').set('bottom', 'l_in + n * (w_r + h_r) + l_out');
        model.component('comp_r').selection('sel_out_d_r').label('Domain of Output');

    model.component('comp_r').selection.create('sel_c_d_r', 'Difference');
        model.component('comp_r').selection('sel_c_d_r').set('add', {'sel_fluid_d_r'});
        model.component('comp_r').selection('sel_c_d_r').set('subtract', {'sel_in_d_r', 'sel_out_d_r'});
        model.component('comp_r').selection('sel_c_d_r').label('Domain of Center');

    model.component('comp_r').selection.create('sel_in_b_r', 'Adjacent');
        model.component('comp_r').selection('sel_in_b_r').label('Boundaries of Input');
        model.component('comp_r').selection('sel_in_b_r').set('input', {'sel_in_d_r'});

    model.component('comp_r').selection.create('sel_c_b_r', 'Adjacent');
        model.component('comp_r').selection('sel_c_b_r').label('Boundaries of Center');
        model.component('comp_r').selection('sel_c_b_r').set('input', {'sel_c_d_r'});

    model.component('comp_r').selection.create('sel_out_b_r', 'Adjacent');
        model.component('comp_r').selection('sel_out_b_r').label('Boundaries of Output');
        model.component('comp_r').selection('sel_out_b_r').set('input', {'sel_out_d_r'});

    model.component('comp_r').selection.create('sel_inlet_r', 'Cylinder');
        model.component('comp_r').selection('sel_inlet_r').set('entitydim', 2);
        model.component('comp_r').selection('sel_inlet_r').set('top', 0);
        model.component('comp_r').selection('sel_inlet_r').set('bottom', 0);
        model.component('comp_r').selection('sel_inlet_r').label('Inlet');

    model.component('comp_r').selection.create('sel_outlet_r', 'Cylinder');
        model.component('comp_r').selection('sel_outlet_r').set('entitydim', 2);
        model.component('comp_r').selection('sel_outlet_r').set('top', 'l_in + n * (w_r + s_r) + l_out');
        model.component('comp_r').selection('sel_outlet_r').set('bottom', 'l_in + n * (w_r + s_r) + l_out');
        model.component('comp_r').selection('sel_outlet_r').label('Outlet');

    model.component('comp_r').selection.create('sel_inlet_c_r', 'Cylinder');
        model.component('comp_r').selection('sel_inlet_c_r').set('entitydim', 2);
        model.component('comp_r').selection('sel_inlet_c_r').set('top', 'l_in');
        model.component('comp_r').selection('sel_inlet_c_r').set('bottom', 'l_in');
        model.component('comp_r').selection('sel_inlet_c_r').label('Inlet of Center');

    model.component('comp_r').selection.create('sel_outlet_c_r', 'Cylinder');
        model.component('comp_r').selection('sel_outlet_c_r').set('entitydim', 2);
        model.component('comp_r').selection('sel_outlet_c_r').set('top', 'l_in + n * (w_r + s_r)');
        model.component('comp_r').selection('sel_outlet_c_r').set('bottom', 'l_in + n * (w_r + s_r)');      
        model.component('comp_r').selection('sel_outlet_c_r').label('Outlet of Center');

    model.component('comp_r').selection.create('sel_int_1', 'Intersection');
        model.component('comp_r').selection('sel_int_1').set('entitydim', 2);
        model.component('comp_r').selection('sel_int_1').set('input', {'sel_ribs_d_b_r' 'sel_c_b_r'});
        model.component('comp_r').selection('sel_int_1').label('Intersection of Boundaries of Ribs & Center');

    model.component('comp_r').selection.create('sel_dif', 'Difference');
        model.component('comp_r').selection('sel_dif').set('entitydim', 2);
        model.component('comp_r').selection('sel_dif').set('add', {'sel_ribs_d_b_r'});
        model.component('comp_r').selection('sel_dif').set('subtract', {'sel_int_1', 'sel_in_b_r', 'sel_c_b_r'});
        model.component('comp_r').selection('sel_dif').label('Boundaries of Ribs for Condition');

    model.component('comp_r').selection.create('sel_int_2', 'Intersection');
        model.component('comp_r').selection('sel_int_2').set('entitydim', 2);
        model.component('comp_r').selection('sel_int_2').set('input', {'sel_fluid_b_r' 'sel_c_b_r'});
        model.component('comp_r').selection('sel_int_2').label('Intersection of Boundaries of Fluid & Center');

    model.component('comp_r').selection.create('sel_c_bc_r', 'Difference');
        model.component('comp_r').selection('sel_c_bc_r').set('entitydim', 2);
        model.component('comp_r').selection('sel_c_bc_r').set('add', {'sel_int_2'});
        model.component('comp_r').selection('sel_c_bc_r').set('subtract', {'sel_ribs_d_b_r'});
        model.component('comp_r').selection('sel_c_bc_r').label('Boundaries of Center for Condition');

    model.component('comp_r').selection.create('sel_in_wall_r', 'Difference');
        model.component('comp_r').selection('sel_in_wall_r').set('entitydim', 2);
        model.component('comp_r').selection('sel_in_wall_r').set('add', {'sel_in_b_r'});
        model.component('comp_r').selection('sel_in_wall_r').set('subtract', {'sel_inlet_r' 'sel_ribs_d_b_r' 'sel_inlet_c_r'});
        model.component('comp_r').selection('sel_in_wall_r').label('Wall of Input');

    model.component('comp_r').selection.create('sel_out_wall_r', 'Difference');
        model.component('comp_r').selection('sel_out_wall_r').set('entitydim', 2);
        model.component('comp_r').selection('sel_out_wall_r').set('add', {'sel_out_b_r'});
        model.component('comp_r').selection('sel_out_wall_r').set('subtract', {'sel_outlet_r' 'sel_outlet_c_r'});
        model.component('comp_r').selection('sel_out_wall_r').label('Wall of Output');

    model.component('comp_r').selection.create('sel_wall_r', 'Difference');
        model.component('comp_r').selection('sel_wall_r').set('entitydim', 2);
        model.component('comp_r').selection('sel_wall_r').set('add', {'sel_fluid_b_r'});
        model.component('comp_r').selection('sel_wall_r').set('subtract', {'sel_inlet_r' 'sel_outlet_r'});
        model.component('comp_r').selection('sel_wall_r').label('Wall');

    model.component('comp_r').selection.create('sel_bc_r', 'Union');
        model.component('comp_r').selection('sel_bc_r').set('entitydim', 2);
        model.component('comp_r').selection('sel_bc_r').set('input', {'sel_dif' 'sel_c_bc_r' 'sel_in_wall_r' 'sel_out_wall_r'});
        model.component('comp_r').selection('sel_bc_r').label('Boundaries of Condition');
        
    model.component('comp_r').selection.create('sel_c_wall_r', 'Difference');
        model.component('comp_r').selection('sel_c_wall_r').label('Wall of Center');
        model.component('comp_r').selection('sel_c_wall_r').set('entitydim', 2);
        model.component('comp_r').selection('sel_c_wall_r').set('add', {'sel_c_b_r'});
        model.component('comp_r').selection('sel_c_wall_r').set('subtract', {'sel_inlet_c_r' 'sel_outlet_c_r'});

% PHYSICS

    model.component('comp_r').physics.create('spf_r', 'LESRBVM', 'geom_s');

        model.component('comp_r').physics('spf_r').selection.named('sel_fluid_d_r');

        model.component('comp_r').physics('spf_r').feature('init1').set('u_init', {'0' '0' 'v_in'});
        model.component('comp_r').physics('spf_r').feature('init1').set('p_init', 'p_0');

        model.component('comp_r').physics('spf_r').create('inl1', 'InletBoundary', 2);
        model.component('comp_r').physics('spf_r').feature('inl1').set('U0in', 'v_in');
            model.component('comp_r').physics('spf_r').feature('inl1').selection.named('sel_inlet_r');

        model.component('comp_r').physics('spf_r').create('out1', 'OutletBoundary', 2);
            model.component('comp_r').physics('spf_r').feature('out1').set('p0', 'p_out');
            model.component('comp_r').physics('spf_r').feature('out1').selection.named('sel_outlet_r');

        model.component('comp_r').physics('spf_r').tag('spf_r');

    model.component('comp_r').physics.create('ht_r', 'HeatTransferInSolidsAndFluids', 'geom_s');

        model.component('comp_r').physics('ht_r').selection.named('sel_task_r');

        model.component('comp_r').physics('ht_r').feature('init1').set('Tinit', 'T_0');

        model.component('comp_r').physics('ht_r').feature('fluid1').selection.named('sel_fluid_d_r');

        model.component('comp_r').physics('ht_r').create('ifl1', 'Inflow', 2);
            model.component('comp_r').physics('ht_r').feature('ifl1').set('Tustr', 'T_in');
            model.component('comp_r').physics('ht_r').feature('ifl1').selection.named('sel_inlet_r');

        model.component('comp_r').physics('ht_r').create('ofl1', 'ConvectiveOutflow', 2);
            model.component('comp_r').physics('ht_r').feature('ofl1').selection.named('sel_outlet_r');

        model.component('comp_r').physics('ht_r').create('temp1', 'TemperatureBoundary', 2);
            model.component('comp_r').physics('ht_r').feature('temp1').set('T0', 'T_w');
            model.component('comp_r').physics('ht_r').feature('temp1').selection.named('sel_bc_r');

        model.component('comp_r').physics('ht_r').tag('ht_r');

    model.component('comp_r').multiphysics.create('nitf_r', 'NonIsothermalFlow', 3);

% PROBES

    model.component('comp_r').probe.create('dom_int_1_r', 'Domain');
        model.component('comp_r').probe('dom_int_1_r').selection.named('sel_c_d_r');
        model.component('comp_r').probe('dom_int_1_r').set('type', 'integral');
        model.component('comp_r').probe('dom_int_1_r').set('expr', 'spf_r.rho * spf_r.U');
        model.component('comp_r').probe('dom_int_1_r').set('probename', 'int_rho_v');
        model.component('comp_r').probe('dom_int_1_r').set('unit', 'kg * m / s');
        model.component('comp_r').probe('dom_int_1_r').label('Domain Probe of Density Mass Rate');

    model.component('comp_r').probe.create('dom_int_2_r', 'Domain');
        model.component('comp_r').probe('dom_int_2_r').selection.named('sel_c_d_r');
        model.component('comp_r').probe('dom_int_2_r').set('type', 'integral');
        model.component('comp_r').probe('dom_int_2_r').set('expr', 'spf_r.rho * spf_r.U * ht_r.T');
        model.component('comp_r').probe('dom_int_2_r').set('probename', 'int_rho_v_T');
        model.component('comp_r').probe('dom_int_2_r').set('unit', 'kg * m * K / s');
        model.component('comp_r').probe('dom_int_2_r').label('Domain Probe of Density rho*v*T');

    model.component('comp_r').probe.create('var_T_m_r', 'GlobalVariable');
        model.component('comp_r').probe('var_T_m_r').set('probename', 'T_r_m');
        model.component('comp_r').probe('var_T_m_r').set('expr', 'int_rho_v_T / int_rho_v');
        model.component('comp_r').probe('var_T_m_r').set('unit', 'K');
        model.component('comp_r').probe('var_T_m_r').set('descractive', true);
        model.component('comp_r').probe('var_T_m_r').set('descr', 'Mass Average Temperature of Rough Wall');
        model.component('comp_r').probe('var_T_m_r').label('Domain Probe of Mass Average Temperature');

    model.component('comp_r').probe.create('bnd_alpha_r', 'Boundary');
        model.component('comp_r').probe('bnd_alpha_r').selection.named('sel_c_wall_r');
        model.component('comp_r').probe('bnd_alpha_r').set('probename', 'alpha_r');
        model.component('comp_r').probe('bnd_alpha_r').set('expr', 'abs(ht_r.ntflux/(T_r_m - T_w))');
        model.component('comp_r').probe('bnd_alpha_r').set('descractive', true);
        model.component('comp_r').probe('bnd_alpha_r').set('descr', 'Average Alpha of Rough Wall');
        model.component('comp_r').probe('bnd_alpha_r').label('Boundary Probe of Average Alpha');

    model.component('comp_r').probe.create('var_l_ref_r', 'GlobalVariable');
        model.component('comp_r').probe('var_l_ref_r').set('probename', 'l_ref_r');
        model.component('comp_r').probe('var_l_ref_r').set('expr', 'root.l_ref_r');
        model.component('comp_r').probe('var_l_ref_r').set('unit', 'm');
        model.component('comp_r').probe('var_l_ref_r').set('descractive', true);
        model.component('comp_r').probe('var_l_ref_r').set('descr', 'Reference Length of Rough Wall');
        model.component('comp_r').probe('var_l_ref_r').set('expr', 'root.l_ref_r');
        model.component('comp_r').probe('var_l_ref_r').label('Global Variable Probe of Reference Length');

    model.component('comp_r').probe.create('bnd_lambda_r', 'Boundary');
        model.component('comp_r').probe('bnd_lambda_r').selection.named('sel_c_wall_r');
        model.component('comp_r').probe('bnd_lambda_r').set('expr', 'ht_r.k_iso');
        model.component('comp_r').probe('bnd_lambda_r').set('probename', 'lambda_r');
        model.component('comp_r').probe('bnd_lambda_r').set('descractive', true);
        model.component('comp_r').probe('bnd_lambda_r').set('descr', 'Avarage Thermal Conductivity of Rough Wall');
        model.component('comp_r').probe('bnd_lambda_r').label('Boundary Probe of Avarege Lambda');

    model.component('comp_r').probe.create('var_Nu_r', 'GlobalVariable');
        model.component('comp_r').probe('var_Nu_r').set('probename', 'Nu_r');
        model.component('comp_r').probe('var_Nu_r').set('expr', 'alpha_r * l_ref_r / lambda_r');
        model.component('comp_r').probe('var_Nu_r').set('descractive', true);
        model.component('comp_r').probe('var_Nu_r').set('descr', 'Nu of Rough Wall');
        model.component('comp_r').probe('var_Nu_r').label('Global Variable Probe of Nu Number');

    model.component('comp_r').probe.create('bnd_p_in_r', 'Boundary');
        model.component('comp_r').probe('bnd_p_in_r').set('probename', 'p_in_r');
        model.component('comp_r').probe('bnd_p_in_r').selection.named('sel_inlet_c_r');
        model.component('comp_r').probe('bnd_p_in_r').set('expr', 'spf_r.p');
        model.component('comp_r').probe('bnd_p_in_r').set('unit', 'Pa');
        model.component('comp_r').probe('bnd_p_in_r').set('descractive', true);
        model.component('comp_r').probe('bnd_p_in_r').set('descr', 'Avarage Inlet Pressure of Rough Wall');
        model.component('comp_r').probe('bnd_p_in_r').label('Boundary Probe of Avarage Inlet Pressure');

    model.component('comp_r').probe.create('bnd_p_out_r', 'Boundary');
        model.component('comp_r').probe('bnd_p_out_r').set('probename', 'p_out_r');
        model.component('comp_r').probe('bnd_p_out_r').selection.named('sel_outlet_c_r');
        model.component('comp_r').probe('bnd_p_out_r').set('expr', 'spf_r.p');
        model.component('comp_r').probe('bnd_p_out_r').set('unit', 'Pa');
        model.component('comp_r').probe('bnd_p_out_r').set('descractive', true);
        model.component('comp_r').probe('bnd_p_out_r').set('descr', 'Outlet Pressure of Rough Wall');
        model.component('comp_r').probe('bnd_p_out_r').label('Boundary Probe of Avarage Outlet Pressure');

    model.component('comp_r').probe.create('dom_p_dyn_r', 'Domain');
        model.component('comp_r').probe('dom_p_dyn_r').selection.named('sel_c_d_r');
        model.component('comp_r').probe('dom_p_dyn_r').set('probename', 'p_dyn_r');
        model.component('comp_r').probe('dom_p_dyn_r').set('expr', '0.5 * spf_r.rho * spf_r.U');
        model.component('comp_r').probe('dom_p_dyn_r').set('unit', 'Pa');
        model.component('comp_r').probe('dom_p_dyn_r').set('descractive', true);
        model.component('comp_r').probe('dom_p_dyn_r').set('descr', 'Avarage Dynamic Pressure of Rough Wall');
        model.component('comp_r').probe('dom_p_dyn_r').label('Boundary Probe of Avarage Dynamic Pressure');

    model.component('comp_r').probe.create('var_zeta_r', 'GlobalVariable');
        model.component('comp_r').probe('var_zeta_r').set('probename', 'zeta_r');
        model.component('comp_r').probe('var_zeta_r').set('expr', 'abs(p_out_r - p_in_r) / p_dyn_r');
        model.component('comp_r').probe('var_zeta_r').set('unit', '1');
        model.component('comp_r').probe('var_zeta_r').set('descractive', true);
        model.component('comp_r').probe('var_zeta_r').set('descr', 'Zeta of Rough Wall');
        model.component('comp_r').probe('var_zeta_r').label('Global Variable Probe of Zeta');

end