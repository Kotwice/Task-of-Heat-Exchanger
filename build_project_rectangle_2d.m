
% PARAMETERS LIST

parameters.d_ht = 10;
parameters.l_st_in = 10 * parameters.d_ht;
parameters.l_st_out = 10 * parameters.d_ht;
parameters.h_r = 0.5;
parameters.s_r = 5;
parameters.w_r = 1;
parameters.n = 10;
parameters.Re = 1e4;
parameters.p_0 = 0;
parameters.T_wall = 1e3;
parameters.T_0 = 298.15;
parameters.T_in = 300;

% EXECUTION CODE

model = create_model(parameters);
model = build_component_smooth(model);
model = build_component_rough(model, parameters);
model = set_studies(model);

mphsave(model, strcat(pwd, '\build_project_rectangle.mph'));

% INITIALIZATION BLOCK

function model = create_model(parameters)

    import com.comsol.model.*
    import com.comsol.model.util.*

    model = ModelUtil.create('Model');

        model.param.set('rho_air', '1.1839 [kg/m^3]', 'Density of Air at 25 Cdeg');
        model.param.set('mu_air', '18.35 [uPa*s]', 'Dynamic Viscosity of Air at 25 Cdeg');

        model.param.set('Re', strcat(num2str(parameters.Re)), 'Reynolds Number');
        model.param.set('d_ht', strcat(num2str(parameters.d_ht), '[mm]'), 'Diameter of Tube');
        model.param.set('h_r', strcat(num2str(parameters.h_r), '[mm]'), 'High of Roughness');
        model.param.set('s_r', strcat(num2str(parameters.s_r), '[mm]'), 'Step Length of Roughness');
        model.param.set('w_r', strcat(num2str(parameters.w_r), '[mm]'), 'Width of Roughness');
        model.param.set('v_z_in', 'mu_air * Re / (d_ht * rho_air)', 'Inlet Velocity');
        model.param.set('p_0', strcat(num2str(parameters.p_0), '[Pa]'), 'Initial Pressure');
        model.param.set('T_wall', strcat(num2str(parameters.T_wall), '[K]'), 'Temperature of Wall');
        model.param.set('T_0', strcat(num2str(parameters.T_0), '[K]'), 'Initial Temperature');
        model.param.set('T_in', strcat(num2str(parameters.T_in ), '[K]'), 'Inlet Temperature');

        model.param.set('l_st_in', strcat(num2str(parameters.l_st_in), ' [mm]'), 'Length of Input Stabilization Domain');
        model.param.set('n', strcat(num2str(parameters.n)), 'Number of Ribs'); 
        model.param.set('l_ht', 'n * (s_r + w_r)', 'Length of Heat Transfer Domain');
        model.param.set('l_st_out', strcat(num2str(parameters.l_st_out), ' [mm]'), 'Length of Output Stabilization Domain');

        model.param.set('h_r_min', '0.05 * d_ht', 'Minimal Value of High of Roughness');
        model.param.set('h_r_step', '1 [mm]', 'Step Value of High of Roughness');
        model.param.set('h_r_max', '0.3 * d_ht', 'Maximal Value of High of Roughness');

        model.param.set('s_r_min', '0.25 * d_ht', 'Minimal Value of Step Length of Roughness');
        model.param.set('s_r_step', '1 [mm]', 'Step Value of Step Length of Roughness');
        model.param.set('s_r_max', 'd_ht', 'Maximal Step Length of High of Roughness');   

        model.param.set('Re_min', '1e4', 'Minimal Reynolds Number');
        model.param.set('Re_step', '5 * 1e3', 'Step Reynolds Number');
        model.param.set('Re_max', '1e5', 'Maximal Reynolds Number');  

end

function model = build_component_smooth (model)

    model.component.create('smooth');
 
    model.component('smooth').label('Smooth Wall');

    % BUILD GEOMETRY

    model.component('smooth').geom.create('geom_smooth', 2);

            model.component('smooth').geom('geom_smooth').axisymmetric(true);

            model.component('smooth').geom('geom_smooth').lengthUnit('mm');

            model.component('smooth').geom('geom_smooth').feature.create('r1', 'Rectangle');                    
                    
                    model.component('smooth').geom('geom_smooth').feature('r1').set('size', {'d_ht / 2' 'l_st_in'});

            model.component('smooth').geom('geom_smooth').create('r2', 'Rectangle');

                    model.component('smooth').geom('geom_smooth').feature('r2').set('size', {'d_ht / 2' 'n * (s_r + w_r)'});
                    model.component('smooth').geom('geom_smooth').feature('r2').set('pos', {'0' 'l_st_in'});

            model.component('smooth').geom('geom_smooth').create('r3', 'Rectangle');

                    model.component('smooth').geom('geom_smooth').feature('r3').set('size', {'d_ht / 2' 'l_st_out'});
                    model.component('smooth').geom('geom_smooth').feature('r3').set('pos', {'0' 'l_st_in + n * (s_r + w_r)'});

            model.component('smooth').geom('geom_smooth').run;
    
    % APPOINTMETN OF SELECTION NAMES

    model.component('smooth').selection.create('selection_in_st_s', 'Box');

            model.component('smooth').selection('selection_in_st_s').set('xmin', 0);
            model.component('smooth').selection('selection_in_st_s').set('xmax', 'd_ht / 8');
            model.component('smooth').selection('selection_in_st_s').set('ymin', '0.01 * l_st_in');
            model.component('smooth').selection('selection_in_st_s').set('ymax', '0.99 * l_st_in');

            model.component('smooth').selection('selection_in_st_s').label('Input Stabilization Domain');

    model.component('smooth').selection.create('selection_ht_s', 'Box');

            model.component('smooth').selection('selection_ht_s').set('xmin', 0);
            model.component('smooth').selection('selection_ht_s').set('xmax', 'd_ht / 8');
            model.component('smooth').selection('selection_ht_s').set('ymin', 'l_st_in + 0.01 * n * (s_r + w_r)');
            model.component('smooth').selection('selection_ht_s').set('ymax', 'l_st_in + 0.99 * n * (s_r + w_r)');

            model.component('smooth').selection('selection_ht_s').label('Heat Transfer Domain');

    model.component('smooth').selection.create('selection_out_st_s', 'Box');

            model.component('smooth').selection('selection_out_st_s').set('xmin', 0);
            model.component('smooth').selection('selection_out_st_s').set('xmax', 'd_ht / 8');
            model.component('smooth').selection('selection_out_st_s').set('ymin', 'l_st_in + n * (s_r + w_r) + 0.01 * l_st_out');
            model.component('smooth').selection('selection_out_st_s').set('ymax', 'l_st_in + n * (s_r + w_r) + 0.99 * l_st_out');

            model.component('smooth').selection('selection_out_st_s').label('Output Stabilization Domain');

    model.component('smooth').selection.create('selection_fluid_s', 'Union');

            model.component('smooth').selection('selection_fluid_s').set('input', {'selection_in_st_s' 'selection_ht_s' 'selection_out_st_s'});
            model.component('smooth').selection('selection_fluid_s').label('Fluid');

    model.component('smooth').selection.create('selection_axis_in_st_s', 'Box');

            model.component('smooth').selection('selection_axis_in_st_s').set('xmin', 0);
            model.component('smooth').selection('selection_axis_in_st_s').set('xmax', 'd_ht / 8');
            model.component('smooth').selection('selection_axis_in_st_s').set('ymin', '0.01 * l_st_in');
            model.component('smooth').selection('selection_axis_in_st_s').set('ymax', '0.99 * l_st_in');

            model.component('smooth').selection('selection_axis_in_st_s').label('Axis of Input Stabilization Domain');

            model.component('smooth').selection('selection_axis_in_st_s').set('entitydim', 1);

    model.component('smooth').selection.create('selection_axis_ht_s', 'Box');

            model.component('smooth').selection('selection_axis_ht_s').set('xmin', 0);
            model.component('smooth').selection('selection_axis_ht_s').set('xmax', 'd_ht / 8');
            model.component('smooth').selection('selection_axis_ht_s').set('ymin', 'l_st_in + 0.01 * n * (s_r + w_r)');
            model.component('smooth').selection('selection_axis_ht_s').set('ymax', 'l_st_in + 0.99 * n * (s_r + w_r)');

            model.component('smooth').selection('selection_axis_ht_s').label('Axis of Heat Transfer Domain');

            model.component('smooth').selection('selection_axis_ht_s').set('entitydim', 1);

    model.component('smooth').selection.create('selection_axis_out_st_s', 'Box');

            model.component('smooth').selection('selection_axis_out_st_s').set('xmin', 0);
            model.component('smooth').selection('selection_axis_out_st_s').set('xmax', 'd_ht / 4');
            model.component('smooth').selection('selection_axis_out_st_s').set('ymin', 'l_st_in + n * (s_r + w_r) + 0.01 * l_st_out');
            model.component('smooth').selection('selection_axis_out_st_s').set('ymax', 'l_st_in + n * (s_r + w_r) + 0.99 * l_st_out');

            model.component('smooth').selection('selection_axis_out_st_s').label('Axis of Output Stabilization Domain');

            model.component('smooth').selection('selection_axis_out_st_s').set('entitydim', 1);

    model.component('smooth').selection.create('selection_axis_s', 'Union');

            model.component('smooth').selection('selection_axis_s').set('entitydim', 1);
            model.component('smooth').selection('selection_axis_s').set('input', {'selection_axis_in_st_s' 'selection_axis_ht_s' 'selection_axis_out_st_s'});
            model.component('smooth').selection('selection_axis_s').label('Axis');

    model.component('smooth').selection.create('selection_inlet_s', 'Box');

            model.component('smooth').selection('selection_inlet_s').set('xmin', 'd_ht / 16');
            model.component('smooth').selection('selection_inlet_s').set('xmax', 'd_ht / 8');
            model.component('smooth').selection('selection_inlet_s').set('ymin', '0');
            model.component('smooth').selection('selection_inlet_s').set('ymax', '0');

            model.component('smooth').selection('selection_inlet_s').label('Inlet');

            model.component('smooth').selection('selection_inlet_s').set('entitydim', 1);

    model.component('smooth').selection.create('selection_inlet_ht_s', 'Box');

            model.component('smooth').selection('selection_inlet_ht_s').set('xmin', 'd_ht / 16');
            model.component('smooth').selection('selection_inlet_ht_s').set('xmax', 'd_ht / 8');
            model.component('smooth').selection('selection_inlet_ht_s').set('ymin', 'l_st_in');
            model.component('smooth').selection('selection_inlet_ht_s').set('ymax', 'l_st_in');

            model.component('smooth').selection('selection_inlet_ht_s').label('Inlet Heat Transfer Domain');

            model.component('smooth').selection('selection_inlet_ht_s').set('entitydim', 1);

    model.component('smooth').selection.create('selection_outlet_s', 'Box');

            model.component('smooth').selection('selection_outlet_s').set('xmin', 'd_ht / 16');
            model.component('smooth').selection('selection_outlet_s').set('xmax', 'd_ht / 8');
            model.component('smooth').selection('selection_outlet_s').set('ymin', 'l_st_in + n * (s_r + w_r) + l_st_out');
            model.component('smooth').selection('selection_outlet_s').set('ymax', 'l_st_in + n * (s_r + w_r) + l_st_out');

            model.component('smooth').selection('selection_outlet_s').label('Outlet');

            model.component('smooth').selection('selection_outlet_s').set('entitydim', 1);

    model.component('smooth').selection.create('selection_outlet_ht_s', 'Box');

            model.component('smooth').selection('selection_outlet_ht_s').set('xmin', 'd_ht / 16');
            model.component('smooth').selection('selection_outlet_ht_s').set('xmax', 'd_ht / 8');
            model.component('smooth').selection('selection_outlet_ht_s').set('ymin', 'l_st_in + n * (s_r + w_r)');
            model.component('smooth').selection('selection_outlet_ht_s').set('ymax', 'l_st_in + n * (s_r + w_r)');

            model.component('smooth').selection('selection_outlet_ht_s').label('Outlet Heat Transfer Domain');

            model.component('smooth').selection('selection_outlet_ht_s').set('entitydim', 1);

    model.component('smooth').selection.create('selection_wall_st_in_s', 'Box');

            model.component('smooth').selection('selection_wall_st_in_s').set('xmin', 'd_ht / 2');
            model.component('smooth').selection('selection_wall_st_in_s').set('xmax', 'd_ht / 2');
            model.component('smooth').selection('selection_wall_st_in_s').set('ymin', '0.01 * l_st_in');
            model.component('smooth').selection('selection_wall_st_in_s').set('ymax', '0.99 * l_st_in');

            model.component('smooth').selection('selection_wall_st_in_s').label('Wall of Input Stabilization Domain');

            model.component('smooth').selection('selection_wall_st_in_s').set('entitydim', 1);

    model.component('smooth').selection.create('selection_wall_ht_s', 'Box');

            model.component('smooth').selection('selection_wall_ht_s').set('xmin', 'd_ht / 2');
            model.component('smooth').selection('selection_wall_ht_s').set('xmax', 'd_ht / 2');
            model.component('smooth').selection('selection_wall_ht_s').set('ymin', 'l_st_in + 0.01 * n * (s_r + w_r)');
            model.component('smooth').selection('selection_wall_ht_s').set('ymax', 'l_st_in + 0.99 * n * (s_r + w_r)');

            model.component('smooth').selection('selection_wall_ht_s').label('Wall of Heat Transfer');

            model.component('smooth').selection('selection_wall_ht_s').set('entitydim', 1);

    model.component('smooth').selection.create('selection_wall_st_out_s', 'Box');

            model.component('smooth').selection('selection_wall_st_out_s').set('xmin', 'd_ht / 2');
            model.component('smooth').selection('selection_wall_st_out_s').set('xmax', 'd_ht / 2');
            model.component('smooth').selection('selection_wall_st_out_s').set('ymin', 'l_st_in + n * (s_r + w_r) + 0.01 * l_st_out');
            model.component('smooth').selection('selection_wall_st_out_s').set('ymax', 'l_st_in + n * (s_r + w_r) + 0.99 * l_st_out');

            model.component('smooth').selection('selection_wall_st_out_s').label('Wall of Output Stabilization Domain');

            model.component('smooth').selection('selection_wall_st_out_s').set('entitydim', 1);

    model.component('smooth').selection.create('selection_wall_s', 'Union');

        model.component('smooth').selection('selection_wall_s').set('entitydim', 1);
        model.component('smooth').selection('selection_wall_s').set('input', {'selection_wall_st_in_s' 'selection_wall_ht_s' 'selection_wall_st_out_s'});
        model.component('smooth').selection('selection_wall_s').label('Wall of Fluid');
    
    model.component('smooth').selection.create('selection_distr_z_in', 'Union');
    
        model.component('smooth').selection('selection_distr_z_in').set('entitydim', 1);
        model.component('smooth').selection('selection_distr_z_in').set('input', {'selection_axis_in_st_s' 'selection_wall_st_in_s'});
    
        model.component('smooth').selection('selection_distr_z_in').label('Distribution z Input Domain');
    
    model.component('smooth').selection.create('selection_distr_r_in', 'Union');

        model.component('smooth').selection('selection_distr_r_in').set('entitydim', 1);
        model.component('smooth').selection('selection_distr_r_in').set('input', {'selection_inlet_s' 'selection_inlet_ht_s'});
    
        model.component('smooth').selection('selection_distr_r_in').label('Distribution r Input Domain');

    model.component('smooth').selection.create('selection_distr_z_ht', 'Union');
    
        model.component('smooth').selection('selection_distr_z_ht').label('Distribution z Heat Transfer Domain');
        model.component('smooth').selection('selection_distr_z_ht').set('entitydim', 1);
    
        model.component('smooth').selection('selection_distr_z_ht').set('input', {'selection_axis_ht_s' 'selection_wall_ht_s'});
    
    model.component('smooth').selection.create('selection_distr_r_ht', 'Union');
    
        model.component('smooth').selection('selection_distr_r_ht').set('entitydim', 1);
        model.component('smooth').selection('selection_distr_r_ht').set('input', {'selection_inlet_ht_s' 'selection_outlet_ht_s'});
    
        model.component('smooth').selection('selection_distr_r_ht').label('Distribution r Heat Transfer Domain');
    
    model.component('smooth').selection.create('selection_disrt_z_out', 'Union');
    
        model.component('smooth').selection('selection_disrt_z_out').set('entitydim', 1);
        model.component('smooth').selection('selection_disrt_z_out').set('input', {'selection_axis_out_st_s' 'selection_wall_st_out_s'});
    
        model.component('smooth').selection('selection_disrt_z_out').label('Distribution z Output Domain');
    
    model.component('smooth').selection.create('selection_disrt_r_out', 'Union');
    
        model.component('smooth').selection('selection_disrt_r_out').set('entitydim', 1);
        model.component('smooth').selection('selection_disrt_r_out').set('input', {'selection_outlet_s' 'selection_outlet_ht_s'});
    
        model.component('smooth').selection('selection_disrt_r_out').label('Distribution r Output Domain');

    model.component('smooth').selection.create('selection_distr_r', 'Union');
    
        model.component('smooth').selection('selection_distr_r').set('entitydim', 1);
        model.component('smooth').selection('selection_distr_r').set('input', {'selection_inlet_s' 'selection_inlet_ht_s', 'selection_outlet_s' 'selection_outlet_ht_s'});
    
        model.component('smooth').selection('selection_distr_r').label('Distribution r');

    % SET OF MATERIAL PROPORTIES

        % AIR PROPORTIES

        model.component('smooth').material.create('mat_smooth_air', 'Common');

                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func.create('eta', 'Piecewise');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func.create('Cp', 'Piecewise');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func.create('rho', 'Analytic');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func.create('k', 'Piecewise');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func.create('cs', 'Analytic');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func.create('an1', 'Analytic');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func.create('an2', 'Analytic');
                model.component('smooth').material('mat_smooth_air').propertyGroup.create('RefractiveIndex', 'Refractive index');
                model.component('smooth').material('mat_smooth_air').propertyGroup.create('NonlinearModel', 'Nonlinear model');
                model.component('smooth').material('mat_smooth_air').propertyGroup.create('idealGas', 'Ideal gas');
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').func.create('Cp', 'Piecewise');
                model.component('smooth').material('mat_smooth_air').label('Air');
                model.component('smooth').material('mat_smooth_air').set('family', 'air');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('eta').set('arg', 'T');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('eta').set('argunit', 'K');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('Cp').set('arg', 'T');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('Cp').set('argunit', 'K');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/R_const[K*mol/J]/T');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('rho').set('argunit', 'Pa,K');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('rho').set('plotargs', {'pA' '0' '1'; 'T' '0' '1'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('k').set('arg', 'T');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('k').set('argunit', 'K');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*R_const[K*mol/J]/0.02897*T)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('cs').set('args', {'T'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('cs').set('argunit', 'K');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('cs').set('fununit', 'm/s');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('cs').set('plotargs', {'T' '273.15' '373.15'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an1').set('expr', '-1/rho(pA,T)*d(rho(pA,T),T)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an1').set('args', {'pA' 'T'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an1').set('argunit', 'Pa,K');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an1').set('fununit', '1/K');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an1').set('plotargs', {'pA' '101325' '101325'; 'T' '273.15' '373.15'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an2').set('funcname', 'muB');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an2').set('expr', '0.6*eta(T)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an2').set('args', {'T'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an2').set('argunit', 'K');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an2').set('fununit', 'Pa*s');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').func('an2').set('plotargs', {'T' '200' '1600'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('thermalexpansioncoefficient', '');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('molarmass', '');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('bulkviscosity', '');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('molarmass', '0.02897[kg/mol]');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('bulkviscosity', 'muB(T)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('ratioofspecificheat', '1.4');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('electricconductivity', {'0[S/m]' '0' '0' '0' '0[S/m]' '0' '0' '0' '0[S/m]'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('heatcapacity', 'Cp(T)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('density', 'rho(pA,T)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('thermalconductivity', {'k(T)' '0' '0' '0' 'k(T)' '0' '0' '0' 'k(T)'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').set('soundspeed', 'cs(T)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').addInput('temperature');
                model.component('smooth').material('mat_smooth_air').propertyGroup('def').addInput('pressure');
                model.component('smooth').material('mat_smooth_air').propertyGroup('RefractiveIndex').set('n', '');
                model.component('smooth').material('mat_smooth_air').propertyGroup('RefractiveIndex').set('ki', '');
                model.component('smooth').material('mat_smooth_air').propertyGroup('RefractiveIndex').set('n', '');
                model.component('smooth').material('mat_smooth_air').propertyGroup('RefractiveIndex').set('ki', '');
                model.component('smooth').material('mat_smooth_air').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('NonlinearModel').set('BA', '(def.gamma+1)/2');
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').func('Cp').label('Piecewise 2');
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').func('Cp').set('arg', 'T');
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').func('Cp').set('argunit', 'K');
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').func('Cp').set('fununit', 'J/(kg*K)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').set('Rs', 'R_const/Mn');
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').set('heatcapacity', 'Cp(T)');
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').set('ratioofspecificheat', '1.4');
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').set('molarmass', '0.02897');
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').addInput('temperature');
                model.component('smooth').material('mat_smooth_air').propertyGroup('idealGas').addInput('pressure');
                model.component('smooth').material('mat_smooth_air').set('groups', {});
                model.component('smooth').material('mat_smooth_air').set('family', 'air');
                model.component('smooth').material('mat_smooth_air').selection.named('selection_fluid_s');

    % SET MULTIPHYSICS
            
        % SET SPF MODULE

        model.component('smooth').physics.create('spf_smooth', 'TurbulentFlowkeps', 'geom_smooth');

                model.component('smooth').physics('spf_smooth').selection.named('selection_fluid_s');

                model.component('smooth').physics('spf_smooth').feature('init1').set('u_init', {'0' '0' 'v_z_in'});

                model.component('smooth').physics('spf_smooth').create('inl1', 'InletBoundary', 1);

                        model.component('smooth').physics('spf_smooth').feature('inl1').selection.named('selection_inlet_s');
                        model.component('smooth').physics('spf_smooth').feature('inl1').set('U0in', 'v_z_in');

                model.component('smooth').physics('spf_smooth').create('out1', 'OutletBoundary', 1);

                        model.component('smooth').physics('spf_smooth').feature('out1').set('p0', 'p_0');
                        model.component('smooth').physics('spf_smooth').feature('out1').selection.named('selection_outlet_s');

        % SET Ht MODULE

        model.component('smooth').physics.create('ht_smooth', 'HeatTransferInFluids', 'geom_smooth');

                model.component('smooth').physics('ht_smooth').selection.named('selection_fluid_s');

                model.component('smooth').physics('ht_smooth').feature('init1').set('Tinit', 'T_0');

                model.component('smooth').physics('ht_smooth').create('temp1', 'TemperatureBoundary', 1);

                        model.component('smooth').physics('ht_smooth').feature('temp1').selection.named('selection_wall_s');
                        model.component('smooth').physics('ht_smooth').feature('temp1').set('T0', 'T_wall');

                model.component('smooth').physics('ht_smooth').create('ifl1', 'Inflow', 1);

                        model.component('smooth').physics('ht_smooth').feature('ifl1').selection.named('selection_inlet_s');
                        model.component('smooth').physics('ht_smooth').feature('ifl1').set('Tustr', 'T_0');

                model.component('smooth').physics('ht_smooth').create('ofl1', 'ConvectiveOutflow', 1);

                        model.component('smooth').physics('ht_smooth').feature('ofl1').selection.named('selection_outlet_s');
            
         model.component('smooth').multiphysics.create('nitf_smooth', 'NonIsothermalFlow', 2);        

    % SET PROBE VARIABLES

    model.component('smooth').probe.create('dom_mr_s', 'Domain');

        model.component('smooth').probe('dom_mr_s').set('intsurface', true);
        model.component('smooth').probe('dom_mr_s').set('intvolume', true);
        model.component('smooth').probe('dom_mr_s').label('Domain Probe of Mass Rate');
        model.component('smooth').probe('dom_mr_s').selection.named('selection_ht_s');
        model.component('smooth').probe('dom_mr_s').set('expr', 'spf.rho * spf.U');
        model.component('smooth').probe('dom_mr_s').set('probename', 'mr_s');

    model.component('smooth').probe.create('bnd_T_s', 'Domain');

        model.component('smooth').probe('bnd_T_s').set('intsurface', true);
        model.component('smooth').probe('bnd_T_s').set('intvolume', true);
        model.component('smooth').probe('bnd_T_s').set('expr', 'spf.rho * spf.U * T');
        model.component('smooth').probe('bnd_T_s').selection.named('selection_ht_s');
        model.component('smooth').probe('bnd_T_s').set('probename', 'T_s_int');
        model.component('smooth').probe('bnd_T_s').label('Domain Probe of Integral Temperature');

    model.component('smooth').probe.create('var_T_s', 'GlobalVariable');

        model.component('smooth').probe('var_T_s').label('Global Variable Probe of Mean Temperature');
        model.component('smooth').probe('var_T_s').set('probename', 'T_s');
        model.component('smooth').probe('var_T_s').set('expr', 'T_s_int / mr_s');
        model.component('smooth').probe('var_T_s').set('descractive', true);
        model.component('smooth').probe('var_T_s').set('descr', 'T_s');

    model.component('smooth').probe.create('bnd_alpha_s', 'Boundary');

            model.component('smooth').probe('bnd_alpha_s').set('intsurface', true);
            model.component('smooth').probe('bnd_alpha_s').label('Boundary Probe of Alpha Smooth Wall');
            model.component('smooth').probe('bnd_alpha_s').set('probename', 'alpha_s');
            model.component('smooth').probe('bnd_alpha_s').selection.named('selection_wall_ht_s');
            model.component('smooth').probe('bnd_alpha_s').set('expr', 'abs(ht.ntflux/(T_s - T_wall))');
            model.component('smooth').probe('bnd_alpha_s').set('descractive', true);
            model.component('smooth').probe('bnd_alpha_s').set('descr', 'alpha_s');

    model.component('smooth').probe.create('bnd_l_ref_s', 'Boundary');

            model.component('smooth').probe('bnd_l_ref_s').set('intsurface', true);
            model.component('smooth').probe('bnd_l_ref_s').label('Boundary Probe of Reference Length Smooth Wall');
            model.component('smooth').probe('bnd_l_ref_s').set('probename', 'l_ref_s');
            model.component('smooth').probe('bnd_l_ref_s').selection.named('selection_wall_ht_s');
            model.component('smooth').probe('bnd_l_ref_s').set('expr', '1 / (2 * (pi * d_ht + l_ht))');
            model.component('smooth').probe('bnd_l_ref_s').set('type', 'integral');
            model.component('smooth').probe('bnd_l_ref_s').set('unit', 'm');
            model.component('smooth').probe('bnd_l_ref_s').set('descractive', true);
            model.component('smooth').probe('bnd_l_ref_s').set('descr', 'l_ref_s');

    model.component('smooth').probe.create('bnd_lambda_s', 'Boundary');

            model.component('smooth').probe('bnd_lambda_s').set('intsurface', true);
            model.component('smooth').probe('bnd_lambda_s').label('Boundary Probe of Lambda Smooth Wall');
            model.component('smooth').probe('bnd_lambda_s').set('probename', 'lambda_s');
            model.component('smooth').probe('bnd_lambda_s').selection.named('selection_wall_ht_s');
            model.component('smooth').probe('bnd_lambda_s').set('expr', 'mat_smooth_air.def.k(T_s)');
            model.component('smooth').probe('bnd_lambda_s').set('descractive', true);
            model.component('smooth').probe('bnd_lambda_s').set('descr', 'lambda_s');

    model.component('smooth').probe.create('bnd_p_in_s', 'Boundary');

            model.component('smooth').probe('bnd_p_in_s').set('intsurface', true);
            model.component('smooth').probe('bnd_p_in_s').label('Boundary Probe of Pressure Inlet Smooth Wall');
            model.component('smooth').probe('bnd_p_in_s').set('probename', 'p_in_s');
            model.component('smooth').probe('bnd_p_in_s').selection.named('selection_inlet_ht_s');
            model.component('smooth').probe('bnd_p_in_s').set('expr', 'p');
            model.component('smooth').probe('bnd_p_in_s').set('descractive', true);
            model.component('smooth').probe('bnd_p_in_s').set('descr', 'p_in_s');

    model.component('smooth').probe.create('bnd_p_out_s', 'Boundary');

            model.component('smooth').probe('bnd_p_out_s').set('intsurface', true);
            model.component('smooth').probe('bnd_p_out_s').label('Boundary Probe of Pressure Outlet Smooth Wall');
            model.component('smooth').probe('bnd_p_out_s').set('probename', 'p_out_s');
            model.component('smooth').probe('bnd_p_out_s').selection.named('selection_outlet_ht_s');
            model.component('smooth').probe('bnd_p_out_s').set('expr', 'p');
            model.component('smooth').probe('bnd_p_out_s').set('descractive', true);
            model.component('smooth').probe('bnd_p_out_s').set('descr', 'p_out_s');

    model.component('smooth').probe.create('dom_dyn_p_s', 'Domain');

            model.component('smooth').probe('dom_dyn_p_s').set('intsurface', true);
            model.component('smooth').probe('dom_dyn_p_s').set('intvolume', true);
            model.component('smooth').probe('dom_dyn_p_s').label('Boundary Probe of Dynamic Pressure Smooth Wall');
            model.component('smooth').probe('dom_dyn_p_s').selection.named('selection_ht_s');
            model.component('smooth').probe('dom_dyn_p_s').set('probename', 'dyn_p_s');
            model.component('smooth').probe('dom_dyn_p_s').set('expr', '0.5 * spf.U^2 * spf.rho');
            model.component('smooth').probe('dom_dyn_p_s').set('descractive', true);
            model.component('smooth').probe('dom_dyn_p_s').set('descr', 'dyn_p_s');

    model.component('smooth').probe.create('zeta_s', 'GlobalVariable');

            model.component('smooth').probe('zeta_s').set('expr', 'abs(p_out_s - p_in_s) / dyn_p_s');
            model.component('smooth').probe('zeta_s').label('Global Variable Probe of Zeta Smooth Wall');
            model.component('smooth').probe('zeta_s').set('descractive', true);
            model.component('smooth').probe('zeta_s').set('descr', 'zeta_s');

    model.component('smooth').probe.create('Nu_s', 'GlobalVariable');

        model.component('smooth').probe('Nu_s').set('expr', 'alpha_s * l_ref_s / lambda_s');
        model.component('smooth').probe('Nu_s').label('Global Variable Probe of Nu Smooth Wall');
        model.component('smooth').probe('Nu_s').set('descractive', true);
        model.component('smooth').probe('Nu_s').set('descr', 'Nu_s');
    
    % SET MESH

    model.component('smooth').mesh.create('mesh_smooth_normal');

        model.component('smooth').mesh('mesh_smooth_normal').automatic(false);
        model.component('smooth').mesh('mesh_smooth_normal').feature('size').set('custom', true);
        model.component('smooth').mesh('mesh_smooth_normal').feature('size1').selection.named('selection_wall_s');
        model.component('smooth').mesh('mesh_smooth_normal').feature('size1').set('custom', true);
        model.component('smooth').mesh('mesh_smooth_normal').feature('size1').set('hminactive', true);
        model.component('smooth').mesh('mesh_smooth_normal').feature('size1').set('hmaxactive', true);
        model.component('smooth').mesh('mesh_smooth_normal').feature('size1').set('hgradactive', true);
        model.component('smooth').mesh('mesh_smooth_normal').feature('size1').set('hnarrowactive', true);
        model.component('smooth').mesh('mesh_smooth_normal').feature('size1').set('hcurveactive', true);
        model.component('smooth').mesh('mesh_smooth_normal').feature('cr1').selection.named('selection_fluid_s');
        model.component('smooth').mesh('mesh_smooth_normal').feature('cr1').selection('boundary').named('selection_wall_s');
        model.component('smooth').mesh('mesh_smooth_normal').feature('bl1').selection.named('selection_fluid_s');
        model.component('smooth').mesh('mesh_smooth_normal').feature('bl1').feature('blp1').selection.named('selection_wall_s');

        model.component('smooth').mesh('mesh_smooth_normal').feature('size').set('hmax', 0.2);
        model.component('smooth').mesh('mesh_smooth_normal').feature('size').set('hmin', 0.01);
        model.component('smooth').mesh('mesh_smooth_normal').feature('size1').set('hmax', 0.15);
        model.component('smooth').mesh('mesh_smooth_normal').feature('size1').set('hmin', 0.001);

        model.component('smooth').mesh('mesh_smooth_normal').label('Mesh Smooth Normal');

    model.component('smooth').mesh.create('mesh_smooth_mapped');

        model.component('smooth').mesh('mesh_smooth_mapped').feature('size').set('custom', true);
        model.component('smooth').mesh('mesh_smooth_mapped').feature('size').set('hmin', 0.01);
        model.component('smooth').mesh('mesh_smooth_mapped').feature('size').set('hmax', 1);
        
        model.component('smooth').mesh('mesh_smooth_mapped').create('map', 'Map');
        
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').selection.geom('geom_smooth', 2);
        
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').selection.named('selection_fluid_s');
        
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').create('distr_r', 'Distribution');
        
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_r').selection.named('selection_distr_r');
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_r').set('type', 'predefined');
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_r').set('elemcount', 100);
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_r').set('elemratio', 100);
        
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').create('distr_z_ht', 'Distribution');
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_z_ht').selection.named('selection_distr_z_ht');
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_z_ht').set('type', 'predefined');
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_z_ht').set('elemcount', 100);
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_z_ht').set('elemratio', 100);
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_z_ht').set('symmetric', true);
        
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').create('distr_z_in', 'Distribution');
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_z_in').selection.named('selection_distr_z_in');
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_z_in').set('type', 'predefined');
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_z_in').set('elemcount', 100);
        
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').create('distr_z_out', 'Distribution');
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_z_out').selection.named('selection_disrt_z_out');
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_z_out').set('type', 'predefined');
            model.component('smooth').mesh('mesh_smooth_mapped').feature('map').feature('distr_z_out').set('elemcount', 100);
        
        model.component('smooth').mesh('mesh_smooth_mapped').label('Mesh Smooth Mapped');

end


function model = build_component_rough (model, parameters)

    model.component.create('rough');

    model.component('rough').label('Rough Wall');

    % BUILD GEOMETRY

    model.component('rough').geom.create('geom_rough', 2);

        model.component('rough').geom('geom_rough').axisymmetric(true);

        model.component('rough').geom('geom_rough').lengthUnit('mm');

        model.component('rough').geom('geom_rough').axisymmetric(true);

        model.component('rough').geom('geom_rough').lengthUnit('mm');

        model.component('rough').geom('geom_rough').create('r1', 'Rectangle');

            model.component('rough').geom('geom_rough').feature('r1').set('size', {'d_ht / 2' 'l_st_in'});
    
        model.component('rough').geom('geom_rough').create('r2', 'Rectangle');
    
            model.component('rough').geom('geom_rough').feature('r2').set('size', {'d_ht / 2 - h_r' 'w_r'});
            model.component('rough').geom('geom_rough').feature('r2').set('pos', {'0' 'l_st_in'});
    
        model.component('rough').geom('geom_rough').create('r3', 'Rectangle');
    
            model.component('rough').geom('geom_rough').feature('r3').set('size', {'h_r' 'w_r'});
            model.component('rough').geom('geom_rough').feature('r3').set('pos', {'d_ht / 2 - h_r' 'l_st_in'});
    
        model.component('rough').geom('geom_rough').create('r4', 'Rectangle');
    
            model.component('rough').geom('geom_rough').feature('r4').set('size', {'d_ht / 2' 's_r'});
            model.component('rough').geom('geom_rough').feature('r4').set('pos', {'0' 'l_st_in + w_r'});
    
        model.component('rough').geom('geom_rough').create('arr1', 'Array');
    
            model.component('rough').geom('geom_rough').feature('arr1').selection('input').set({'r2' 'r3' 'r4'});
            model.component('rough').geom('geom_rough').feature('arr1').set('displ', {'0' 'w_r + s_r'});
            model.component('rough').geom('geom_rough').feature('arr1').set('fullsize', {'1' 'n'});
    
        model.component('rough').geom('geom_rough').create('r5', 'Rectangle');
    
            model.component('rough').geom('geom_rough').feature('r5').set('size', {'d_ht / 2' 'l_st_out'});
            model.component('rough').geom('geom_rough').feature('r5').set('pos', {'0' 'l_st_in + n * (w_r + s_r)'});
    
        model.component('rough').geom('geom_rough').run;

    % APPOINTMETN OF SELECTION NAMES

    model.component('rough').selection.create('selection_task_r', 'Box');

        model.component('rough').selection('selection_task_r').label('Task Domain');

    model.component('rough').selection.create('selection_in_st_r', 'Box');

        model.component('rough').selection('selection_in_st_r').set('xmin', 0);
        model.component('rough').selection('selection_in_st_r').set('xmax', 'd_ht / 8');
        model.component('rough').selection('selection_in_st_r').set('ymin', '0.1 * l_st_in');
        model.component('rough').selection('selection_in_st_r').set('ymax', '0.9 * l_st_in');

        model.component('rough').selection('selection_in_st_r').label('Input Stabilization Domain');

    model.component('rough').selection.create('selection_ht_r', 'Box');

        model.component('rough').selection('selection_ht_r').set('xmin', 0);
        model.component('rough').selection('selection_ht_r').set('xmax', 'd_ht / 8');
        model.component('rough').selection('selection_ht_r').set('ymin', 'l_st_in + 0.5 * w_r');
        model.component('rough').selection('selection_ht_r').set('ymax', 'l_st_in + (n - 0.5) * (s_r + w_r)');

        model.component('rough').selection('selection_ht_r').label('Heat Transfer Domain');

    model.component('rough').selection.create('selection_out_st_r', 'Box');

        model.component('rough').selection('selection_out_st_r').set('xmin', 0);
        model.component('rough').selection('selection_out_st_r').set('xmax', 'd_ht / 8');
        model.component('rough').selection('selection_out_st_r').set('ymin', 'l_st_in + n * (s_r + w_r) + 0.1 * l_st_out');
        model.component('rough').selection('selection_out_st_r').set('ymax', 'l_st_in + n * (s_r + w_r) + 0.9 * l_st_out');

        model.component('rough').selection('selection_out_st_r').label('Output Stabilization Domain');

    model.component('rough').selection.create('selection_fluid_r', 'Union');

        model.component('rough').selection('selection_fluid_r').set('input', {'selection_in_st_r' 'selection_ht_r' 'selection_out_st_r'});
        model.component('rough').selection('selection_fluid_r').label('Fluid');

    model.component('rough').selection.create('selection_ribs_domain', 'Difference');

        model.component('rough').selection('selection_ribs_domain').set('add', {'selection_task_r'});
        model.component('rough').selection('selection_ribs_domain').set('subtract', {'selection_fluid_r'});

        model.component('rough').selection('selection_ribs_domain').label('Ribs Domain');

    model.component('rough').selection.create('selection_ht_boundary', 'Adjacent');

        model.component('rough').selection('selection_ht_boundary').set('entitydim', 2);
        model.component('rough').selection('selection_ht_boundary').set('input', {'selection_ht_r'});

        model.component('rough').selection('selection_ht_boundary').label('Heat Transfer Boundary');

    model.component('rough').selection.create('selecion_ribs_boundary', 'Adjacent');

        model.component('rough').selection('selecion_ribs_boundary').set('entitydim', 2);
        model.component('rough').selection('selecion_ribs_boundary').set('input', {'selection_ribs_domain'});

        model.component('rough').selection('selecion_ribs_boundary').label('Ribs Boundary');

    model.component('rough').selection.create('selection_inlet_r', 'Box');

        model.component('rough').selection('selection_inlet_r').set('xmin', 'd_ht / 8');
        model.component('rough').selection('selection_inlet_r').set('xmax', 'd_ht / 4');
        model.component('rough').selection('selection_inlet_r').set('ymin', '0');
        model.component('rough').selection('selection_inlet_r').set('ymax', '0');

        model.component('rough').selection('selection_inlet_r').label('Inlet');

        model.component('rough').selection('selection_inlet_r').set('entitydim', 1);

    model.component('rough').selection.create('selection_outlet_r', 'Box');

        model.component('rough').selection('selection_outlet_r').set('xmin', 'd_ht / 8');
        model.component('rough').selection('selection_outlet_r').set('xmax', 'd_ht / 4');
        model.component('rough').selection('selection_outlet_r').set('ymin', 'l_st_in + n * (s_r + w_r) + l_st_out');
        model.component('rough').selection('selection_outlet_r').set('ymax', 'l_st_in + n * (s_r + w_r) + l_st_out');

        model.component('rough').selection('selection_outlet_r').label('Outlet');

        model.component('rough').selection('selection_outlet_r').set('entitydim', 1);

    model.component('rough').selection.create('selection_inlet_ht_r', 'Box');

        model.component('rough').selection('selection_inlet_ht_r').set('entitydim', 1);
        
        model.component('rough').selection('selection_inlet_ht_r').set('xmin', 'd_ht / 16');
        model.component('rough').selection('selection_inlet_ht_r').set('xmax', 'd_ht / 8');
        model.component('rough').selection('selection_inlet_ht_r').set('ymin', 'l_st_in');
        model.component('rough').selection('selection_inlet_ht_r').set('ymax', 'l_st_in');

        model.component('rough').selection('selection_inlet_ht_r').label('Inlet in Heat Transfer Domain');

    model.component('rough').selection.create('selection_outlet_ht_r', 'Box');

        model.component('rough').selection('selection_outlet_ht_r').set('entitydim', 1);

        model.component('rough').selection('selection_outlet_ht_r').set('xmin', 'd_ht / 16');
        model.component('rough').selection('selection_outlet_ht_r').set('xmax', 'd_ht / 8');
        model.component('rough').selection('selection_outlet_ht_r').set('ymax', 'l_st_in + n * (s_r + w_r)');
        model.component('rough').selection('selection_outlet_ht_r').set('ymin', 'l_st_in + n * (s_r + w_r)');

        model.component('rough').selection('selection_outlet_ht_r').label('Outlet in Heat Transfer Domain');
        
    model.component('rough').selection.create('selection_inter_task_1', 'Box');

        model.component('rough').selection('selection_inter_task_1').set('xmin', '0');
        model.component('rough').selection('selection_inter_task_1').set('xmax', '0');
        model.component('rough').selection('selection_inter_task_1').set('ymin', '0.1 * l_st_in');
        model.component('rough').selection('selection_inter_task_1').set('ymax', '0.9 * (l_st_in + n * (s_r + w_r) + l_st_out)');

        model.component('rough').selection('selection_inter_task_1').label('Intermediate Bounadry of Task Domian 1');

        model.component('rough').selection('selection_inter_task_1').set('entitydim', 1);

    model.component('rough').selection.create('selection_inter_task_2', 'Box');

        model.component('rough').selection('selection_inter_task_2').set('xmin', 'd_ht / 16');
        model.component('rough').selection('selection_inter_task_2').set('xmax', 'd_ht / 8');
        model.component('rough').selection('selection_inter_task_2').set('ymin', '0.1 * l_st_in');
        model.component('rough').selection('selection_inter_task_2').set('ymax', '0.9 * (l_st_in + n * (s_r + w_r) + l_st_out)');

        model.component('rough').selection('selection_inter_task_2').label('Intermediate Bounadry of Task Domian 2');

        model.component('rough').selection('selection_inter_task_2').set('entitydim', 1);

    model.component('rough').selection.create('selection_axis_r', 'Difference');

        model.component('rough').selection('selection_axis_r').set('entitydim', 1);       
        model.component('rough').selection('selection_axis_r').set('add', {'selection_inter_task_1'});
        model.component('rough').selection('selection_axis_r').set('subtract', {'selection_inter_task_2'});

        model.component('rough').selection('selection_axis_r').label('Axis');

    model.component('rough').selection.create('selection_wall_ht_1', 'Box');

        model.component('rough').selection('selection_wall_ht_1').set('entitydim', 1);

        model.component('rough').selection('selection_wall_ht_1').set('xmin', 'd_ht / 2 - h_r');
        model.component('rough').selection('selection_wall_ht_1').set('xmax', 'd_ht / 2 + h_r');
        model.component('rough').selection('selection_wall_ht_1').set('ymin', 'l_st_in + 0.5 * w_r');
        model.component('rough').selection('selection_wall_ht_1').set('ymax', 'l_st_in + (n - 0.5) * (s_r + w_r)');

        model.component('rough').selection('selection_wall_ht_1').label('Intermediate Bounadry of Heat Transfer 1');

    model.component('rough').selection.create('selection_wall_ht_2', 'Difference');

        model.component('rough').selection('selection_wall_ht_2').set('entitydim', 1);
        
        model.component('rough').selection('selection_wall_ht_2').set('add', {'selection_wall_ht_1'});
        model.component('rough').selection('selection_wall_ht_2').set('subtract', {'selection_ht_boundary' 'selection_inter_task_2'});

        model.component('rough').selection('selection_wall_ht_2').label('Intermediate Boundary of Heat Transfer 2');

    model.component('rough').selection.create('selection_wall_ht_3', 'Box');

        model.component('rough').selection('selection_wall_ht_3').set('entitydim', 1);

        model.component('rough').selection('selection_wall_ht_3').set('xmin', 0);
        model.component('rough').selection('selection_wall_ht_3').set('xmax', 'd_ht / 2 - h_r / 2');
        model.component('rough').selection('selection_wall_ht_3').set('ymin', 'l_st_in');
        model.component('rough').selection('selection_wall_ht_3').set('ymax', 'l_st_in + n * (s_r + w_r)');

        model.component('rough').selection('selection_wall_ht_3').label('Intermediate Boundary of Heat Transfer 3');

    model.component('rough').selection.create('selection_wall_ht_4', 'Difference');

        model.component('rough').selection('selection_wall_ht_4').set('entitydim', 1);

        model.component('rough').selection('selection_wall_ht_4').set('add', {'selection_ht_boundary'});
        model.component('rough').selection('selection_wall_ht_4').set('subtract', {'selection_wall_ht_3'});

        model.component('rough').selection('selection_wall_ht_4').label('Intermediate Boundary of Heat Transfer 4');

    model.component('rough').selection.create('selection_wall_st_in_r', 'Box');

        model.component('rough').selection('selection_wall_st_in_r').set('entitydim', 1);

        model.component('rough').selection('selection_wall_st_in_r').set('xmin', 'd_ht / 2');
        model.component('rough').selection('selection_wall_st_in_r').set('xmax', 'd_ht / 2');
        model.component('rough').selection('selection_wall_st_in_r').set('ymin', '0.1 * l_st_in');
        model.component('rough').selection('selection_wall_st_in_r').set('ymax', '0.9 * l_st_in');

        model.component('rough').selection('selection_wall_st_in_r').label('Wall of Input Stabilization Domain');

    model.component('rough').selection.create('selection_wall_st_out_r', 'Box');

        model.component('rough').selection('selection_wall_st_out_r').set('entitydim', 1);

        model.component('rough').selection('selection_wall_st_out_r').set('xmin', 'd_ht / 2');
        model.component('rough').selection('selection_wall_st_out_r').set('xmax', 'd_ht / 2');
        model.component('rough').selection('selection_wall_st_out_r').set('ymin', 'l_st_in + n * (s_r + w_r) + 0.1 * l_st_out');
        model.component('rough').selection('selection_wall_st_out_r').set('ymax', 'l_st_in + n * (s_r + w_r) + 0.9 * l_st_out');

        model.component('rough').selection('selection_wall_st_out_r').label('Wall of Output Stabilization Domain');

        model.component('rough').selection.create('selection_temp_cond_r', 'Union');

            model.component('rough').selection('selection_temp_cond_r').set('entitydim', 1);

            model.component('rough').selection('selection_temp_cond_r').set('input', {'selection_wall_ht_2' 'selection_wall_ht_4' 'selection_wall_st_in_r' 'selection_wall_st_out_r'});

            model.component('rough').selection('selection_temp_cond_r').label('Wall of Temperature Condition');

        model.component('rough').selection.create('selection_steam_ribs_bound', 'Difference');

            model.component('rough').selection('selection_steam_ribs_bound').set('entitydim', 1);

            model.component('rough').selection('selection_steam_ribs_bound').set('add', {'selecion_ribs_boundary'});
            model.component('rough').selection('selection_steam_ribs_bound').set('subtract', {'selection_temp_cond_r'});

            model.component('rough').selection('selection_steam_ribs_bound').label('Boundary of Streamlined Ribs');

        model.component('rough').selection.create('selection_wall_fluid_ht', 'Union');

            model.component('rough').selection('selection_wall_fluid_ht').set('entitydim', 1);
            
            model.component('rough').selection('selection_wall_fluid_ht').set('input', {'selection_steam_ribs_bound', 'selection_wall_ht_4'});
            
            model.component('rough').selection('selection_wall_fluid_ht').label('Wall Fluid in Heat Transfer Domain');
    
        model.component('rough').selection.create('selection_wall_fluid', 'Union');
    
            model.component('rough').selection('selection_wall_fluid').set('entitydim', 1);

            model.component('rough').selection('selection_wall_fluid').set('input', {'selection_wall_st_in_r', 'selection_wall_fluid_ht', 'selection_wall_st_out_r'});
            
            model.component('rough').selection('selection_wall_fluid').label('Wall of Fluid');

    % SET OF MATERIAL PROPORTIES

        % AIR PROPORTIES

        model.component('rough').material.create('mat_rough_air', 'Common');

            model.component('rough').material('mat_rough_air').propertyGroup('def').func.create('eta', 'Piecewise');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func.create('Cp', 'Piecewise');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func.create('rho', 'Analytic');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func.create('k', 'Piecewise');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func.create('cs', 'Analytic');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func.create('an1', 'Analytic');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func.create('an2', 'Analytic');
            model.component('rough').material('mat_rough_air').propertyGroup.create('RefractiveIndex', 'Refractive index');
            model.component('rough').material('mat_rough_air').propertyGroup.create('NonlinearModel', 'Nonlinear model');
            model.component('rough').material('mat_rough_air').propertyGroup.create('idealGas', 'Ideal gas');
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').func.create('Cp', 'Piecewise');
            model.component('rough').material('mat_rough_air').label('Air');
            model.component('rough').material('mat_rough_air').set('family', 'air');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('eta').set('arg', 'T');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('eta').set('argunit', 'K');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('Cp').set('arg', 'T');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('Cp').set('argunit', 'K');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/R_const[K*mol/J]/T');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('rho').set('argunit', 'Pa,K');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('rho').set('plotargs', {'pA' '0' '1'; 'T' '0' '1'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('k').set('arg', 'T');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('k').set('argunit', 'K');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*R_const[K*mol/J]/0.02897*T)');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('cs').set('args', {'T'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('cs').set('argunit', 'K');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('cs').set('fununit', 'm/s');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('cs').set('plotargs', {'T' '273.15' '373.15'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an1').set('expr', '-1/rho(pA,T)*d(rho(pA,T),T)');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an1').set('args', {'pA' 'T'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an1').set('argunit', 'Pa,K');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an1').set('fununit', '1/K');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an1').set('plotargs', {'pA' '101325' '101325'; 'T' '273.15' '373.15'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an2').set('funcname', 'muB');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an2').set('expr', '0.6*eta(T)');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an2').set('args', {'T'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an2').set('argunit', 'K');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an2').set('fununit', 'Pa*s');
            model.component('rough').material('mat_rough_air').propertyGroup('def').func('an2').set('plotargs', {'T' '200' '1600'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('thermalexpansioncoefficient', '');
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('molarmass', '');
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('bulkviscosity', '');
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('molarmass', '0.02897[kg/mol]');
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('bulkviscosity', 'muB(T)');
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('ratioofspecificheat', '1.4');
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('electricconductivity', {'0[S/m]' '0' '0' '0' '0[S/m]' '0' '0' '0' '0[S/m]'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('heatcapacity', 'Cp(T)');
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('density', 'rho(pA,T)');
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('thermalconductivity', {'k(T)' '0' '0' '0' 'k(T)' '0' '0' '0' 'k(T)'});
            model.component('rough').material('mat_rough_air').propertyGroup('def').set('soundspeed', 'cs(T)');
            model.component('rough').material('mat_rough_air').propertyGroup('def').addInput('temperature');
            model.component('rough').material('mat_rough_air').propertyGroup('def').addInput('pressure');
            model.component('rough').material('mat_rough_air').propertyGroup('RefractiveIndex').set('n', '');
            model.component('rough').material('mat_rough_air').propertyGroup('RefractiveIndex').set('ki', '');
            model.component('rough').material('mat_rough_air').propertyGroup('RefractiveIndex').set('n', '');
            model.component('rough').material('mat_rough_air').propertyGroup('RefractiveIndex').set('ki', '');
            model.component('rough').material('mat_rough_air').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
            model.component('rough').material('mat_rough_air').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
            model.component('rough').material('mat_rough_air').propertyGroup('NonlinearModel').set('BA', '(def.gamma+1)/2');
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').func('Cp').label('Piecewise 2');
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').func('Cp').set('arg', 'T');
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').func('Cp').set('argunit', 'K');
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').func('Cp').set('fununit', 'J/(kg*K)');
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').set('Rs', 'R_const/Mn');
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').set('heatcapacity', 'Cp(T)');
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').set('ratioofspecificheat', '1.4');
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').set('molarmass', '0.02897');
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').addInput('temperature');
            model.component('rough').material('mat_rough_air').propertyGroup('idealGas').addInput('pressure');
            model.component('rough').material('mat_rough_air').set('groups', {});
            model.component('rough').material('mat_rough_air').set('family', 'air');

        % STEEL PROPORTIES

        model.component('rough').material.create('mat_rough_steel', 'Common');

            model.component('rough').material('mat_rough_steel').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
            model.component('rough').material('mat_rough_steel').label('Steel AISI 4340');
            model.component('rough').material('mat_rough_steel').set('family', 'steel');
            model.component('rough').material('mat_rough_steel').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
            model.component('rough').material('mat_rough_steel').propertyGroup('def').set('electricconductivity', {'4.032e6[S/m]' '0' '0' '0' '4.032e6[S/m]' '0' '0' '0' '4.032e6[S/m]'});
            model.component('rough').material('mat_rough_steel').propertyGroup('def').set('thermalexpansioncoefficient', {'12.3e-6[1/K]' '0' '0' '0' '12.3e-6[1/K]' '0' '0' '0' '12.3e-6[1/K]'});
            model.component('rough').material('mat_rough_steel').propertyGroup('def').set('heatcapacity', '475[J/(kg*K)]');
            model.component('rough').material('mat_rough_steel').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
            model.component('rough').material('mat_rough_steel').propertyGroup('def').set('density', '7850[kg/m^3]');
            model.component('rough').material('mat_rough_steel').propertyGroup('def').set('thermalconductivity', {'44.5[W/(m*K)]' '0' '0' '0' '44.5[W/(m*K)]' '0' '0' '0' '44.5[W/(m*K)]'});
            model.component('rough').material('mat_rough_steel').propertyGroup('Enu').set('youngsmodulus', '205e9[Pa]');
            model.component('rough').material('mat_rough_steel').propertyGroup('Enu').set('poissonsratio', '0.28');
            model.component('rough').material('mat_rough_steel').set('groups', {});
            model.component('rough').material('mat_rough_steel').set('family', 'steel');
            model.component('rough').material('mat_rough_air').selection.named('selection_fluid_r');
            model.component('rough').material('mat_rough_steel').selection.named('selection_ribs_domain');

    % SET MULTIPHYSICS
            
        % SET SPF MODULE
        
        model.component('rough').physics.create('spf_rough', 'TurbulentFlowkeps', 'geom_rough');

            model.component('rough').physics('spf_rough').selection.named('selection_fluid_r');

            model.component('rough').physics('spf_rough').feature('init1').set('u_init', {'0' '0' 'v_z_in'});

            model.component('rough').physics('spf_rough').create('inl1', 'InletBoundary', 1);

                    model.component('rough').physics('spf_rough').feature('inl1').selection.named('selection_inlet_r');
                    model.component('rough').physics('spf_rough').feature('inl1').set('U0in', 'v_z_in');

            model.component('rough').physics('spf_rough').create('out1', 'OutletBoundary', 1);

                    model.component('rough').physics('spf_rough').feature('out1').set('p0', 'p_0');
                    model.component('rough').physics('spf_rough').feature('out1').selection.named('selection_outlet_r');

        % SET HT MODULE

        model.component('rough').physics.create('ht_rough', 'HeatTransferInSolidsAndFluids', 'geom_rough');

            model.component('rough').physics('ht_rough').selection.named('selection_task_r');

            model.component('rough').physics('ht_rough').feature('fluid1').selection.named('selection_fluid_r');

            model.component('rough').physics('ht_rough').feature('init1').set('Tinit', 'T_0');

            model.component('rough').physics('ht_rough').create('temp1', 'TemperatureBoundary', 1);

                    model.component('rough').physics('ht_rough').feature('temp1').selection.named('selection_temp_cond_r');
                    model.component('rough').physics('ht_rough').feature('temp1').set('T0', 'T_wall');

            model.component('rough').physics('ht_rough').create('ifl1', 'Inflow', 1);

                    model.component('rough').physics('ht_rough').feature('ifl1').selection.named('selection_inlet_r');
                    model.component('rough').physics('ht_rough').feature('ifl1').set('Tustr', 'T_0');

            model.component('rough').physics('ht_rough').create('ofl1', 'ConvectiveOutflow', 1);

                    model.component('rough').physics('ht_rough').feature('ofl1').selection.named('selection_outlet_r');

        model.component('rough').multiphysics.create('nitf_rough', 'NonIsothermalFlow', 2);        

    % SET PROBE VARIABLES

    model.component('rough').probe.create('dom_mr_r', 'Domain');

        model.component('rough').probe('dom_mr_r').set('intsurface', true);
        model.component('rough').probe('dom_mr_r').set('intvolume', true);
        model.component('rough').probe('dom_mr_r').label('Domain Probe of Mass Rate');
        model.component('rough').probe('dom_mr_r').selection.named('selection_ht_r');
        model.component('rough').probe('dom_mr_r').set('expr', 'spf2.rho * spf2.U');
        model.component('rough').probe('dom_mr_r').set('probename', 'mr_r');

    model.component('rough').probe.create('bnd_T_r', 'Domain');

        model.component('rough').probe('bnd_T_r').set('intsurface', true);
        model.component('rough').probe('bnd_T_r').set('intvolume', true);
        model.component('rough').probe('bnd_T_r').set('expr', 'spf2.rho * spf2.U * T2');
        model.component('rough').probe('bnd_T_r').selection.named('selection_ht_r');
        model.component('rough').probe('bnd_T_r').set('probename', 'T_r_int');
        model.component('rough').probe('bnd_T_r').label('Domain Probe of Integral Temperature');

    model.component('rough').probe.create('var_T_r', 'GlobalVariable');

        model.component('rough').probe('var_T_r').label('Global Variable Probe of Mean Temperature');
        model.component('rough').probe('var_T_r').set('probename', 'T_r');
        model.component('rough').probe('var_T_r').set('expr', 'T_r_int / mr_r');
        model.component('rough').probe('var_T_r').set('descractive', true);
        model.component('rough').probe('var_T_r').set('descr', 'T_r');

    model.component('rough').probe.create('bnd_alpha_r', 'Boundary');

        model.component('rough').probe('bnd_alpha_r').set('intsurface', true);
        model.component('rough').probe('bnd_alpha_r').label('Boundary Probe of Alpha Rough Wall');
        model.component('rough').probe('bnd_alpha_r').set('probename', 'alpha_r');
        model.component('rough').probe('bnd_alpha_r').selection.named('selection_wall_fluid_ht');
        model.component('rough').probe('bnd_alpha_r').set('expr', 'abs(ht2.ntflux/(T_r - T_wall))');
        model.component('rough').probe('bnd_alpha_r').set('descractive', true);
        model.component('rough').probe('bnd_alpha_r').set('descr', 'alpha_r');

    model.component('rough').probe.create('bnd_l_ref_r', 'Boundary');

        model.component('rough').probe('bnd_l_ref_r').set('intsurface', true);
        model.component('rough').probe('bnd_l_ref_r').label('Boundary Probe of Reference Length Rough Wall');
        model.component('rough').probe('bnd_l_ref_r').set('probename', 'l_ref_r');
        model.component('rough').probe('bnd_l_ref_r').selection.named('selection_wall_fluid_ht');
        model.component('rough').probe('bnd_l_ref_r').set('expr', '1 / (2 * (pi * d_ht + l_ht))');
        model.component('rough').probe('bnd_l_ref_r').set('type', 'integral');
        model.component('rough').probe('bnd_l_ref_r').set('unit', 'm');
        model.component('rough').probe('bnd_l_ref_r').set('descractive', true);
        model.component('rough').probe('bnd_l_ref_r').set('descr', 'l_ref_r');

    model.component('rough').probe.create('bnd_lambda_r', 'Boundary');

        model.component('rough').probe('bnd_lambda_r').set('intsurface', true);
        model.component('rough').probe('bnd_lambda_r').label('Boundary Probe of Lambda Rough Wall');
        model.component('rough').probe('bnd_lambda_r').set('probename', 'lambda_r');
        model.component('rough').probe('bnd_lambda_r').selection.named('selection_wall_fluid_ht');
        model.component('rough').probe('bnd_lambda_r').set('expr', 'mat_rough_air.def.k(T_r)');
        model.component('rough').probe('bnd_lambda_r').set('descractive', true);
        model.component('rough').probe('bnd_lambda_r').set('descr', 'lambda_r');

    model.component('rough').probe.create('bnd_p_in_r', 'Boundary');

        model.component('rough').probe('bnd_p_in_r').set('intsurface', true);
        model.component('rough').probe('bnd_p_in_r').label('Boundary Probe of Pressure Inlet Rough Wall');
        model.component('rough').probe('bnd_p_in_r').set('probename', 'p_in_r');
        model.component('rough').probe('bnd_p_in_r').selection.named('selection_inlet_ht_r');
        model.component('rough').probe('bnd_p_in_r').set('expr', 'p2');
        model.component('rough').probe('bnd_p_in_r').set('descractive', true);
        model.component('rough').probe('bnd_p_in_r').set('descr', 'p_in_r');

    model.component('rough').probe.create('bnd_p_out_r', 'Boundary');

        model.component('rough').probe('bnd_p_out_r').set('intsurface', true);
        model.component('rough').probe('bnd_p_out_r').label('Boundary Probe of Pressure Outlet Rough Wall');
        model.component('rough').probe('bnd_p_out_r').set('probename', 'p_out_r');
        model.component('rough').probe('bnd_p_out_r').selection.named('selection_outlet_ht_r');
        model.component('rough').probe('bnd_p_out_r').set('expr', 'p2');
        model.component('rough').probe('bnd_p_out_r').set('descractive', true);
        model.component('rough').probe('bnd_p_out_r').set('descr', 'p_out_r');

    model.component('rough').probe.create('dom_dyn_p_r', 'Domain');

        model.component('rough').probe('dom_dyn_p_r').set('intsurface', true);
        model.component('rough').probe('dom_dyn_p_r').set('intvolume', true);
        model.component('rough').probe('dom_dyn_p_r').label('Boundary Probe of Dynamic Pressure Rough Wall');
        model.component('rough').probe('dom_dyn_p_r').selection.named('selection_ht_r');
        model.component('rough').probe('dom_dyn_p_r').set('probename', 'dyn_p_r');
        model.component('rough').probe('dom_dyn_p_r').set('expr', '0.5 * spf2.U^2 * spf2.rho');
        model.component('rough').probe('dom_dyn_p_r').set('descractive', true);
        model.component('rough').probe('dom_dyn_p_r').set('descr', 'dyn_p_r');

    model.component('rough').probe.create('zeta_r', 'GlobalVariable');

        model.component('rough').probe('zeta_r').set('expr', 'abs(p_out_r - p_in_r) / dyn_p_r');
        model.component('rough').probe('zeta_r').label('Global Variable Probe of Zeta Rough Wall');
        model.component('rough').probe('zeta_r').set('descractive', true);
        model.component('rough').probe('zeta_r').set('descr', 'zeta_r');

    model.component('rough').probe.create('Nu_r', 'GlobalVariable');

        model.component('rough').probe('Nu_r').set('expr', 'alpha_r * l_ref_r / lambda_r');
        model.component('rough').probe('Nu_r').label('Global Variable Probe of Nu Rough Wall');
        model.component('rough').probe('Nu_r').set('descractive', true);
        model.component('rough').probe('Nu_r').set('descr', 'Nu_r');

    model.component('rough').probe.create('J', 'GlobalVariable');

        model.component('rough').probe('J').set('expr', '(rough.Nu_r / smooth.Nu_s) ^1.4 / (rough.zeta_r / smooth.zeta_s) ^0.4');
        model.component('rough').probe('J').label('Global Variable Probe of Criterion Rough Wall');
        model.component('rough').probe('J').set('descractive', true);
        model.component('rough').probe('J').set('descr', 'J');

    % SET MESH

    model.component('rough').mesh.create('mesh_rough_normal');

        model.component('rough').mesh('mesh_rough_normal').automatic(false);
        model.component('rough').mesh('mesh_rough_normal').feature('size').set('custom', true);
        model.component('rough').mesh('mesh_rough_normal').feature('size1').selection.named('selection_fluid_r');
        model.component('rough').mesh('mesh_rough_normal').feature('size2').selection.named('selection_wall_fluid');
        model.component('rough').mesh('mesh_rough_normal').feature('size2').set('custom', true);
        model.component('rough').mesh('mesh_rough_normal').feature('size2').set('hmaxactive', true);
        model.component('rough').mesh('mesh_rough_normal').feature('size2').set('hminactive', true);
        model.component('rough').mesh('mesh_rough_normal').feature('size2').set('hgradactive', true);
        model.component('rough').mesh('mesh_rough_normal').feature('size2').set('hcurveactive', true);
        model.component('rough').mesh('mesh_rough_normal').feature('size2').set('hnarrowactive', true);
        model.component('rough').mesh('mesh_rough_normal').feature('cr1').selection.named('selection_task_r');
        model.component('rough').mesh('mesh_rough_normal').feature('cr1').selection('boundary').named('selection_wall_fluid');
        model.component('rough').mesh('mesh_rough_normal').feature('bl1').selection.named('selection_fluid_r');
        model.component('rough').mesh('mesh_rough_normal').feature('bl1').feature('blp1').selection.named('selection_wall_fluid');
        model.component('rough').mesh('mesh_rough_normal').feature('size1').selection.named('selection_task_r');
        model.component('rough').mesh('mesh_rough_normal').feature('size1').set('custom', true);
        model.component('rough').mesh('mesh_rough_normal').feature('size1').set('hmaxactive', true);
        model.component('rough').mesh('mesh_rough_normal').feature('size1').set('hminactive', true);
        model.component('rough').mesh('mesh_rough_normal').feature('size1').set('hgradactive', true);
        model.component('rough').mesh('mesh_rough_normal').feature('size1').set('hcurveactive', true);
        model.component('rough').mesh('mesh_rough_normal').feature('size1').set('hnarrowactive', true);

        model.component('rough').mesh('mesh_rough_normal').feature('size').set('hmax', 15);
        model.component('rough').mesh('mesh_rough_normal').feature('size').set('hmin', 0.1);
        model.component('rough').mesh('mesh_rough_normal').feature('size1').set('hmax', 0.2);
        model.component('rough').mesh('mesh_rough_normal').feature('size1').set('hmin', 0.01);
        model.component('rough').mesh('mesh_rough_normal').feature('size2').set('hmax', 0.15);
        model.component('rough').mesh('mesh_rough_normal').feature('size2').set('hmin', 0.001);

        model.component('rough').mesh('mesh_rough_normal').label('Mesh Rough Normal');

end

function model = set_studies (model)

    % SET PARAMETRIC SWEEP WITH OPTIMIZATION TASK

    model.study.create('std_sweep_opt_normal');

        model.study('std_sweep_opt_normal').setGenPlots(false);

        model.study('std_sweep_opt_normal').create('stat', 'Stationary');

            model.study('std_sweep_opt_normal').label('Study of Combine Sweep and Optimization Normal');

            model.study('std_sweep_opt_normal').feature('stat').activate('spf_smooth', true);
            model.study('std_sweep_opt_normal').feature('stat').activate('ht_smooth', true);
            model.study('std_sweep_opt_normal').feature('stat').activate('spf_rough', true);
            model.study('std_sweep_opt_normal').feature('stat').activate('ht_rough', true);
            model.study('std_sweep_opt_normal').feature('stat').activate('nitf_smooth', true);
            model.study('std_sweep_opt_normal').feature('stat').activate('nitf_rough', true);

            model.study('std_sweep_opt_normal').feature('stat').setIndex('mesh', 'mesh_smooth_mapped', 1);
            model.study('std_sweep_opt_normal').feature('stat').setIndex('mesh', 'mesh_rough_normal', 3);

        model.study('std_sweep_opt_normal').create('param_Re', 'Parametric');

            model.study('std_sweep_opt_normal').feature('param_Re').setIndex('pname', 'Re', 0);
            model.study('std_sweep_opt_normal').feature('param_Re').setIndex('plistarr', 'range(Re_min, Re_step, Re_max)', 0);
            model.study('std_sweep_opt_normal').feature('param_Re').setIndex('punit', 'mm', 0);
            model.study('std_sweep_opt_normal').feature('param_Re').label('Parametric Sweep Re');

            model.study('std_sweep_opt_normal').feature('param_Re').set('accumtableall', false);
            model.study('std_sweep_opt_normal').feature('param_Re').set('useaccumtable', false);

        model.study('std_sweep_opt_normal').create('opt', 'Optimization');

            model.study('std_sweep_opt_normal').feature('opt').setIndex('optobj', 'rough.J', 0);
            model.study('std_sweep_opt_normal').feature('opt').set('objectivetype', 'maximization');
        
            model.study('std_sweep_opt_normal').feature('opt').setIndex('pname', 'h_r', 0);
            model.study('std_sweep_opt_normal').feature('opt').setIndex('pname', 's_r', 1);

            model.study('std_sweep_opt_normal').feature('opt').setIndex('initval', '0.5[mm]', 0);
            model.study('std_sweep_opt_normal').feature('opt').setIndex('initval', '1[mm]', 1);
        
            model.study('std_sweep_opt_normal').feature('opt').setIndex('lbound', 'h_r_min', 0);
            model.study('std_sweep_opt_normal').feature('opt').setIndex('ubound', 'h_r_max', 0);
            model.study('std_sweep_opt_normal').feature('opt').setIndex('lbound', 's_r_min', 1);
            model.study('std_sweep_opt_normal').feature('opt').setIndex('ubound', 's_r_max', 1);
        
            model.study('std_sweep_opt_normal').feature('opt').set('probesel', 'none');
            model.study('std_sweep_opt_normal').feature('opt').set('showindobj', false);
            model.study('std_sweep_opt_normal').feature('opt').set('useobjtable', true);
            model.study('std_sweep_opt_normal').feature('opt').set('plotobj', true);
            model.study('std_sweep_opt_normal').feature('opt').set('window', 'new');
            model.study('std_sweep_opt_normal').feature('opt').set('useconstrtable', false);

end