import com.comsol.model.*
import com.comsol.model.util.*

model = mphopen('calculate_project.mph');

model = plotting (model);

mphsave(model, strcat(pwd, '\postproccesing_project.mph'));

function model = plotting (model)

    model.result('pg1').label('Probe Plot Parametric Sweep');
    model.result('pg2').label('Optimization Task');

    model.result.table('tbl1').label('Probe Table Parametric Sweep');
    model.result.table('tbl3').label('Objective Table');

    model.view.create('view_plot', 2);

        model.view('view_plot').label('View Sweep Grid');
        model.view('view_plot').axis.set('xmin', 's_r_min*1e3');
        model.view('view_plot').axis.set('xmax', 's_r_max*1e3');
        model.view('view_plot').axis.set('ymin', 'h_r_min*1e3');
        model.view('view_plot').axis.set('ymax', 'h_r_max*1e3');

    model.result.numerical('gev3').set('data', 'dset4');

        model.result.table.create('tbl_zeta_r', 'Table');
    
            model.result.table('tbl_zeta_r').comments('Global Variable Probe of Zeta Rough Wall');
            model.result.numerical('gev3').set('table', 'tbl_zeta_r');
            model.result.numerical('gev3').setResult;
            model.result.table('tbl_zeta_r').label('Table of Zeta Rough Wall');
    
        model.result.export.create('tbl_zeta_r_export', 'tbl_zeta_r', 'Table');
    
            model.result.export('tbl_zeta_r_export').label('Table of Zeta Rough Wall');
            model.result.export('tbl_zeta_r_export').set('filename', 'param_sweep_zeta_r_normal.txt');
            model.result.export('tbl_zeta_r_export').run;
    
        model.result.create('pg_zeta_r', 'PlotGroup2D');
    
            model.result('pg_zeta_r').label('Parametric Sweep Zeta Rough Wall');
    
            model.result('pg_zeta_r').set('data', 'none');
    
            model.result('pg_zeta_r').create('tbls_zeta_r', 'TableSurface');
    
                model.result('pg_zeta_r').feature('tbls_zeta_r').label('Table Surface');
    
                model.result('pg_zeta_r').feature('tbls_zeta_r').set('source', 'table');
                model.result('pg_zeta_r').feature('tbls_zeta_r').set('table', 'tbl_zeta_r');   
    
                model.result('pg_zeta_r').feature('tbls_zeta_r').set('plotdata', 'manual');
                model.result('pg_zeta_r').feature('tbls_zeta_r').set('xaxisdata', 's_r');

                model.result('pg_zeta_r').set('titletype', 'manual');
                model.result('pg_zeta_r').set('title', 'Table Surface of Zeta_r=Zeta_r(h_r, s_r)');
                
                model.result('pg_zeta_r').feature('tbls_zeta_r').set('colortable', 'HeatCamera');    
    
            model.result('pg_zeta_r').set('xlabelactive', true);
            model.result('pg_zeta_r').set('ylabelactive', true);
            model.result('pg_zeta_r').set('xlabel', 's_r [mm]');
            model.result('pg_zeta_r').set('ylabel', 'h_r [mm]');
            
            model.result('pg_zeta_r').set('view', 'view_plot');
            model.result('pg_zeta_r').set('showlegendsmaxmin', true);
            model.result('pg_zeta_r').set('showlegendsunit', true);
    
            model.result('pg_zeta_r').run;

    model.result.numerical('gev4').set('data', 'dset4');

    model.result.table.create('tbl_Nu_r', 'Table');

        model.result.table('tbl_Nu_r').comments('Global Variable Probe of Nu Rough Wall');
        model.result.numerical('gev4').set('table', 'tbl_Nu_r');
        model.result.numerical('gev4').setResult;
        model.result.table('tbl_Nu_r').label('Table of Nu Rough Wall');

    model.result.export.create('tbl_Nu_r_export', 'tbl_Nu_r', 'Table');

        model.result.export('tbl_Nu_r_export').label('Table of Nu Roudh Wall');
        model.result.export('tbl_Nu_r_export').set('filename', 'param_sweep_Nu_r_normal.txt');
        model.result.export('tbl_Nu_r_export').run;

    model.result.create('pg_Nu_r', 'PlotGroup2D');

        model.result('pg_Nu_r').label('Parametric Sweep Nu Rough Wall');

        model.result('pg_Nu_r').set('data', 'none');

        model.result('pg_Nu_r').create('tbls_Nu_r', 'TableSurface');

            model.result('pg_Nu_r').feature('tbls_Nu_r').label('Table Surface');

            model.result('pg_Nu_r').feature('tbls_Nu_r').set('source', 'table');
            model.result('pg_Nu_r').feature('tbls_Nu_r').set('table', 'tbl_Nu_r');   

            model.result('pg_Nu_r').feature('tbls_Nu_r').set('plotdata', 'manual');
            model.result('pg_Nu_r').feature('tbls_Nu_r').set('xaxisdata', 's_r');

            model.result('pg_Nu_r').set('titletype', 'manual');
            model.result('pg_Nu_r').set('title', 'Table Surface of Nu_r=Nu_r(h_r, s_r)');
            
            model.result('pg_Nu_r').feature('tbls_Nu_r').set('colortable', 'HeatCamera');    

        model.result('pg_Nu_r').set('xlabelactive', true);
        model.result('pg_Nu_r').set('ylabelactive', true);
        model.result('pg_Nu_r').set('xlabel', 's_r [mm]');
        model.result('pg_Nu_r').set('ylabel', 'h_r [mm]');
        
        model.result('pg_Nu_r').set('view', 'view_plot');
        model.result('pg_Nu_r').set('showlegendsmaxmin', true);
        model.result('pg_Nu_r').set('showlegendsunit', true);

        model.result('pg_Nu_r').run;

    model.result.numerical('gev5').set('data', 'dset3');

    model.result.table.create('tbl_criterion', 'Table');

        model.result.table('tbl_criterion').comments('Global Variable Probe of Criterion');
        model.result.numerical('gev5').set('table', 'tbl_criterion');
        model.result.numerical('gev5').setResult;
        model.result.table('tbl_criterion').label('Table of Criterion');

    model.result.export.create('tbl_criterion_export', 'tbl_criterion', 'Table');

        model.result.export('tbl_criterion_export').label('Table of Criterion');
        model.result.export('tbl_criterion_export').set('filename', 'param_sweep_J_normal.txt');
        model.result.export('tbl_criterion_export').run;

    model.result.create('pg_criterion', 'PlotGroup2D');

        model.result('pg_criterion').label('Parametric Sweep Criterion');

        model.result('pg_criterion').set('data', 'none');

        model.result('pg_criterion').create('tbls_J', 'TableSurface');

            model.result('pg_criterion').feature('tbls_J').label('Table Surface');

            model.result('pg_criterion').feature('tbls_J').set('source', 'table');
            model.result('pg_criterion').feature('tbls_J').set('table', 'tbl_criterion');   

            model.result('pg_criterion').feature('tbls_J').set('plotdata', 'manual');
            model.result('pg_criterion').feature('tbls_J').set('xaxisdata', 's_r');

            model.result('pg_criterion').set('titletype', 'manual');
            model.result('pg_criterion').set('title', 'Table Surface of J=J(h_r, s_r)');
            
            model.result('pg_criterion').feature('tbls_J').set('colortable', 'HeatCamera');    

        model.result('pg_criterion').set('xlabelactive', true);
        model.result('pg_criterion').set('ylabelactive', true);
        model.result('pg_criterion').set('xlabel', 's_r [mm]');
        model.result('pg_criterion').set('ylabel', 'h_r [mm]');
        
        model.result('pg_criterion').set('view', 'view_plot');
        model.result('pg_criterion').set('showlegendsmaxmin', true);
        model.result('pg_criterion').set('showlegendsunit', true);

        model.result('pg_criterion').run;

        model.component('rough').view.create('view_mesh', 'geom_rough');

            model.component('rough').view('view_mesh').label('View Mesh');
            model.component('rough').view('view_mesh').axis.set('xmin', '-h_r');
            model.component('rough').view('view_mesh').axis.set('xmax', 'd_ht / 2 + h_r');
            model.component('rough').view('view_mesh').axis.set('ymin', 'l_in_st - s_r / 2');
            model.component('rough').view('view_mesh').axis.set('ymax', 'l_in_st + s_r * 2');

        model.result.create('pg_mesh_smooth_normal', 'PlotGroup2D');

            model.result('pg_mesh_smooth_normal').set('data', 'none');

            model.result('pg_mesh_smooth_normal').set('titletype', 'manual');
            model.result('pg_mesh_smooth_normal').set('title', 'Mesh of Smooth Wall');

            model.result('pg_mesh_smooth_normal').create('mesh_smooth', 'Mesh');

                model.result('pg_mesh_smooth_normal').feature('mesh_smooth').set('data', 'dset1');

                model.result('pg_mesh_smooth_normal').feature('mesh_smooth').set('colortablerev', false);
                model.result('pg_mesh_smooth_normal').feature('mesh_smooth').set('colorlegend', false);

                model.result('pg_mesh_smooth_normal').feature('mesh_smooth').set('elemcolor', 'gray');
                model.result('pg_mesh_smooth_normal').feature('mesh_smooth').set('wireframecolor', 'black');
                model.result('pg_mesh_smooth_normal').feature('mesh_smooth').set('resolution', 'norefine');                

                model.result('pg_mesh_smooth_normal').set('view', 'view_mesh');
                model.result('pg_mesh_smooth_normal').label('Mesh Smooth Wall Normal');

                model.result('pg_mesh_smooth_normal').run;

        model.result.create('pg_mesh_rough_normal', 'PlotGroup2D');

            model.result('pg_mesh_rough_normal').set('data', 'none');

            model.result('pg_mesh_rough_normal').set('titletype', 'manual');
            model.result('pg_mesh_rough_normal').set('title', 'Mesh of Rough Wall');

            model.result('pg_mesh_rough_normal').create('mesh_rough', 'Mesh');

                model.result('pg_mesh_rough_normal').feature('mesh_rough').set('data', 'dset2');

                model.result('pg_mesh_rough_normal').feature('mesh_rough').set('colortablerev', false);
                model.result('pg_mesh_rough_normal').feature('mesh_rough').set('colorlegend', false);

                model.result('pg_mesh_rough_normal').feature('mesh_rough').set('elemcolor', 'gray');
                model.result('pg_mesh_rough_normal').feature('mesh_rough').set('wireframecolor', 'black');
                model.result('pg_mesh_rough_normal').feature('mesh_rough').set('resolution', 'norefine');                

                model.result('pg_mesh_rough_normal').set('view', 'view_mesh');
                model.result('pg_mesh_rough_normal').label('Mesh Rough Wall Normal');

                model.result('pg_mesh_rough_normal').run;

end