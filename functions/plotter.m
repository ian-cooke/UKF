function plotter(figs, tspan, type, data, filter, filter_2, dataset, dataset_2)
%PLOTTER General plotter for Project 1, basically a giany switch statement
scale = 1000; % for conversion to km
cov_bound = 3;
nan_color = 'g-';
switch type
    %%%%%%%%%%
    % state  %
    %%%%%%%%%%
    case 'state'
        % 1) plot spacecraft state
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        ylabels = {'$x$ [km]', '$y$ [km]', '$z$ [km]', '$\dot{x}$ [km/s]', '$\dot{y}$ [km/s]', '$\dot{z}$ [km/s]'};
        for i = 1:6
            ax = subplot(2,3,i);
            xlabel('Time [hr]')
            ylabel(ylabels{i}, 'Interpreter', 'latex')
            hold on
            grid on
            plot(tspan./3600, dataset.xhat(i,:)./scale, 'Linewidth', 1.5)
            plot(tspan./3600, dataset.xhat_nan(i, :)./scale, 'Linewidth', 1.5)
            plot(tspan./3600, dataset.xhat(i,:)./scale +cov_bound*sqrt(reshape(dataset.Phat(i, i, :), 1, []))./scale, 'r--', 'Linewidth', 1)
            plot(tspan./3600, dataset.xhat(i,:)./scale -cov_bound*sqrt(reshape(dataset.Phat(i, i, :), 1, []))./scale, 'r--', 'Linewidth', 1)
            linker(i) = ax;
        end
        linkaxes(linker, 'x')
        mtit(figs(1), sprintf('Spacecraft State w/ 3-sigma Bounds Estimated by %s', filter));
        linker = [];
        
        %%%%%%%%%%%%%%
        % state_nan  %
        %%%%%%%%%%%%%%
    case 'state_nan'
        % 1) plot spacecraft dstate
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        ylabels = {'$ x$ [km]', '$ y$ [km]', '$ z$ [km]', '$\dot{x}$ [km/s]', '$\dot{y}$ [km/s]', '$\dot{z}$ [km/s]'};
        for i = 1:6
            ax = subplot(2,3,i);
            xlabel('Time [hr]')
            ylabel(ylabels{i}, 'Interpreter', 'latex')
            hold on
            grid on
            plot(tspan./3600, dataset.xhat_nan(i,:)./scale, 'Linewidth', 1.5)
            plot(tspan./3600, dataset.xhat_nan(i,:)./scale +cov_bound*sqrt(reshape(dataset.Phat_nan(i, i, :), 1, []))./scale, 'r--', 'Linewidth', 1)
            plot(tspan./3600, dataset.xhat_nan(i,:)./scale -cov_bound*sqrt(reshape(dataset.Phat_nan(i, i, :), 1, []))./scale, 'r--', 'Linewidth', 1)
            linker(i) = ax;
        end
        linkaxes(linker, 'x')
        mtit(figs(1), sprintf('Spacecraft State w/ 3-sigma Bounds Estimated by %s', filter));
        linker = [];
        
        
        %%%%%%%%%%
        % dstate %
        %%%%%%%%%%
    case 'dstate'
        % 1) plot spacecraft dstate
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        ylabels = {'$\delta x$ [km]', '$\delta y$ [km]', '$\delta z$ [km]', '$\delta\dot{x}$ [km/s]', '$\delta\dot{y}$ [km/s]', '$\delta\dot{z}$ [km/s]'};
        for i = 1:6
            ax = subplot(2,3,i);
            xlabel('Time [hr]')
            ylabel(ylabels{i}, 'Interpreter', 'latex')
            hold on
            grid on
            plot(tspan./3600, dataset.dxhat(i,:)./scale, 'Linewidth', 1.5)
            plot(tspan./3600, dataset.dxhat_nan(i,:)./scale, nan_color, 'Linewidth', 1.5)
            plot(tspan./3600, cov_bound*sqrt(reshape(dataset.Phat(i, i, :), 1, []))./scale, 'r--', 'Linewidth', 1)
            plot(tspan./3600, -cov_bound*sqrt(reshape(dataset.Phat(i, i, :), 1, []))./scale, 'r--', 'Linewidth', 1)
            linker(i) = ax;
        end
        linkaxes(linker, 'x')
        mtit(figs(1), sprintf('Spacecraft State Error w/ 3-sigma Bounds Estimated by %s', filter));
        linker = [];
        
        
        %%%%%%%%%%%%%%%
        % dstate_true %
        %%%%%%%%%%%%%%%
    case 'dstate_true'
        % 1) plot spacecraft dstate_true
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        ylabels = {'$\Delta x$ [km]', '$\Delta y$ [km]', '$\Delta z$ [km]', '$\Delta\dot{x}$ [km/s]', '$\Delta\dot{y}$ [km/s]', '$\Delta\dot{z}$ [km/s]'};
        for i = 1:6
            ax = subplot(2,3,i);
            xlabel('Time [hr]')
            ylabel(ylabels{i}, 'Interpreter', 'latex')
            hold on
            grid on
            plot(tspan./3600, dataset.deltax(i,:)./scale, 'Linewidth', 1.5)
            plot(tspan./3600, dataset.deltax_nan(i,:)./scale, nan_color, 'Linewidth', 1.5)
            linker(i) = ax;
        end
        linkaxes(linker, 'x')
        mtit(figs(1), sprintf('Spacecraft State Error Relative to xhat Estimated by %s', filter));
        linker = [];
        
        %%%%%%%%%%%%%%
        % dstate_nan %
        %%%%%%%%%%%%%%
    case 'dstate_nan'
        % 1) plot spacecraft dstate
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        ylabels = {'$\delta x$ [km]', '$\delta y$ [km]', '$\delta z$ [km]', '$\delta\dot{x}$ [km/s]', '$\delta\dot{y}$ [km/s]', '$\delta\dot{z}$ [km/s]'};
        for i = 1:6
            ax = subplot(2,3,i);
            xlabel('Time [hr]')
            ylabel(ylabels{i}, 'Interpreter', 'latex')
            hold on
            grid on
            plot(tspan./3600, dataset.dxhat_nan(i,:)./scale, 'Linewidth', 1.5)
            plot(tspan./3600, cov_bound*sqrt(reshape(dataset.Phat_nan(i, i, :), 1, []))./scale, 'r--', 'Linewidth', 1)
            plot(tspan./3600, -cov_bound*sqrt(reshape(dataset.Phat_nan(i, i, :), 1, []))./scale, 'r--', 'Linewidth', 1)
            linker(i) = ax;
        end
        linkaxes(linker, 'x')
        mtit(figs(1), sprintf('Spacecraft State Error w/ 3-sigma Bounds Estimated by %s', filter));
        linker = [];
        
        %%%%%%%%%
        % resid %
        %%%%%%%%%
    case 'resid'
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        switch data
            case 'both'
                ylabels = {'$\delta\rho$ [km]', '$\delta\dot{\rho}$ [km/s]'};
                for i = 1:2
                    ax = subplot(2,1,i);
                    xlabel('Time [hr]')
                    ylabel(ylabels{i}, 'Interpreter', 'latex')
                    hold on
                    grid on
                    if i == 1
                        title(sprintf('Range Residual Estimated by %s', filter))
                        plot(tspan./3600, dataset.resid(1, :)./scale, 'Linewidth', 1.5)
                        plot(tspan./3600, dataset.resid(3, :)./scale, 'Linewidth', 1.5)
                        plot(tspan./3600, dataset.resid(5, :)./scale, 'Linewidth', 1.5)
                    else
                        title(sprintf('Range Rate Residual Estimated by %s', filter))
                        plot(tspan./3600, dataset.resid(2, :)./scale, 'Linewidth', 1.5)
                        plot(tspan./3600, dataset.resid(4, :)./scale, 'Linewidth', 1.5)
                        plot(tspan./3600, dataset.resid(6, :)./scale, 'Linewidth', 1.5)
                        legend('Station 101', 'Station 337', 'Station 394', 'Location', 'EastOutside')
                    end
                    linker(i) = ax;
                end
                linkaxes(linker, 'x')
                linker = [];
            otherwise
                hold on
                grid on
                xlabel('Time [hr]')
                if strcmp(data, 'rho')
                    title(sprintf('Range Residual Estimated by %s', filter))
                else
                    title(sprintf('Range Rate Residual Estimated by %s', filter))
                end
                plot(tspan./3600, dataset.resid(1, :)./scale, 'Linewidth', 1.5)
                plot(tspan./3600, dataset.resid(2, :)./scale, 'Linewidth', 1.5)
                plot(tspan./3600, dataset.resid(3, :)./scale, 'Linewidth', 1.5)
                legend('Station 101', 'Station 337', 'Station 394', 'Location', 'EastOutside')
        end
        
        
        %%%%%%%%%
        % trace %
        %%%%%%%%%
    case 'trace'
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        ylabels = {'tr$(P_{pos})$ [km$^2$]', 'tr$(P_{vel})$ [km$^2$/s$^2$]'};
        for i = 1:2
            ax = subplot(2,1,i);
            set(ax, 'YScale', 'log')
            xlabel('Time [hr]')
            ylabel(ylabels{i}, 'Interpreter', 'latex')
            grid on
            hold on
            if i == 1
                title(sprintf('S/C Position Covariance Trace Estimated by %s', filter))
                semilogy(tspan./3600, dataset.trace_pos./(scale^2), 'Linewidth', 1.5)
                semilogy(tspan./3600, dataset.trace_pos_nan./(scale^2), nan_color, 'Linewidth', 1.5)
            else
                title(sprintf('S/C Velocity Covariance Trace Estimated by %s', filter))
                semilogy(tspan./3600, dataset.trace_vel./(scale^2), 'Linewidth', 1.5)
                semilogy(tspan./3600, dataset.trace_vel_nan./(scale^2), nan_color, 'Linewidth', 1.5)
            end
            linker(i) = ax;
        end
        linkaxes(linker, 'x')
        
        
        %%%%%%%%%%%%%
        % trace_nan %
        %%%%%%%%%%%%%
    case 'trace_nan'
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        ylabels = {'tr$(P_{pos})$ [km$^2$]', 'tr$(P_{vel})$ [km$^2$/s$^2$]'};
        for i = 1:2
            ax = subplot(2,1,i);
            set(ax, 'YScale', 'log')
            xlabel('Time [hr]')
            ylabel(ylabels{i}, 'Interpreter', 'latex')
            grid on
            hold on
            if i == 1
                title(sprintf('S/C Position Covariance Trace, When Observation Available, Estimated by %s', filter))
                semilogy(tspan./3600, dataset.trace_pos_nan./(scale^2), 'Linewidth', 1.5)
            else
                title(sprintf('S/C Velocity Covariance Trace, When Observation Available, Estimated by %s', filter))
                semilogy(tspan./3600, dataset.trace_vel_nan./(scale^2), 'Linewidth', 1.5)
            end
            linker(i) = ax;
        end
        linkaxes(linker, 'x')
        
        
        %%%%%%%%%
        % cond  %
        %%%%%%%%%
    case 'cond'
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        hold on
        grid on
        xlabel('Time [hr]')
        ylabel('log(cond($\hat{P}^+$))', 'Interpreter', 'latex')
        title('Log Of Posterior Covariance Matrix Condition Number vs. Time')
        plot(tspan, log10(dataset.mat_cond), 'Linewidth', 1.5)
        plot(tspan, log10(dataset.mat_cond_nan), nan_color, 'Linewidth', 1.5)
        plot(tspan, ones(1,length(tspan))*log10(1/eps), 'r--', 'Linewidth', 1.5)
        leg = legend('log(cond($\hat{P}^+$))', 'log($\frac{1}{\epsilon_{machine}}$)', 'Location', 'best');
        set(leg, 'Interpreter', 'latex');
        
        
        %%%%%%%%%%%
        % ellipse %
        %%%%%%%%%%%
    case 'ellipse'
        figure(figs(1))
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        
        subplot(1,2,1);
        hold on
        grid on
        xlabel('s/c position $x$ [km]', 'interpreter', 'latex')
        ylabel('s/c position $y$ [km]', 'interpreter', 'latex')
        zlabel('s/c position $z$ [km]', 'interpreter', 'latex')
        title(sprintf('Spacecraft Position Covariance Ellipsoids Estimated by %s', filter))
        surf(dataset.xyz_ellipse_r.x./scale, dataset.xyz_ellipse_r.y./scale, dataset.xyz_ellipse_r.z./scale);
        %shading interp;
        colormap default;
        axis equal;
        alpha(0.7);
        
        subplot(1,2,2);
        hold on
        grid on
        xlabel('s/c velocity $\dot{x}$ [km]', 'interpreter', 'latex')
        ylabel('s/c position $\dot{y}$ [km]', 'interpreter', 'latex')
        zlabel('s/c position $\dot{z}$ [km]', 'interpreter', 'latex')
        title(sprintf('Spacecraft Position Covariance Ellipsoids Estimated by %s', filter))
        surf(dataset.xyz_ellipse_v.x./scale, dataset.xyz_ellipse_v.y./scale, dataset.xyz_ellipse_v.z./scale);
        %shading interp;
        colormap default;
        axis equal;
        alpha(0.7);
        
        
        %%%%%%%%%%%%%
        % iteration %
        %%%%%%%%%%%%%
    case 'iteration'
        figure(figs(1))
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        
        ax = subplot(1,2,1);
        set(ax, 'YScale', 'log')
        hold on
        grid on
        xlabel('Iteration No.')
        xticks(0:1:3)
        ylabel('tr($P_{0,r}^+$) [km$^2$]', 'interpreter', 'latex')
        title(sprintf('Trace Of Position Covariance vs. Iteration by %s', filter))
        semilogy(0:1:3, dataset.trace_pos_iter./(scale), '-ro', 'Linewidth', 1.5)
        linker(1) = ax;
        
        ax = subplot(1,2,2);
        set(ax, 'YScale', 'log')
        hold on
        grid on
        xlabel('Iteration No.')
        xticks(0:1:3)
        ylabel('tr($P_{0,v}^+$) [km$^2$/s$^2$]', 'interpreter', 'latex')
        title(sprintf('Trace Of Velocity Covariance vs. Iteration by %s', filter))
        semilogy(0:1:3, dataset.trace_vel_iter./(scale^2), '-bo', 'Linewidth', 1.5)
        linker(2) = ax;
        linkaxes(linker, 'x');
        
    %%%%%%%%%%%%%
    % rms_resid %
    %%%%%%%%%%%%%
    case 'rms_resid'
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        ax = subplot(2,1,1);
        set(ax, 'XScale', 'log')
        set(ax, 'Yscale', 'log')
        hold on
        grid on
        xlabel('$\sigma_u$ [m/s$^2$]', 'Interpreter', 'latex')
        ylabel('RMS Value [km]', 'Interpreter', 'latex')
        title('Range Residual RMS vs. $\sigma_u$', 'Interpreter', 'Latex')
        loglog(tspan, dataset(1,:)./scale, 'Linewidth', 1.5)
        linker(1) = ax;
        
        ax = subplot(2,1,2);
        set(ax, 'XScale', 'log')
        set(ax, 'Yscale', 'log')
        hold on
        grid on
        xlabel('$\sigma_u$ [m/s$^2$]', 'Interpreter', 'latex')
        ylabel('RMS Value [km/s]', 'Interpreter', 'latex')
        title('Range Rate Residual RMS vs. $\sigma_u$', 'Interpreter', 'Latex')
        loglog(tspan, dataset(2,:)./scale, 'Linewidth', 1.5)
        linker(2) = ax;
        linkaxes(linker, 'x');
        
    %%%%%%%%%%
    % rms_3d %
    %%%%%%%%%%
    case 'rms_3d'
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        ax = subplot(2,1,1);
        set(ax, 'XScale', 'log')
        set(ax, 'YScale', 'log')
        hold on
        grid on
        xlabel('$\sigma_u \; [m/s^2]$', 'Interpreter', 'latex')
        ylabel('RMS Value [km]', 'Interpreter', 'latex')
        title('S/C Position Error RMS vs. $\sigma_u$', 'Interpreter', 'Latex')
       	loglog(tspan, dataset(1,:)./scale, 'Linewidth', 1.5)
        linker(1) = ax;
        
        ax = subplot(2,1,2);
        set(ax, 'XScale', 'log')
        set(ax, 'Yscale', 'log')
        hold on
        grid on
        xlabel('$\sigma_u \; [m/s^2]$', 'interpreter', 'latex')
        ylabel('RMS Value [km/s]', 'Interpreter', 'latex')
        title('Spacecraft Velocity Error RMS vs. $\sigma_u$', 'Interpreter', 'Latex')
        loglog(tspan, dataset(2,:)./scale, 'Linewidth', 1.5)
        linker(2) = ax;
        linkaxes(linker, 'x');
        
    %%%%%%%%%%%%%%%%%%%%%
    % compare_ckf_other %
    %%%%%%%%%%%%%%%%%%%%%
    case 'compare_ckf_other'
        % dstate
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        ylabels = {'$\delta\hat{x}^+_{x,ckf} - \delta\hat{x}^+_{x,batch}$ [km]', '$\delta\hat{x}^+_{y,ckf} - \delta\hat{x}^+_{y,batch}$ [km]' ,'$\delta\hat{x}^+_{z,ckf} - \delta\hat{x}^+_{z,batch}$ [km]', '$\delta\hat{x}^+_{\dot{x},ckf} - \delta\hat{x}^+_{\dot{x},batch}$ [km/s]', '$\delta\hat{x}^+_{\dot{y},ckf} - \delta\hat{x}^+_{\dot{y},batch}$ [km/s]', '$\delta\hat{x}^+_{\dot{z},ckf} - \delta\hat{x}^+_{\dot{z},batch}$ [km/s]'};
        for i = 1:6
            ax = subplot(2,3,i);
            hold on
            grid on
            xlabel('Time t [hr]')
            ylabel(ylabels{i}, 'Interpreter', 'latex')
            plot(tspan./3600, (dataset.dxhat(i,:) - dataset_2.dxhat(i,:))./scale, 'Linewidth', 1.5)
            plot(tspan./3600, (dataset.dxhat_nan(i,:) - dataset_2.dxhat_nan(i,:))./scale, nan_color, 'Linewidth', 1.5)
            linker(i) = ax;
        end
        linkaxes(linker, 'x')
        mtit(figs(1), sprintf('Difference In State Errors Between %s and %s', filter, filter_2));
        linker = [];
        
        % bounds
        figure(figs(2))
        set(figs(2), 'defaultaxesfontsize', 16)
        set(figs(2), 'Position', get(0, 'Screensize'));
        ylabels = {'$\sqrt{P^+_{x,ckf}} - \sqrt{P^+_{x,batch}}$ [km]', '$\sqrt{P^+_{y,ckf}} - \sqrt{P^+_{y,batch}}$ [km]' ,'$\sqrt{P^+_{z,ckf}} - \sqrt{P^+_{z,batch}}$ [km]',...
            '$\sqrt{P^+_{\dot{x},ckf}} - \sqrt{P^+_{\dot{x},batch}} [km/s]$', '$\sqrt{P^+_{\dot{y},ckf}} - \sqrt{P^+_{\dot{y},batch}}$ [km/s]', '$\sqrt{P^+_{\dot{z},ckf}} - \sqrt{P^+_{\dot{z},batch}}$ [km/s]'};
        for i = 1:6
            ax = subplot(2,3,i);
            hold on
            grid on
            xlabel('Time t [hr]')
            ylabel(ylabels{i}, 'Interpreter', 'latex')
            plot(tspan./3600, (sqrt(reshape(dataset.Phat(i, i, :), 1, [])) - sqrt(reshape(dataset_2.Phat(i, i, :), 1, [])))./scale, 'Linewidth', 1.5)
            plot(tspan./3600, (sqrt(reshape(dataset.Phat_nan(i, i, :), 1, [])) - sqrt(reshape(dataset_2.Phat_nan(i, i, :), 1, [])))./scale, nan_color, 'Linewidth', 1.5)
            linker(i) = ax;
        end
        linkaxes(linker, 'x')
        mtit(figs(1), sprintf('Difference In State Error Bounds Between %s and %s', filter, filter_2));
        linker = [];
        
        % residuals
        figure(figs(3))
        set(figs(3), 'defaultaxesfontsize', 16)
        set(figs(3), 'Position', get(0, 'Screensize'));
        ylabels = {'$\delta\rho_{ckf} - \delta\rho_{batch}$ [km]', '$\delta\dot{\rho}_{ckf} - \delta\dot{\rho}_{batch}$ [km/s]'};
        
        ax = subplot(2,1,1);
        hold on
        grid on
        xlabel('Time t [hr]')
        ylabel(ylabels{1}, 'Interpreter', 'latex')
        plot(tspan./3600, (dataset.resid(1,:) - dataset_2.resid(1,:))./scale, 'Linewidth', 1.5)
        plot(tspan./3600, (dataset.resid(3,:) - dataset_2.resid(3,:))./scale, 'Linewidth', 1.5)
        plot(tspan./3600, (dataset.resid(5,:) - dataset_2.resid(5,:))./scale, 'Linewidth', 1.5)
        linker(1) = ax;
        
        ax = subplot(2,1,2);
        hold on
        grid on
        xlabel('Time t [hr]')
        ylabel(ylabels{2}, 'Interpreter', 'latex')
        plot(tspan./3600, (dataset.resid(2,:) - dataset_2.resid(2,:))./scale, 'Linewidth', 1.5)
        plot(tspan./3600, (dataset.resid(4,:) - dataset_2.resid(4,:))./scale, 'Linewidth', 1.5)
        plot(tspan./3600, (dataset.resid(6,:) - dataset_2.resid(6,:))./scale, 'Linewidth', 1.5)
        legend('Station 1', 'Station 2', 'Station 3', 'Location', 'EastOutside')
        linker(2) = ax;
        linkaxes(linker, 'x')
        mtit(figs(1), sprintf('Difference In Post-Fit Measurement Residuals Between %s and %s', filter, filter_2));
        linker = [];
        
    %%%%%%%%%%%%%%%%
    % compare_cond %
    %%%%%%%%%%%%%%%%
    case 'compare_cond'
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        hold on
        grid on
        xlabel('Time [hr]')
        ylabel('log(cond($\hat{P}^+$))', 'Interpreter', 'latex')
        title('Log Of Posterior Covariance Matrix Condition Number vs. Time')
        plot(tspan, log10(dataset.mat_cond), 'Linewidth', 1.5)
        plot(tspan, log10(dataset_2.mat_cond), 'Linewidth', 1.5)
        plot(tspan, ones(1,length(tspan))*log10(1/eps), 'r--', 'Linewidth', 1.5)
        leg = legend('log(cond($\hat{P}_{CKF}^+$))', 'log(cond($\hat{R}_{SRIF}^+$))', 'log($\frac{1}{\epsilon_{machine}}$)', 'Location', 'best');
        set(leg, 'Interpreter', 'latex');
        
    %%%%%%%%%%%%%%%%%%%
    % compare_ekf_ukf %
    %%%%%%%%%%%%%%%%%%%
    case 'compare_ekf_ukf'
        figure(figs(1))
        set(figs(1), 'defaultaxesfontsize', 16)
        set(figs(1), 'Position', get(0, 'Screensize'));
        ax = subplot(2,1,1);
        hold on
        grid on
        xlabel('$\Delta x_0$ Scale Factor', 'Interpreter', 'latex')
        ylabel('3D RMS, Position [m]', 'Interpreter', 'latex')
       	plot(tspan, dataset(1, :), 'Linewidth', 1.5)
        plot(tspan, dataset(2, :), 'Linewidth', 1.5)
        linker(1) = ax;
        leg = legend('EKF', 'UKF');
        set(leg, 'Location', 'best', 'Interpreter', 'latex');
        
        ax = subplot(2,1,2);
        hold on
        grid on
        xlabel('$\Delta x_0$ Scale Factor', 'Interpreter', 'latex')
        ylabel('3D RMS, Velocity [m/s]', 'Interpreter', 'latex')
       	plot(tspan, dataset_2(1, :), 'Linewidth', 1.5)
        plot(tspan, dataset_2(2, :), 'Linewidth', 1.5)
        linker(1) = ax;
        leg = legend('EKF', 'UKF');
        set(leg, 'Location', 'best', 'Interpreter', 'latex');
        
        linkaxes(linker, 'x')
        mtit(figs(1), sprintf('3D RMS Values for %s and %s verus Initial Perturbation Scale Factor', filter, filter_2));
        linker = [];
        
        
end

end

