function [] = plotDGMs(dgms)
    %:param dgms: A cell array of persistence diagrams returned
    % from ripser
    concat_dgms = [];
    for ii = 1:length(dgms)
        concat_dgms = [concat_dgms; dgms{ii}(:)];
    end
    has_inf = sum(isinf(concat_dgms)) > 0;
    finite_dgms = concat_dgms(isfinite(concat_dgms));
    
    ax_min = min(finite_dgms);
    ax_max = max(finite_dgms);
    x_r = ax_max - ax_min;
    buffer = x_r/5;
    x_down = ax_min - buffer/2;
    x_up = ax_max + buffer;
    y_down = x_down;
    y_up = x_up;
    yr = y_up - y_down;
    

    hold on;
    
    if has_inf
        b_inf = y_down + yr*0.95;
        for ii = 1:length(dgms)
            dgms{ii}(isinf(dgms{ii})) = b_inf;
        end
    end
    
    labels = {};
    for ii = 1:length(dgms)
        scatter(dgms{ii}(:, 1), dgms{ii}(:, 2), 'fill');
        labels{end+1} = sprintf('H%i', ii-1);
    end
    
    % Plot diagonal and horizontal infinity line if applicable
    plot([x_down, x_up], [x_down, x_up], '--');
    labels{end+1} = 'diag';
    if has_inf
        plot([x_down, x_up], [b_inf, b_inf], '--');
        labels{end+1} = '\infty';
    end
    legend(labels);
    xlim([x_down, x_up]);
    ylim([y_down, y_up]);
    axis equal;
end

