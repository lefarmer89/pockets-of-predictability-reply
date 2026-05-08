function figure2()
%FIGURE2  Render Figure 2: correlations between excess-return forecasts
% and Coibion-Gorodnichenko forecast errors, for the nine main
% forecasting series and three SPF target variables (gy, ue, ip). Bars
% show 95% HAC confidence intervals.
%
% Reads:  data/quarterly_avg_forecast_revsanderrs.csv
%         results/daily/forecastResults_1_2.5yE_1yDM_1S_pm.mat
% Writes: gcf (saved to .eps by saveFigureScript).

folder = fullfile('results', 'daily');

% Load the forecast-error panel and forecasted excess returns.
data1 = readtable(fullfile('data', 'quarterly_avg_forecast_revsanderrs.csv'));
data2 = load(fullfile(folder, 'forecastResults_1_2.5yE_1yDM_1S_pm.mat'), ...
    'yF2MatS', 'dateVecS');

varNames = categorical(cellstr(data1.forvar));
data_3tb = data1(varNames == '3tb', :);
data_aaa = data1(varNames == 'aaa', :);
data_gy  = data1(varNames == 'gy',  :);
data_ue  = data1(varNames == 'ue',  :);
data_c   = data1(varNames == 'c',   :);
data_ip  = data1(varNames == 'ip',  :);

quarterlyYears    = data_3tb.yyyy;
quarterlyQuarters = data_3tb.qq;

% yF2MatS columns: 1-7 = dp tbl tsp rvar mv pc erL; 8-10 = comb1/2/3.
% Figure 2 omits erL.
erForecasts = data2.yF2MatS(:, [1:6, 8:10]);
erForecasts(any(month(data2.dateVecS) == [2 3 5 6 8 9 11 12], 2), :) = NaN;
dailyDates  = data2.dateVecS;

cgFEDaily = NaN(size(dailyDates, 1), 6);
for t = 1:numel(quarterlyYears)
    dateInd = (year(dailyDates) == quarterlyYears(t)) & ...
              (quarter(dailyDates) == quarterlyQuarters(t));
    cgFEDaily(dateInd, :) = repmat([data_3tb.ind_fe(t), data_aaa.ind_fe(t), ...
        data_gy.ind_fe(t), data_ue.ind_fe(t), data_c.ind_fe(t), data_ip.ind_fe(t)], ...
        [sum(dateInd), 1]);
end

% Correlations and HAC standard errors.
corrMatFE   = NaN(9, 6);
corrSEMatFE = NaN(9, 6);
for ii = 1:9
    for jj = 1:6
        keepInd = ~isnan(erForecasts(:, ii)) & ~isnan(cgFEDaily(:, jj));
        corrMatFE(ii, jj) = corr(erForecasts(:, ii), cgFEDaily(:, jj), ...
            'rows', 'complete');
        temp = regstats2Fast(erForecasts(:, ii), cgFEDaily(:, jj), 'linear', ...
            {'beta', 'hac'});
        corrSEMatFE(ii, jj) = temp.hac.se(2) * ...
            std(cgFEDaily(keepInd, jj), 'omitnan') / ...
            std(erForecasts(keepInd, ii), 'omitnan');
    end
end

select1 = [3 4 6];   % gy, ue, ip
select2 = 1:9;

corrMatFEUB = corrMatFE + 1.96 * corrSEMatFE;
corrMatFELB = corrMatFE - 1.96 * corrSEMatFE;

figure(2)
b = bar(corrMatFE(select2, select1)', 'grouped');
hold on
[ngroups, nbars] = size(corrMatFE(select2, select1)');
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i, :) = b(i).XEndPoints;
end
errorbar(x', corrMatFE(select2, select1)', ...
    corrMatFEUB(select2, select1)' - corrMatFELB(select2, select1)', ...
    'k', 'linestyle', 'none')
hold off
set(gca, 'xticklabel', {'gy', 'ue', 'ip'})
xlabel('Coibion-Gorodnichenko Variables')
ylabel('Correlation')
legend('dp', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3', ...
    'Orientation', 'horizontal', 'Location', 'best')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
end
