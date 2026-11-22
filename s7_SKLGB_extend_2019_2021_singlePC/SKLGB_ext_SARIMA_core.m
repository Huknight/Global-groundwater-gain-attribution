% SOURCE-PRESERVING EXTRACTION from SARIMA_global_forecast.m
% Original per-grid estimator, validation search, forecast calls,
% NRMSE/TDA and failure handling are unchanged.
% The wrapper supplies testYears=[2019 2020 2021] in place of
% [2022 2023 2024], with identical original time-index formula.
% This function is called with 2000 selected master rows only.
function [pred36, valPred36, nrmse3, tda3, bestP3, valBest3, ...
    rmse3, testRange3, success3, failureCode3, failureID3] = ...
    SKLGB_ext_SARIMA_core(Y, p_list, d, q, S, P, D_seasonal, ...
    Q, testYears, horizon)

    Y = double(Y(:));
    nRound = numel(testYears);
    pred36 = NaN(nRound * horizon, 1);
    valPred36 = NaN(nRound * horizon, 1);
    nrmse3 = NaN(nRound, 1);
    tda3 = NaN(nRound, 1);
    bestP3 = NaN(nRound, 1);
    valBest3 = NaN(nRound, 1);
    rmse3 = NaN(nRound, 1);
    testRange3 = NaN(nRound, 1);
    success3 = false(nRound, 1);
    failureCode3 = uint8(zeros(nRound, 1));
    failureID3 = cell(nRound, 1);

    for r = 1:nRound
        testYear = testYears(r);
        trainEndIdx = (testYear - 2005) * 12;
        valIdx = (trainEndIdx + 1):(trainEndIdx + horizon);
        testIdx = (trainEndIdx + horizon + 1): ...
            (trainEndIdx + 2 * horizon);
        predIdx = ((r - 1) * horizon + 1):(r * horizon);

        Y_train = Y(1:trainEndIdx);
        Y_val = Y(valIdx);
        Y_test = Y(testIdx);
        bestValNRMSE = Inf;
        bestModel = [];
        bestValPred = [];
        bestP = NaN;
        firstValidationError = '';

        for k = 1:numel(p_list)
            p = p_list(k);
            maxLag = max([p, q, S * P, S * Q]);
            if maxLag >= numel(Y_train)
                continue;
            end
            try
                model = makeSARIMAModel(p, d, q, S, P, D_seasonal, Q);
                estModel = estimate(model, Y_train, 'Display', 'off');
                Y_val_pred = forecastSARIMANumeric(estModel, Y_train, horizon);
                valScore = calcNRMSE(Y_val, Y_val_pred);
                if isfinite(valScore) && valScore < bestValNRMSE
                    bestValNRMSE = valScore;
                    bestModel = estModel;
                    bestValPred = Y_val_pred;
                    bestP = p;
                end
            catch ME
                if isempty(firstValidationError)
                    firstValidationError = ME.identifier;
                end
            end
        end

        if isempty(bestModel)
            failureCode3(r) = 1;
            failureID3{r} = firstValidationError;
            continue;
        end

        bestP3(r) = bestP;
        valBest3(r) = bestValNRMSE;
        valPred36(predIdx) = bestValPred;

        % Use the observed validation year only as forecast history.
        Y_train_val = Y(1:trainEndIdx + horizon);
        try
            Y_test_pred = forecastSARIMANumeric(bestModel, Y_train_val, horizon);
        catch ME
            failureCode3(r) = 2;
            failureID3{r} = ME.identifier;
            continue;
        end

        if any(~isfinite(Y_test_pred))
            failureCode3(r) = 3;
            continue;
        end

        pred36(predIdx) = Y_test_pred;
        testNRMSE = calcNRMSE(Y_test, Y_test_pred);
        if ~isfinite(testNRMSE)
            failureCode3(r) = 4;
            continue;
        end

        nrmse3(r) = testNRMSE;
        rmse3(r) = sqrt(mean((Y_test_pred - Y_test).^2));
        testRange3(r) = max(Y_test) - min(Y_test);
        tda3(r) = calcTDA(Y_test, Y_test_pred);
        success3(r) = true;
    end
end

function model = makeSARIMAModel(p, d, q, S, P, D_seasonal, Q)
    if D_seasonal == 1
        model = arima('ARLags', 1:p, 'D', d, 'MALags', 1:q, ...
            'Seasonality', S, 'SARLags', S * (1:P), ...
            'SMALags', S * (1:Q));
    else
        model = arima('ARLags', 1:p, 'D', d, 'MALags', 1:q, ...
            'SARLags', S * (1:P), 'SMALags', S * (1:Q));
    end
end

function yhat = forecastSARIMANumeric(model, historyY, horizon)
    yhat = forecast(model, horizon, 'Y0', historyY(:));
    yhat = double(yhat(:));
    if numel(yhat) ~= horizon
        error('SKLGB:ForecastLength', ...
            'Forecast returned %d months; expected %d.', numel(yhat), horizon);
    end
end

function nrmse = calcNRMSE(obs, pred)
    obs = double(obs(:));
    pred = double(pred(:));
    if isempty(obs) || numel(obs) ~= numel(pred) || ...
            any(~isfinite(obs)) || any(~isfinite(pred))
        nrmse = NaN;
        return;
    end
    rmse = sqrt(mean((pred - obs).^2));
    denom = max(obs) - min(obs);
    if denom == 0
        if rmse == 0
            nrmse = 0;
        else
            nrmse = Inf;
        end
    else
        nrmse = rmse / denom;
    end
end

function tda = calcTDA(obs, pred)
    obs = double(obs(:));
    pred = double(pred(:));
    if numel(obs) ~= numel(pred) || numel(obs) < 2 || ...
            any(~isfinite(obs)) || any(~isfinite(pred))
        tda = NaN;
        return;
    end
    actualDiff = diff(obs);
    predDiff = diff(pred);
    tda = mean(sign(actualDiff) == sign(predDiff)) * 100;
end
