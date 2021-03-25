%%% 
%%% calcRF.m
%%%
%%% Convenience function to calculate RFs.
%%%
function [G,eps] = calcRF (X,Y,Nrf)

  %%% Total number of samples
  Nt = length(X);
  
  %%% Construct matrix of lagged predictor variable X
  S = zeros(Nt-Nrf+1,Nrf);
  for p = 1:Nt-Nrf+1
    for q = 1:Nrf
      S(p,q) = X(p-q+Nrf);
    end
  end
  
  %%% Predictand - only latter part of time series is used because we need
  %%% at least Nrf samples of X prior to each point in Y
  P = Y(Nrf:Nt);
  
  %%% Least squares regression to determine the Green's function (the RF)
  G = S \ P;
  
  %%% Error in reconstructed time series
  eps = P - S*G;
  
end

