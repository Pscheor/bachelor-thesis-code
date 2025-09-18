function [dprime, c, nHit, nMiss, nCR, nFA] = computeSDT(codes)
    codes(codes == 4) = [];

    % Zählwerte
    nHit  = sum(codes == 0);
    nMiss = sum(codes == 1);
    nFA   = sum(codes == 2);
    nCR   = sum(codes == 3);

    nSignal = nHit + nMiss;
    nNoise  = nFA  + nCR;

    % Raten berechnen
    if nSignal > 0
        hitRate = min(max(nHit / nSignal, 1/(2*nSignal)), 1 - 1/(2*nSignal));
    else
        hitRate = NaN;
    end

    if nNoise > 0
        faRate = min(max(nFA / nNoise, 1/(2*nNoise)), 1 - 1/(2*nNoise));
    else
        faRate = NaN;
    end

    % SDT-Metriken berechnen
    if ~isnan(hitRate) && ~isnan(faRate)
        dprime = norminv(hitRate) - norminv(faRate);
        c = -0.5 * (norminv(hitRate) + norminv(faRate));
    else
        dprime = NaN;
        c = NaN;
    end
end

%%Option A: So lassen
% Die Korrektur (siehe zeile 15) schützt vor extremen Verzerrungen (z. B. unendlichem d′)
% 
% Standard in vielen SDT-Implementierungen (Macmillan & Creelman-Logik)