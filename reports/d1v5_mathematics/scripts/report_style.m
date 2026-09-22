function S = report_style
%REPORT_STYLE shared colors, fonts, and line widths for every report figure.
% Categorical order and the ordinal blue ramp follow a colorblind-validated
% palette (light print surface).
    S.cat = hex2rgb({'#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300','#4a3aa7','#e34948'});
    S.seq = hex2rgb({'#86b6ef','#5598e7','#2a78d6','#1c5cab','#104281'}); % ordinal: coarse -> fine
    S.ink = hex2rgb({'#0b0b0b'}); S.ink2 = hex2rgb({'#52514e'});
    S.axis = hex2rgb({'#c3c2b7'}); S.grid = hex2rgb({'#e1e0d9'}); S.band = hex2rgb({'#eef3fb'});
    S.fs = 9; S.lw = 1.4;
end
function rgb = hex2rgb(h)
    rgb = zeros(numel(h),3);
    for k = 1:numel(h)
        s = h{k}(2:end); rgb(k,:) = [hex2dec(s(1:2)) hex2dec(s(3:4)) hex2dec(s(5:6))]/255;
    end
end
