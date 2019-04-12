function [A] = A_RectangleWin_CT(R, filter)

% R must be normalised between [-1, 1]
switch filter
    case 'ram-lak'
        % Do nothing
        A = 1;
    case 'shepp-logan'
        % be careful not to divide by 0:
        A = sin(pi.*R./2)./(pi.*R./2);
    case 'cosine'
        A = cos(pi.*R./2);
    case 'hamming'
        alpha = 0.54;
        A = alpha + (1-alpha) .* cos(pi.*R);
    case 'hann'
        A = (1+cos(pi.*R))./ 2;
    case 'rect'
        A = double(R<=1);
end
%A = A.*(abs(R)<=1);

return
