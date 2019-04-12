% LineMask.m
%
% Returns the indicator of the domain in 2D fourier space for the 
% specified line geometry.
% Usage :  [M,Mh,mi,mhi] = LineMask(L,N)
%
% Written by : Justin Romberg
% Created : 1/26/2004
% Revised : 12/2/2004

function [M,Mh,mi,mhi, lineList] = myLineMask(L,N)


thc = linspace(0, pi-pi/L, L);
%thc = linspace(pi/(2*L), pi-pi/(2*L), L);

M = sparse(N);
lineList = cell(L, 1);

% full mask
for ll = 1:L
    MM = sparse(N);
	if ((thc(ll) <= pi/4) | (thc(ll) > 3*pi/4))
		yr = round(tan(thc(ll))*(-N/2:N/2-1))+N/2+1;
        if max(yr)>N
            [v, p] = max(yr);
            yr(p) = 1;
        end
        
    	for nn = 1:N
            M(yr(nn),nn) = ll*10;
            MM(yr(nn),nn) = 1;
        end
  else 
		xc = round(cot(thc(ll))*(-N/2:N/2-1))+N/2+1;
        if max(xc)>N
            [v, p] = max(xc);
            xc(p) = 1;
        end
        
		for nn = 1:N
			M(nn,xc(nn)) = ll*10;
            MM(nn,xc(nn)) = 1;
		end
    end
    
    [xx, yy] = find(MM>0); 
    lineList{ll} = [xx, yy];
    index_line = lineList{ll};
    if (length(index_line)> N)
        lineList{ll} = index_line(2:end, :);
	end
end


% upper half plane mask (not including origin)
Mh = sparse(N);
Mh = M;
Mh(N/2+2:N,:) = 0;
Mh(N/2+1,N/2+1:N) = 0;


M = ifftshift(M);
mi = find(M);
Mh = ifftshift(Mh);
mhi = find(Mh);
