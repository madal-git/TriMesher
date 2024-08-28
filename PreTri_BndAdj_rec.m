function [NH, NL]=PreTri_BndAdj_rec(H,L,NH,NL,s)
% function [NH, NL]=Rec_nd_density[H,L,NH,NL] tunes the node density on a
% rectangular boundaries in order to make sure the node spacing in the
% vertical and horizontal directions are almost the same.
% NH is the total number of nodes in the vertical direction.
% NL is the total number of nodes in the horizontal direction.
% H is the height of the rectangle.
% L is the length of the rectangle.
% s is the control parameter. s=0----Condense the loose spacing;
%  s=1----Loose the condense spacing;
%  s=2----Average spacing;
C=0.2; % Criteria for determing whether tuning spacing is needed.

dx=L/(NL-1);
dy=H/(NH-1);

if dx>dy
    dl=dx;
    ds=dy;
else
    dl=dy;
    ds=dx;
end

if (dl-ds)/ds<=C % The spacing in x  and y direction is OK, no need for tuning
    ;
else
    if dx<dy % The spacing in y direction is much larger than x direction
        if s==0
            while (dl-ds)/ds>C
                NH=NH+1;
                dy=H/(NH-1);
                if dx>dy
                    dl=dx;
                    ds=dy;
                else
                    dl=dy;
                    ds=dx;
                end
            end
        elseif s==1
            while (dl-ds)/ds>C
                NL=NL-1;
                dx=L/(NL-1);
                if dx>dy
                    dl=dx;
                    ds=dy;
                else
                    dl=dy;
                    ds=dx;
                end
            end
        elseif s==2
            while (dl-ds)/ds>C
                NL=NL-1;
                dx=L/(NL-1);
                NH=NH+1;
                dy=H/(NH-1);
                if dx>dy
                    dl=dx;
                    ds=dy;
                else
                    dl=dy;
                    ds=dx;
                end
            end
        else
            error('The tag for tuning spacing is incorrect');
        end
    else % The spacing in x direction is much larger than y direction
        if s==0
            while (dl-ds)/ds>C
                NL=NL+1;
                dx=L/(NL-1);
                if dx>dy
                    dl=dx;
                    ds=dy;
                else
                    dl=dy;
                    ds=dx;
                end
            end
        elseif s==1
            while (dl-ds)/ds>C
                NH=NH-1;
                dy=L/(NH-1);
                if dx>dy
                    dl=dx;
                    ds=dy;
                else
                    dl=dy;
                    ds=dx;
                end
            end
        elseif s==2
            while (dl-ds)/ds>C
                NL=NL+1;
                dx=L/(NL-1);
                NH=NH-1;
                dy=H/(NH-1);
                if dx>dy
                    dl=dx;
                    ds=dy;
                else
                    dl=dy;
                    ds=dx;
                end
            end
        else
            error('The tag for tuning spacing is incorrect');
        end
    end
end