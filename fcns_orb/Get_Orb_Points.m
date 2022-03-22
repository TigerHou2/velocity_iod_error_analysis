function pos = Get_Orb_Points(r,v,mu,res,dur,start_day)
%GET_ORB_POINTS Summary of this function goes here
%   Detailed explanation goes here

[a,e,i,omg,w] = Get_Orb_Params(r,v,mu);
e = norm(e);

pos = [];

if a > 0
    for f = 0:pi/4:2*pi
        r1 = a*(1-e^2) / (1+e*cos(f));
        x = r1 * (cos(omg)*cos(w+f)-sin(omg)*sin(w+f)*cos(i));
        y = r1 * (sin(omg)*cos(w+f)+cos(omg)*sin(w+f)*cos(i));
        z = r1 * sin(w+f)*sin(i);
        pos = [pos; [x, y, z, f, 1]];
    end
    % adaptive sampling
    while true
        refine = [];
        for m = 1:size(pos,1)
            if m == size(pos,1)
                n = 1;
                f = (pos(m,4) + pos(n,4) + 2*pi) / 2;
            else
                n = m + 1;
                f = (pos(m,4) + pos(n,4)) / 2;
            end
            if pos(m,5) == 0 && pos(n,5) == 0
                continue
            end
            rc = a*(1-e^2) / (1+e*cos(f));
            x = rc * (cos(omg)*cos(w+f)-sin(omg)*sin(w+f)*cos(i));
            y = rc * (sin(omg)*cos(w+f)+cos(omg)*sin(w+f)*cos(i));
            z = rc * sin(w+f)*sin(i);
            rc = [x y z];
            r1 = pos(m,1:3);
            r2 = pos(n,1:3);
            v1 = rc - r1;
            v2 = r2 - rc;
            turn_angle = atan2d(norm(cross(v1,v2)), dot(v1,v2));
            if abs(turn_angle) > 3
                refine = [[rc, f, 1, m]; refine];
                pos(m,5) = 1;
                pos(n,5) = 1;
            else
                pos(m,5) = 0;
                pos(n,5) = 0;
            end
        end
        if isempty(refine)
            break
        end
        for q = 1:size(refine,1)
            p = refine(q,6);
            pos = [pos(1:p,:); refine(q,1:5); pos(p+1:end,:)];
        end
    end
else
    for i = start_day-dur:dur/res:start_day+dur
        r1 = TimeProp_V4(r, v, mu, i);
        pos = [pos; r1];
    end
end

end
