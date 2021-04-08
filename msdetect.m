function [sac,v] = msdetect(data,VFAC,MINDUR,srate)

% vars
slength = 2;

data(1,:) = tn_rmvlinenoise(data(1,:),50,2,srate);
data(2,:) = tn_rmvlinenoise(data(2,:),50,2,srate);


data(1,:) = smoothdata(data(1,:),2,'movmean',5);
data(2,:) = smoothdata(data(2,:),2,'movmean',5);

% create temp data
temp_data(1,:) = data(1,:);
temp_data(2,:) = data(2,:);

% transform into velocity space/differentiation
v = [];
for idx = (slength+1):(size(temp_data,2)-slength)
    v(:,idx) = sum(-temp_data(:,idx-[1:slength]) + temp_data(:,idx+[1:slength]),2) ./ (2*sum(1:slength));   
end

% compute SD based on median estimator
medx = median(v(1,:));
msdx = sqrt(median((v(1,:)-medx).^2) );
medy = median(v(2,:));
msdy = sqrt(median((v(2,:)-medy).^2) );

% elliptical threshold
radiusx = VFAC * msdx;
radiusy = VFAC * msdy;

% find possible saccadic events
test = (v(1,:)/radiusx).^2 + (v(2,:)/radiusy).^2;
indx = find(test > 1);



% detect MSs here
N    = length(indx); 
sac  = [];
nsac = 0;
dur  = 1;
a    = 1;
k    = 1;

while k<N
    if indx(k+1)-indx(k)==1
        dur = dur + 1;
    else
        % duration > MINDUR?
        if (dur >= MINDUR) 
            nsac = nsac + 1;
            b = k;
            sac = [sac;[indx(a),indx(b),repmat(0,1,5)]];
        end
        a = k+1;
        dur = 1;
    end
    k = k + 1;
end
% last saccade: duration > MINDUR?
if (dur >= MINDUR) 
    nsac = nsac + 1;
    b = k;
    sac = [sac;[indx(a),indx(b),repmat(0,1,5)]];
end

% create sac vector with saccade information
v = v';
data = data';

if ( nsac>0 )
    % Compute peak velocity, horiztonal and vertical components
    for s = 1:nsac
        % Onset and offset for saccades
        a = sac(s,1);
        b = sac(s,2);
        idx = a:b;
        % Saccade peak velocity (vpeak)
        vpeak = max( sqrt( v(idx,1).^2 + v(idx,2).^2 ) );
        sac(s,3) = vpeak;
        % Saccade vector (dx,dy)
        dx = data(b,1)-data(a,1);
        dy = data(b,2)-data(a,2);
        sac(s,4:5) = [dx,dy];
        % Saccade amplitude (dX,dY)
        minx = min(data(idx,1));
        maxx = max(data(idx,1));
        miny = min(data(idx,2));
        maxy = max(data(idx,2));
        [~,ix1] = min(data(idx,1));
        [~,ix2] = max(data(idx,1));
        [~,iy1] = min(data(idx,2));
        [~,iy2] = max(data(idx,2));
        dX = sign(ix2-ix1)*(maxx-minx);
        dY = sign(iy2-iy1)*(maxy-miny);
        sac(s,6:7) = [dX,dY];
    end
end

end



