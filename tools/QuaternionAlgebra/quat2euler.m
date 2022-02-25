function euler=quat2euler(quat)
%Copyright (C) 2022 by Frida Viset

q0=quat(1,:); q1=quat(2,:); q2=quat(3,:); q3=quat(4,:);
roll=atan2(2.*(q0.*q1+q2.*q3),1-2.*(q1.^2+q2.^2));
pitch=asin(2.*(q0.*q2-q3.*q1));
yaw=atan2(2.*(q0.*q3+q1.*q2),1-2.*(q2.^2+q3.^2));

for t=2:length(roll)
    

    
    if roll(t)-roll(t-1)>5
        roll(t:end)=roll(t:end)-2*pi;
    elseif roll(t)-roll(t-1)<-5
        roll(t:end)=roll(t:end)+2*pi;
    end
    
    if yaw(t)-yaw(t-1)>5
        yaw(t:end)=yaw(t:end)-2*pi;
    elseif yaw(t)-yaw(t-1)<-5
        yaw(t:end)=yaw(t:end)+2*pi;
    end
    
    if pitch(t)-pitch(t-1)>5
        pitch(t:end)=pitch(t:end)-2*pi;
    elseif pitch(t)-pitch(t-1)<-5
        pitch(t:end)=pitch(t:end)+2*pi;
    end
    
    
end

euler=[roll; pitch; yaw];
end