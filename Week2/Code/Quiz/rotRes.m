
function rotRes = rotRes(A,axis,radRot)

%radRot = degree * pi / 180;
if axis == 'x'
    rotRes = A * [1 0 0; 0 cos(radRot) -sin(radRot); 0 sin(radRot) cos(radRot)];
    display('RotX')
elseif axis == 'y'
    rotRes = A * [cos(radRot) 0 sin(radRot); 0 1 0; -sin(radRot) 0 cos(radRot)]; 
    display('RotY')
elseif axis == 'z'
    rotRes = A * [cos(radRot) -sin(radRot) 0 ; sin(radRot) cos(radRot) 0; 0 0 1];    
    display('RotZ')
end

end

%%