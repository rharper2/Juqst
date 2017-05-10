
function xtheta( angle )
#returns the PTSM representation for X rotation of angle
#   Returns 
#   (1     0       0              0  
#    0     1                      0
#    0     0   cos(theta)   -sin(theta)
#    0     0   sin(theta)    cos(theta)    )

 [1 0 0 0;0 1 0 0;0 0 cos(angle) -sin(angle);0 0 sin(angle) cos(angle)];

end

function ytheta( angle )
#returns the PTSM representation for Y rotation of angle
#   Returns 
#   (1     0             0              0  
#    0     cos(theta)    0            sin(theta)
#    0        0          1              0
#    0    -sin(theta)    0            cos(theta)    )

 [1 0 0 0;0 cos(angle) 0 sin(angle); 0 0 1 0; 0 -sin(angle) 0 cos(angle)];

end


function generateClifford( clifford,deltaX,deltaY,deltaX2)
# generateClifford(Clifford, deltaX,deltaY,deltaX2), 
#   generates one qubit Cliffords, from 1..24 using
#   elementary rotation gates, allowing an over-rotation to be supplied
#   deltaX2 happens where there is a second xgate needed (e..g in the
#   I-H-S^2 gate (gate 15).
#   
#   Reorderd to match the symplectic generated cliffords 
    if clifford == 24
        return  [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
    elseif clifford ==14
        return  xtheta(pi/2 + deltaX)*ytheta(pi/2+deltaY);
        #case 3 % 1-1-S^2
    elseif clifford == 5
        return ytheta(-pi/2 + deltaY)*xtheta(-pi/2+deltaX);
        #case 4 % x-1-1
    elseif clifford == 22 
        return xtheta(pi+deltaX);
        #case 5 % x-1-s
    elseif clifford == 13 
        return xtheta(-pi/2 + deltaX)*ytheta(-pi/2+deltaY);
        #case 6 % x-i-s^2
    elseif clifford == 8 
        return ytheta(-pi/2+deltaY)*xtheta(pi/2+deltaX);
        #case 7 % y-1-1
    elseif clifford == 23 
        return ytheta(pi+deltaY);
        #case 8 % y-1-s
    elseif clifford ==  15 
        return xtheta(pi/2+deltaX)*ytheta(-pi/2+deltaY);
        # case 9 % y-1-s^2
    elseif clifford == 7
        return ytheta(pi/2+deltaY)*xtheta(pi/2+deltaX);
        #case 10 % z-1-1
    elseif clifford == 21 
        return ytheta(pi+deltaY)*xtheta(pi+deltaX);
        #case 11 % z-i-s
    elseif clifford == 16 
        return xtheta(-pi/2+deltaX)*ytheta(pi/2+deltaY);
        #case 12 % z-i-s^2
    elseif clifford == 6 
        return ytheta(pi/2+deltaY)*xtheta(-pi/2+deltaX);
        #case 13 % i-h-i
    elseif clifford == 4 
        return xtheta(pi+deltaX)*ytheta(pi/2+deltaY);
        #case 14 % i-h-s
    elseif clifford == 12 
        return xtheta(-pi/2+deltaX);
        #case 15 % i-h-s^2
    elseif clifford == 20 
        return xtheta(-pi/2+deltaX2)*ytheta(-pi/2+deltaY)*xtheta(pi/2+deltaX);
        #case 16 % x-h-i
    elseif clifford == 1 
        return ytheta(-pi/2+deltaY);
        #case 17 % x-h-s
    elseif clifford == 10 
        return xtheta(pi/2+deltaX);
        #case 18 % x-h-s^2
    elseif clifford == 19 
        return xtheta(pi/2+deltaX2)*ytheta(pi/2+deltaY)*xtheta(pi/2+deltaX);
        #case 19 % y-h-1
    elseif clifford == 3 
        return xtheta(pi+deltaX)*ytheta(-pi/2+deltaY);
        #case 20 % y-h-s
    elseif clifford == 9 
        return ytheta(pi+deltaY)*xtheta(pi/2+deltaX);
        #case 21 % y-h-s^2
    elseif clifford == 18 
        return xtheta(pi/2+deltaX2)*ytheta(-pi/2+deltaY)*xtheta(pi/2+deltaX);
        #case 22 % z-h-i
    elseif clifford == 2 
        return ytheta(pi/2+deltaY);
        #case 23 % z-h-s
    elseif clifford == 11 
        return ytheta(pi+deltaY)*xtheta(-pi/2+deltaX);
        #case 24 % z-h-s^2
    elseif clifford == 17
        return xtheta(pi/2+deltaX2)*ytheta(-pi/2+deltaY)*xtheta(-pi/2+deltaX);
    end

end




