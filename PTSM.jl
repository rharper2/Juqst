using Match


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
    @match clifford begin
        1 => [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
        # case 2 % 1-1-S 
        2 =>  xtheta(pi/2 + deltaX)*ytheta(pi/2+deltaY);
        #case 3 % 1-1-S^2
        3 => ytheta(-pi/2 + deltaY)*xtheta(-pi/2+deltaX);
        #case 4 % x-1-1
        4 => xtheta(pi+deltaX);
        #case 5 % x-1-s
        5 => xtheta(-pi/2 + deltaX)*ytheta(-pi/2+deltaY);
        #case 6 % x-i-s^2
        6 => ytheta(-pi/2+deltaY)*xtheta(pi/2+deltaX);
        #case 7 % y-1-1
        7 => ytheta(pi+deltaY);
        #case 8 % y-1-s
         8 => xtheta(pi/2+deltaX)*ytheta(-pi/2+deltaY);
        # case 9 % y-1-s^2
        9 => ytheta(pi/2+deltaY)*xtheta(pi/2+deltaX);
        #case 10 % z-1-1
        10 => ytheta(pi+deltaY)*xtheta(pi+deltaX);
        #case 11 % z-i-s
        11 => xtheta(-pi/2+deltaX)*ytheta(pi/2+deltaY);
        #case 12 % z-i-s^2
        12 => ytheta(pi/2+deltaY)*xtheta(-pi/2+deltaX);
        #case 13 % i-h-i
        13 => xtheta(pi+deltaX)*ytheta(pi/2+deltaY);
        #case 14 % i-h-s
        14 => xtheta(-pi/2+deltaX);
        #case 15 % i-h-s^2
        15 => xtheta(-pi/2+deltaX2)*ytheta(-pi/2+deltaY)*xtheta(pi/2+deltaX);
        #case 16 % x-h-i
        16 => ytheta(-pi/2+deltaY);
        #case 17 % x-h-s
        17 => xtheta(pi/2+deltaX);
        #case 18 % x-h-s^2
        18 => xtheta(pi/2+deltaX2)*ytheta(pi/2+deltaY)*xtheta(pi/2+deltaX);
        #case 19 % y-h-1
        19 => xtheta(pi+deltaX)*ytheta(-pi/2+deltaY);
        #case 20 % y-h-s
        20 => ytheta(pi+deltaY)*xtheta(pi/2+deltaX);
        #case 21 % y-h-s^2
        21 => xtheta(pi/2+deltaX2)*ytheta(-pi/2+deltaY)*xtheta(pi/2+deltaX);
        #case 22 % z-h-i
        22 => ytheta(pi/2+deltaY);
        #case 23 % z-h-s
        23 => ytheta(pi+deltaY)*xtheta(-pi/2+deltaX);
        #case 24 % z-h-s^2
        24 => xtheta(pi/2+deltaX2)*ytheta(-pi/2+deltaY)*xtheta(-pi/2+deltaX);
    end


end




