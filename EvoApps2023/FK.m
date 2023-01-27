function [A] = FK(angle)   
% forward kinematics
    th1 = angle(1);
    th2 = angle(2);
    th3 = angle(3);
    th4 = angle(4);
    th5 = angle(5);
    th6 = angle(6);

    o1 = pi/2;
    d1 = 0.1519;
    a1 = 0;
    
    o2 = 0;
    d2 = 0;
    a2 = -0.24365;

    o3 = 0;
    d3 = 0;
    a3 = -0.21325;
    
    o4 = pi/2;
    d4 = 0.11235;
    a4 = 0;

    o5 = -pi/2;
    d5 = 0.08535;
    a5 = 0;
    
    o6 = 0;
    d6 = 0.0819;
    a6 = 0;

    A1 = [ cos(th1)  -sin(th1) * cos(o1)     sin(th1) * sin(o1)     a1 * cos(th1) ;
           sin(th1)   cos(th1) * cos(o1)    -cos(th1) * sin(o1)     a1 * sin(th1) ;
           0         sin(o1)               cos(o1)               d1           ;
           0         0                     0                     1          ] ;

    A2 = [ cos(th2)  -sin(th2) * cos(o2)     sin(th2) * sin(o2)     a2 * cos(th2) ;
           sin(th2)   cos(th2) * cos(o2)    -cos(th2) * sin(o2)     a2 * sin(th2) ;
           0         sin(o2)               cos(o2)               d2           ;
           0         0                     0                     1          ] ;

    A3 = [ cos(th3)  -sin(th3) * cos(o3)     sin(th3) * sin(o3)     a3 * cos(th3) ;
           sin(th3)   cos(th3) * cos(o3)    -cos(th3) * sin(o3)     a3 * sin(th3) ;
           0         sin(o3)               cos(o3)               d3           ;
           0         0                     0                     1          ] ;

    A4 = [ cos(th4)  -sin(th4) * cos(o4)     sin(th4) * sin(o4)     a4 * cos(th4) ;
           sin(th4)   cos(th4) * cos(o4)    -cos(th4) * sin(o4)     a4 * sin(th4) ;
           0         sin(o4)               cos(o4)               d4           ;
           0         0                     0                     1          ] ;

    A5 = [ cos(th5)  -sin(th5) * cos(o5)     sin(th5) * sin(o5)     a5 * cos(th5) ;
           sin(th5)   cos(th5) * cos(o5)    -cos(th5) * sin(o5)     a5 * sin(th5) ;
           0         sin(o5)               cos(o5)               d5           ;
           0         0                     0                     1          ] ;

    A6 = [ cos(th6)  -sin(th6) * cos(o6)     sin(th6) * sin(o6)     a6 * cos(th6) ;
           sin(th6)   cos(th6) * cos(o6)    -cos(th6) * sin(o6)     a6 * sin(th6) ;
           0         sin(o6)               cos(o6)               d6           ;
           0         0                     0                     1          ] ;
    
    A = A1 * A2 * A3 * A4 * A5 * A6;
end

