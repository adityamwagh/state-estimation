f = [
 
                                                     uPrev(7);
                                                     uPrev(8);
                                                     uPrev(9);
 -cos(uPrev(4))*cos(uPrev(6))*sin(uPrev(5))*(uPrev(10) + n(1) - angVel(1));
  cos(uPrev(5))*sin(uPrev(4))*sin(uPrev(6))*(uPrev(11) + n(2) - angVel(2));
(cos(uPrev(6))*sin(uPrev(5))*(uPrev(12) + n(3) - angVel(3)))/cos(uPrev(5));
             -cos(uPrev(4))*cos(uPrev(5))*(uPrev(13) - acc(1) + n(4));
             -cos(uPrev(5))*sin(uPrev(4))*(uPrev(14) - acc(2) + n(5));
                sin(uPrev(5))*(uPrev(15) - acc(3) + n(6)) - 981/100;
                                                   n(7);
                                                   n(8);
                                                   n(9);
                                                   n(10);
                                                   n(11);
                                                   n(12);
];

U_t = [
                                   0,                                   0,                                     0,                        0,                        0,           0, 0, 0, 0, 0, 0, 0;
                                   0,                                   0,                                     0,                        0,                        0,           0, 0, 0, 0, 0, 0, 0;
                                   0,                                   0,                                     0,                        0,                        0,           0, 0, 0, 0, 0, 0, 0;
-cos(uPrev(4))*cos(uPrev(6))*sin(uPrev(5)),                                   0,                                     0,                        0,                        0,           0, 0, 0, 0, 0, 0, 0;
                                   0, cos(uPrev(5))*sin(uPrev(4))*sin(uPrev(6)),                                     0,                        0,                        0,           0, 0, 0, 0, 0, 0, 0;
                                   0,                                   0, (cos(uPrev(6))*sin(uPrev(5)))/cos(uPrev(5)),                        0,                        0,           0, 0, 0, 0, 0, 0, 0;
                                   0,                                   0,                                     0, -cos(uPrev(4))*cos(uPrev(5)),                        0,           0, 0, 0, 0, 0, 0, 0;
                                   0,                                   0,                                     0,                        0, -cos(uPrev(5))*sin(uPrev(4)),           0, 0, 0, 0, 0, 0, 0;
                                   0,                                   0,                                     0,                        0,                        0, sin(uPrev(5)), 0, 0, 0, 0, 0, 0;
                                   0,                                   0,                                     0,                        0,                        0,           0, 1, 0, 0, 0, 0, 0;
                                   0,                                   0,                                     0,                        0,                        0,           0, 0, 1, 0, 0, 0, 0;
                                   0,                                   0,                                     0,                        0,                        0,           0, 0, 0, 1, 0, 0, 1;
                                   0,                                   0,                                     0,                        0,                        0,           0, 0, 0, 0, 1, 0, 0;
                                   0,                                   0,                                     0,                        0,                        0,           0, 0, 0, 0, 0, 1, 0;
                                   0,                                   0,                                     0,                        0,                        0,           0, 0, 0, 0, 0, 0, 0;
];
