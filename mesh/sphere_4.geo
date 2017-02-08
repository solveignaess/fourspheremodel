Mesh.Algorithm = 5;

Function CheeseHole 
  // taken from t5.geo

  p1 = newp; Point(p1) = {x,  y,  z,  lcar3} ;
  p2 = newp; Point(p2) = {x+r,y,  z,  lcar3} ;
  p3 = newp; Point(p3) = {x,  y+r,z,  lcar3} ;
  p4 = newp; Point(p4) = {x,  y,  z+r,lcar3} ;
  p5 = newp; Point(p5) = {x-r,y,  z,  lcar3} ;
  p6 = newp; Point(p6) = {x,  y-r,z,  lcar3} ;
  p7 = newp; Point(p7) = {x,  y,  z-r,lcar3} ;

  c1 = newreg; Circle(c1) = {p2,p1,p7};
  c2 = newreg; Circle(c2) = {p7,p1,p5};
  c3 = newreg; Circle(c3) = {p5,p1,p4};
  c4 = newreg; Circle(c4) = {p4,p1,p2};
  c5 = newreg; Circle(c5) = {p2,p1,p3};
  c6 = newreg; Circle(c6) = {p3,p1,p5};
  c7 = newreg; Circle(c7) = {p5,p1,p6};
  c8 = newreg; Circle(c8) = {p6,p1,p2};
  c9 = newreg; Circle(c9) = {p7,p1,p3};
  c10 = newreg; Circle(c10) = {p3,p1,p4};
  c11 = newreg; Circle(c11) = {p4,p1,p6};
  c12 = newreg; Circle(c12) = {p6,p1,p7};

  l1 = newreg; Line Loop(l1) = {c5,c10,c4};   Ruled Surface(newreg) = {l1};
  l2 = newreg; Line Loop(l2) = {c9,-c5,c1};   Ruled Surface(newreg) = {l2};
  l3 = newreg; Line Loop(l3) = {c12,-c8,-c1}; Ruled Surface(newreg) = {l3};
  l4 = newreg; Line Loop(l4) = {c8,-c4,c11};  Ruled Surface(newreg) = {l4};
  l5 = newreg; Line Loop(l5) = {-c10,c6,c3};  Ruled Surface(newreg) = {l5};
  l6 = newreg; Line Loop(l6) = {-c11,-c3,c7}; Ruled Surface(newreg) = {l6};
  l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12};Ruled Surface(newreg) = {l7};
  l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};  Ruled Surface(newreg) = {l8};

  theloops[t] = newreg ; 
  Surface Loop(theloops[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};

  // surf_sphere = newreg;
  // Surface Loop(surf_sphere) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};

  //thehole = newreg ; 
  //Volume(thehole) = theloops[t] ;
  // Volume(thehole) = surf_sphere ;

Return

lcar3 = 1.5;
x = 0.; y= 0.; z=0.; r=6.0; t=1;
Call CheeseHole ;
whitemattersurf = news;
Physical Surface (whitemattersurf) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1}; 
whitemattervolume = newreg;
Volume(whitemattervolume) = {theloops[1]};
whitemattervol = newreg;
Physical Volume (whitemattervol) = whitemattervolume;

lcar3 = .08;
x = 0.; y= 0.; z=0.; r=7.9; t=2;
Call CheeseHole ;
graymattersurf = news;
Physical Surface (graymattersurf) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1}; 
graymattervolume = newreg;
Volume(graymattervolume) = {theloops[2], theloops[1]};
graymattervol = newreg;
Physical Volume (graymattervol) = graymattervolume;

lcar3 = .03;
x = 0.; y= 0.; z=0.; r=8.; t=3;
Call CheeseHole ;
csfsurf = news;
Physical Surface (csfsurf) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1}; 
csfvolume = newreg;
Volume(csfvolume) = {theloops[3], theloops[2]};
csfvol = newreg;
Physical Volume (csfvol) = csfvolume;

lcar3 = .25;
x = 0.; y= 0.; z=0.; r=8.5; t=4;
Call CheeseHole ;
skullsurf = news;
Physical Surface (skullsurf) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1}; 
skullvolume = newreg;
Volume(skullvolume) = {theloops[4], theloops[3]};
skullvol = newreg;
Physical Volume (skullvol) = skullvolume;

lcar3 = .25;
x = 0.; y= 0.; z=0.; r=9.; t=5;
Call CheeseHole ;
scalpsurf = news;
Physical Surface (scalpsurf) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1}; 
scalpvolume = newreg;
Volume(scalpvolume) = {theloops[5], theloops[4]};
scalpvol = newreg;
Physical Volume (scalpvol) = scalpvolume;

// lcar3 = 10.;
// x = 0.; y= 0.; z=0.; r=50.; t=5;
// Call CheeseHole ;
// airsurf = news;
// Physical Surface (airsurf) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1}; 
// airvolume = newreg;
// Volume(airvolume) = {theloops[5], theloops[4]};
// airvol = newreg;
// Physical Volume (airvol) = airvolume;

