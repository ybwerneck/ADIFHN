# start value for H
h1 = 117/360.0
# end value for H
h2 = 227/360.0
# creating the palette by specifying H,S,V
set palette model HSV functions (1-gray)*(h2-h1)+h1,1,0.68

