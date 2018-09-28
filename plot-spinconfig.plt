
set size ratio 1
unset key

set palette model RGB defined ( -1 'grey', 1 'black' )
unset colorbox

plot 'spinconf_n8000_T20.000000' u 1:2:($3<0.0001?-1:1) w points palette pt 7 ps 1
