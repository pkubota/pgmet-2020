'reinit'
'set display color white'
'c'
'open saida_GridB.ctl'
'set gxout shaded'
it=1
while(it<=3000   )
'rgbset_local.gs'
'set clevs    0.996   0.9965   0.997    0.9975   0.998     0.999  0.9995     1.0000  1.0005 1.001 '
'set ccols 21      22        24      26       28        29      49        0         0     44      42  41    '
'd h(t='it')'
'!sleep 1'
say it
it=it+20
endwhile

