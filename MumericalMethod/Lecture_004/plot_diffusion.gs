'reinit'
'set display color white'
'c'
'open ImplicitLinearDiffusion1D.ctl'

it=0
while(it<=1000)
'c'
'd a(t=1)'
'd a(t='it')'
if(it=999)
  'q pos'  
  it=0
endif
it=it+1
endwhile
