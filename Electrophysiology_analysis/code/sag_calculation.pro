PRO Sag_calculation,data,sweep_length,Sag_Ratio=Sag_Ratio
Episode=N_elements(data(*,0))
Sample=N_elements(data(0,*))+1.0000
Sample_interval=Sweep_length/Sample

window,/free,title='Sag Ratio Calculation',xsize=1900,ysize=1000,xpos=0,ypos=0


line_a=data(1,*)
line_a=smooth(line_a,1/Sample_interval)
plot,line_a
Base_Start=100
Plots, Base_Start,[-180,100]
Base_end=100/Sample_interval-100
Plots, Base_End,[-180,100]
sag_Start=100/Sample_interval+100
Plots, sag_Start,[-180,100]
sag_end=300/Sample_interval-100
Plots, sag_End,[-180,100]

flat_end=700/Sample_interval-100
Plots, flat_end,[-180,100]
flat_start=600/Sample_interval+100
Plots, flat_start,[-180,100]

base=mean(line_a(Base_Start:Base_end))
sagtemp=min(line_a(sag_Start:sag_end))
flat=mean(line_a(flat_Start:flat_end))

Sag_Ratio=(sagtemp-base)/(flat-base)
Print,'Sag Ratio ',Sag_Ratio

cursor,xn,yn, /DOWN 
wdelete

Return
END