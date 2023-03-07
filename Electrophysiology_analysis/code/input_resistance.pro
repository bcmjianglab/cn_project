PRO Input_Resistance,data,Holding,First_Step,Step_Interval,Sweep_length,RMP=RMP,Rinput=Rinput

  Episode=N_elements(data(*,0))
  Sample=N_elements(data(0,*))+1.0000
  Sample_interval=Sweep_length/Sample

window,/free,title='Input Resistance and Resting Membrane Potential Calculation',xsize=1900,ysize=1000,xpos=0,ypos=0

n=Round((holding*1.0-First_Step)/Step_Interval) 
temp=Step_Interval*(n)+First_Step
n=n+1
if temp gt holding then begin
  n=n-1
endif





if n gt  5 then begin
  ;------------------------ n gt 5-----------------------------
  initial = median(reform(data(n-4,*)),1/ Sample_interval)


  plot,initial


  current = fltarr(5)
  voltage = fltarr(5)
  holdStart=599/Sample_interval
  holdEnd=698/Sample_interval

  for j = n-5 , n-1 do begin

    line_v = median(reform(data(j+1,*)),1/ Sample_interval)
    oplot, line_v
    current(j-n+5) = Step_Interval*(j)+First_step
    voltage(j-n+5) = mean(line_v(holdStart:holdEnd))

  endfor


  Plots, holdStart,[-180,100]
  Plots, holdEnd,[-180,100]
  wait,2
  plot, current,Voltage,linestyle=1,psym=4

  Fit=linfit(current,Voltage)
  Voltage1=Fit(1)*current+Fit(0)
  oPlot,current,Voltage1,color = 255

  Fit2=ladfit(current,Voltage)

  Voltage2=Fit2(1)*current+Fit2(0)
  oPlot,current,Voltage2,color='90EE90'x


  wait,0.2


endif else begin
  ;------------------------ n lt 5-----------------------------
  current = fltarr(n)
  voltage = fltarr(n)
  holdStart=599/Sample_interval
  holdEnd=698/Sample_interval
  line_v = median(reform(data(1,*)),1/ Sample_interval)
  plot, line_v
  for j = 0 , n-1 do begin
    line_v = median(reform(data(j+1,*)),1/ Sample_interval)
    oplot, line_v

    current(j) = Step_Interval*(j)+First_Step
    voltage(j) = mean(line_v(holdStart:holdEnd))

  endfor


  Plots, holdStart,[-180,100]

  Plots, holdEnd,[-180,100]
  wait,0.2

  plot, current,Voltage,linestyle=1,psym=4

  Fit=linfit(current,Voltage)
  Voltage1=Fit(1)*current+Fit(0)
  oPlot,current,Voltage1,color = 255

  Fit2=ladfit(current,Voltage)

  Voltage2=Fit2(1)*current+Fit2(0)
  oPlot,current,Voltage2,color='90EE90'x

  wait,0.2

endelse

print, current,voltage

  RMP=Fit2(0)
  Rinput=Fit2(1)*1000
  Print,'RMP(mV)',RMP
  Print,'Rm(MO) ',Rinput

cursor,xn,yn, /DOWN 
wdelete

Return

END