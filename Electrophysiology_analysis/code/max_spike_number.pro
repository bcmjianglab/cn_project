PRO Max_Spike_Number,data,holding,First_Step,Step_Interval,Sweep_length,Max_Spike_Num=Max_Spike_Num

Episode=N_elements(data(*,0))
Sample=N_elements(data(0,*))+1.0000
Sample_interval=Sweep_length/Sample

Step=fltarr(Episode-1)
Spike_Num=fltarr(Episode-1)

Step_Start=99/Sample_interval
Step_End=701/Sample_interval

For j= 1,Episode-1 do begin
line_a=reform(data(j,*));
line_a=line_a(Step_Start:Step_End)
smooth_step=1/Sample_interval
line_a=smooth(line_a,smooth_step)
line_da=shift(line_a,-1)-shift(line_a,1)
base=mean(line_da)+stddev(line_da)
line_da=smooth(line_da,1/Sample_interval)
line_da= line_da((1/Sample_interval):(600/Sample_interval))
Temppoint=10/sample_interval
Temp=WHERE(line_da eq max(line_da(0:temppoint)), /NULL)
Temp=Temp[0]


If temp eq 0 then begin
temp=1
endif

line_da(0:Temp-1)=0
If max(line_da) gt 0.5 then begin
  line_da[WHERE(line_da le base, /NULL)] = -100
  line_da[WHERE(line_da gt base, /NULL)] = 100
  line_da[WHERE(line_da eq -100, /NULL)] = 0 
  line_da[WHERE(line_da eq 100, /NULL)] = 1
endif else begin 
 line_da[WHERE(line_da ne 0, /NULL)] = 0
endelse

  For i=0, 599/Sample_interval-1 do begin
  If line_da(i) eq 0 and line_da(i+1) eq 1 then begin
  line_da(i)=1
  endif else begin
  line_da(i)=0
  Endelse
  Endfor
Step(j-1)=First_Step+Step_Interval*(j-1)-holding
if Step(j-1) le 0 then begin
  Spike_Num(j-1)=0
endif else begin
Spike_Num(j-1)=total(line_da)
Endelse
Endfor


Print,'max',max(Spike_Num)
Max_Spike_Num=max(Spike_Num)
return
END