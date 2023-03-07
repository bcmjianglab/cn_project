PRO Action_Potiential_Analysis,data,Holding,First_Step,Step_Interval,Sweep_length,$
    Rheo_Base=Rheo_Base,Spike_Delay=Spike_Delay,Threshold=Threshold,Overshoot=Overshoot,$
    APamp=APamp,AHP=AHP,Depotime=Depotime, Repotime= Repotime,half_width=half_width,RheoBase_Num=RheoBase_Num


Episode=N_elements(data(*,0))
Sample=N_elements(data(0,*))+1.0000 
Sample_interval=Sweep_length/Sample

window,/free,title='Action Potiential Analysis',xsize=1900,ysize=1000,xpos=0,ypos=0
win_plot=!d.window
!p.MULTI=[0,1,1,0,0];
Temp=(Holding-First_Step)/Step_Interval+1
plot,data(0,*),data(temp,*),yrange=[-180,100]
xyouts,100,900,'LEFT Click and RIGHT Click to select the RheoBase Sweep',charsize=2,/device

         
For i=temp, Episode-1 do begin

  plots,data(0,*),data(i,*)
  cursor,xn,yn, /DOWN
           if !mouse.button eq 4 then begin
           goto,flag1
           endif                    
Endfor
flag1:

xyouts,100,900,'LEFT Click and RIGHT Click to select the RheoBase Sweep',charsize=2,color=0,/device
xyouts,100,900,'LEFT Click and RIGHT Click to select the AP Sweep',charsize=2,/device
Rheo_Base=First_Step+Step_Interval*(i-1)-Holding
Print,'RheoBase(pA)',Rheo_Base
RheoBase_Num=i

for m=i, Episode-1 do begin
  plots,data(0,*),data(m,*)
  cursor,xn,yn, /DOWN
  if !mouse.button eq 4 then begin
    goto,flag3
  endif
endfor

flag3:
xyouts,100,900,'LEFT Click and RIGHT Click to select the AP Sweep',charsize=2,color=0,/device
win_plot=!d.window
!p.MULTI=[0,1,2,0,0]
line_v=data(m,*)
;line_v=data(i,*)
line_v=smooth(line_v,5);,edge=1)
line_dv=shift(line_v,0,-1)-shift(line_v,0,1)
line_dv=smooth(line_dv,5);,edge=1)
line_ddv=shift(line_dv,0,-1)-shift(line_dv,0,1)
line_ddv=smooth(line_ddv,5);,edge=1)
line_dddv=shift(line_ddv,0,-1)-shift(line_ddv,0,1)
plot,line_v
plot,line_dddv
xyouts,100,900,'Click to Zoom Action Potential Area',charsize=2,/device
cursor,x1,yn, /DOWN
cursor,x2,yn, /DOWN
  If x1 GT x2 Then Begin
  X2Temp=x2
  x2=x1
  x2=x2Temp               
  Endif

plot,line_v(x1:x2)
plot,line_dddv(x1:x2)
xyouts,100,900,'Click to Zoom Action Potential Area Again',charsize=2,Color=255,/device
cursor,x3,yn, /DOWN
cursor,x4,yn, /DOWN
  If x3 GT x4 Then Begin
  X4Temp=x4
  x4=x3
  x3=x4Temp               
  Endif

plot,line_v(x1+x3:x1+x4)
plot,line_dddv(x1+x3:x1+x4)
temp1=line_dddv(x1+x3:x1+x4)
xyouts,100,900,'Click to Zoom Action Potential Area Again',charsize=2,Color=0,/device
xyouts,100,900,'LEFT Click to NEXT',charsize=2,/device

point1=where(temp1 eq min(temp1))
point1=point1[0]
temp2=temp1(0:point1-1)
Point_Max=where(temp2 eq Max(temp2))
xr=Point_Max[0]
Plots,xr,[-180,100],color=255
cursor,xn,yn, /DOWN

xyouts,100,900,'LEFT Click to NEXT or RIGHT to select Threshold Point',charsize=2,/device
if !mouse.button eq 2 then begin
   xr=round(xn) 
endif else begin
   goto,flag2
endelse

flag2:


Threshold=line_v(xr+x1+x3)
Overshoot=max(line_v(x1:x2),index)
APamp=Overshoot-Threshold
   singleAP= line_v(x1:x2)
   overshoot_point=where(singleAP eq Max(singleAP))
   ahp_start= x1 + overshoot_point[0]
   AHP=min(line_v(ahp_start:x2),index)-Threshold

;AHP=min(line_v(x1:x2),index)-Threshold
Spike_Delay=[xr+x1+x3]*Sample_Interval-100
Print,'Spike_Delay(ms)',Spike_Delay
Print,'Threshold(mV)',Threshold
Print,'AP_Amplitude(mV)',APamp
Print,'AHP(mV)',AHP
half=APamp/2+Threshold
;print,half
;cursor,xn,yn, /DOWN
;win_plot=!d.window
!p.MULTI=[0,1,1,0,0]
!p.MULTI=[0,1,1,0,0]
plot,line_v(x1+x3:x1+x4)
TempAP=line_v(x1+x3:x1+x4)
Plots,xr,[-180,100],color='90EE90'x
Plots,[0,2000],[half,half],color='90EE90'x
Plots,[0,2000],[threshold,threshold],color='90EE90'x
  xyouts,100,900,'Check AP and Click to NEXT',charsize=2,/device
;  cursor,xn,yn, /DOWN
  AP_Start=xr
  Plots,AP_Start,[-180,100],color=255
  Overshoot_point=where(tempAP eq max(TempAP))
  Overshoot_point=Overshoot_point[0]
  Plots,Overshoot_point,[-180,100],color=255
  
  TempAP1=TempAP
  For i=0, round(x4-x3-2) do begin
  If TempAP1(i) le half and TempAP1(i+1) gt half then begin
  TempAP1(i)=1
  endif else begin
  TempAP1(i)=0
  Endelse
  Endfor
  AP_Half_Start=where(TempAP1 eq 1)
  AP_Half_Start=AP_Half_Start[0]
  x=[AP_Half_Start,AP_Half_Start+1]
  y=[tempAP(AP_Half_Start),tempAP(AP_Half_Start+1)]
  xx=(half-y[0])/(y[1]-y[0])
  AP_Half_Start=x[0]+xx
  
  Plots,AP_Half_Start,[-180,100],color=255
  
  TempAP2=TempAP
  For i=0, round(x4-x3-2) do begin
  If TempAP2(i) ge half and TempAP2(i+1) lt half then begin
  TempAP2(i)=1
  endif else begin
  TempAP2(i)=0
  Endelse
  Endfor
  AP_Half_End=where(TempAP2 eq 1)
  AP_Half_End=AP_Half_End[0]
  
  x=[ AP_Half_End, AP_Half_End+1]
  y=[tempAP(AP_Half_End),tempAP( AP_Half_End+1)]
  xx=(half-y[0])/(y[1]-y[0])
  AP_Half_End=x[0]+xx
  
  Plots,AP_Half_End,[-180,100],color=255

  TempAP3=TempAP
  For i=0, round(x4-x3-2) do begin
  If TempAP3(i) ge Threshold and TempAP3(i+1) lt Threshold then begin
  TempAP3(i)=1
  endif else begin
  TempAP3(i)=0
  Endelse
  Endfor
  AP_End=where(TempAP3 eq 1)
  AP_End=AP_End[0]
  x=[AP_End,AP_End+1]
  y=[tempAP(AP_End),tempAP(AP_End+1)]
  xx=(Threshold-y[0])/(y[1]-y[0])
  AP_End=x[0]+xx
  if AHP gt 0 then begin
   cursor,x_end,yn, /DOWN 
   AP_End=x_end
  endif
  
  
  Plots,AP_End,[-180,100],color=255
  
  Depotime=(Overshoot_point*1.0000-AP_Start*1.000)*Sample_interval
  Print,'Rise_phase Time(ms)',Depotime
  Repotime=(AP_End-Overshoot_point*1.0000)*Sample_interval
  Print,'Decay_phase Time(ms)',Repotime
  half_width=(AP_Half_End-AP_Half_Start)*Sample_interval
  Print,'Half Width(ms)',Half_width 
  ;xyouts,100,900,'Click to NEXT',charsize=2,/device 
cursor,xn,yn, /DOWN 
wdelete

Return

END