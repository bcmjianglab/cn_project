PRO Spiking_frequency_adaption, data,Holding,First_Step,Step_Interval,Sweep_length,Rheo_Base,$ 
AP_Num = AP_Num,Adapt_Index_Initial = Adapt_Index_Initial,Adapt_Index_End = Adapt_Index_End,AP2_Diff = AP2_Diff,$
AP3_Diff = AP3_Diff, AP_End_Diff = AP_End_Diff


Episode=N_elements(data(*,0))
Sample=N_elements(data(0,*))+1.0000
Sweep_length=1000
Sample_interval=Sweep_length/Sample

  IIReobase=2*Rheo_base
 ;print, IIReobase
 n=Round((IIReobase*1.0+holding-First_Step)/Step_Interval)+1
 ;print,n
 n=n+1
 
 if n gt Episode-1 then begin
  
  n=Episode-1
 
 endif
window,/free,title='adaption Calculation',xsize=1900,ysize=1000,xpos=0,ypos=0 
 plot,data(n,*)
 wait,0.2
 
 
 Step_Start=100/Sample_interval
 Step_End=700/Sample_interval
 
 line_v=reform(data(n,*))
 Base=line_v(Step_Start:Step_End)
 Base=mean(base)+2.53*stddev(base)
 Origin=data(n,*)
 line_v=data(n,*)
 smooth_step=1/Sample_interval
 smooth_step=5
 line_v=smooth(line_v,smooth_step)
 line_v[WHERE(line_v lt base, /NULL)] = base
 line_v=line_v-base
 line_dv=shift(line_v,0,-1)-shift(line_v,0,1)
 line_dv=smooth(line_dv,smooth_step);,edge=1)


 For m=0, Sample-3 do begin
   If line_dv(m) ge 0 and line_dV(m+1) lt 0 then begin
     line_dv(m)=1
   endif else begin
     line_dv(m)=0
   Endelse
 Endfor

 ;!p.MULTI=[0,1,1,0,0]
 plot,Origin
 orign=reform(origin)
 dv=shift(origin,0,-1)-shift(origin,0,1)
 ;plot,dv
 marker=where(line_dv eq 1)
 marker=marker(where( marker gt 101/Sample_interval))
 marker=marker(where( marker lt 701/Sample_interval))
 ;print,N_elements(marker)
 num1=N_elements(marker)
 
;print,marker
 For k=0,num1-1 do begin
   plots,marker(k),[-180,100],thick=2,color='90EE90'x
   ;print,max(orign[(marker(k)-5):(marker(k)+5)])
 Endfor

 if num1 eq 1 then begin
   Adapt1='NaN'  
   Adapt2='NaN'
   print,Adapt1
   print,Adapt2
   AP1=max(orign[marker(0)-5:marker(0)+5])
   AP2='NaN'
   AP3='NaN'
   AP_END= 'NaN'  
   AP2_Diff='NaN'
   AP3_Diff='NaN'
   AP_End_Diff='NaN'   
 endif
 
  if num1 eq 2 then begin
   Adapt1='NaN'
   Adapt2='NaN'
   ;print,Adapt1
   ;print,Adapt2
 AP1=max(orign[marker(0)-5:marker(0)+5])
 AP2=max(orign[marker(1)-5:marker(1)+5])
 AP3='NaN'
 AP_END= 'NaN'
 AP2_Diff= AP1-AP2
 AP3_Diff= 'NaN'
 AP_End_Diff= 'NaN'
  endif
  
  if num1 ge 3 then begin
    Adapt1=[marker(2)-marker(1)]*1.0000/[marker(1)-marker(0)]
    Adapt2=[marker(num1-1)-marker(num1-2)]*1.0000/[marker(1)-marker(0)]
    print,Adapt1
    print,Adapt2
    AP1=max(orign[marker(0)-5:marker(0)+5])
    AP2=max(orign[marker(1)-5:marker(1)+5])
    AP3=max(orign[marker(2)-5:marker(2)+5])
    AP_END=max(orign[marker(num1-1)-5:marker(num1-1)+5])
    AP2_Diff= AP1-AP2
    AP3_Diff= AP1-AP3
    AP_End_Diff= AP1-AP_END
  endif 
  
print,'AP_Num@2xRheobase',';','Adapt_Index_Initial',';','Adapt_Index_End',';','AP2_Diff',';','AP3_Diff',';','AP_End_Diff'
print,num1,Adapt1,Adapt2,AP2_Diff,AP3_Diff,AP_End_Diff, FORMAT= '(10(A,";"))'
cursor,xn,yn, /DOWN 
wdelete

AP_Num = num1
Adapt_Index_Initial = Adapt1
Adapt_Index_End = Adapt2
AP2_Diff = AP2_Diff
AP3_Diff = AP3_Diff
AP_End_Diff = AP_End_Diff
Return

END