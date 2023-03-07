PRO Ephysi_analysis


;  holding = -100.0 ; define the holding current (pA)
;  First_Step = -200.0 ; define relative difference between injected currents of the first sweep and holding current (pA)
;  Step_Interval= 20.0 ; the step differences between steps(pA)
;  Sweep_length=1000 ; length of a sweep (ms)
;  First_Step =  holding + First_Step ;
  
  file=Dialog_Pickfile(PATH='..\ephysi\',FILTER = '*.asc') 
      asc=read_ascii(file,DATA_START=4)
      Print,file
      data = asc.(0)
help,data

close,1
openr,1,file
meta=strarr(3)
readf,1,meta

holding = (STRSPLIT(meta[0], /EXTRACT))[1]*1.0 ; define the holding current (pA)
First_Step = (STRSPLIT(meta[2], /EXTRACT))[1]*1.0; define relative difference between injected currents of the first sweep and holding current (pA)
Step_Interval= (STRSPLIT(meta[1], /EXTRACT))[1]*1.0; the step differences between steps(pA)
Sweep_length=1000 ; length of a sweep (ms)
First_Step =  holding + First_Step ;


tau_calculation,data,Holding,First_Step,Step_Interval,Sweep_length,tau=tau
     
sag_calculation,data,sweep_length,Sag_Ratio=Sag_Ratio
      
rebound_fit, data,Sweep_length,rebound = rebound

Input_Resistance,data,Holding,First_Step,Step_Interval,Sweep_length,RMP=RMP,Rinput=Rinput

Action_Potiential_Analysis,data,Holding,First_Step,Step_Interval,Sweep_length,$
        Rheo_Base=Rheo_Base,Spike_Delay=Spike_Delay,Threshold=Threshold,Overshoot=Overshoot,$
        APamp=APamp,AHP=AHP,Depotime=Depotime, Repotime= Repotime,half_width=half_width,RheoBase_Num=RheoBase_Num

Max_Spike_Number,data,holding,First_Step,Step_Interval,Sweep_length,Max_Spike_Num=Max_Spike_Num

Spiking_frequency_adaption, data,Holding,First_Step,Step_Interval,Sweep_length,Rheo_Base,$ 
        AP_Num = AP_Num,Adapt_Index_Initial = Adapt_Index_Initial,Adapt_Index_End = Adapt_Index_End,AP2_Diff = AP2_Diff,$
        AP3_Diff = AP3_Diff, AP_End_Diff = AP_End_Diff

if Max_Spike_Num lt AP_Num then AP_Num = Max_Spike_Num 
if AHP gt 0 then AHP = 0
      DIR=file_DIRname(file,/MARK_DIRECTORY)
      filename=file_basename(file,'.asc')
      file=DIR+filename+'.csv'
      print,file
      print,'#########################  DONE  ###########################'
      openw,file_id,file,/get_lun
      printf,file_id,'Holding(pA),',Holding
      printf,file_id,'Step(pA),',Step_Interval
      Printf,file_id,'First_Current(pA),',First_Step
      printf,file_id,'Tau(ms),',Tau
      Printf,file_id,'Sag_Ratio,',Sag_Ratio
      Printf,file_id,'Rebound(mV),',rebound
      Printf,file_id,'Rm(MO),',Rinput
      Printf,file_id,'RMP(mV),',RMP
      Printf,file_id,'RheoBase(pA),',Rheo_Base
      Printf,file_id,'Spike_Delay(ms),',Spike_Delay
      Printf,file_id,'Threshold(mV),',Threshold
      Printf,file_id,'AP_Amplitude(mV),',APamp
      Printf,file_id,'AHP(mV),',AHP
      Printf,file_id,'Depolarization_Time(ms),',Depotime
      Printf,file_id,'Repolarization_Time(ms),',Repotime
      Printf,file_id,'Half_Width(ms),',Half_width
      Printf,file_id,'Max_Spike_Number,',Max_Spike_Num
      Printf,file_id,'AP_Num@2XRheobase,',AP_Num
      Printf,file_id,'Adapt_Index_Initial,',Adapt_Index_Initial
      Printf,file_id,'Adapt_Index_End,',Adapt_Index_End
      Printf,file_id,'AP2_Diff,',AP2_Diff
      Printf,file_id,'AP3_Diff,',AP3_Diff
      Printf,file_id,'AP_End_Diff,',AP_End_Diff
      free_lun,file_id

END
