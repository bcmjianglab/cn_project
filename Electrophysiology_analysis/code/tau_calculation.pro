PRO gfunct, X, A, F, pder
;Junzhan Data
  bx = EXP(A[1] * X)
  F = A[0] * bx + A[2]


  IF N_PARAMS() GE 4 THEN $
    pder = [[bx], [A[0] * X * bx], [replicate(1.0, N_ELEMENTS(X))]]
END

pro Tau_calculation,data,Holding,First_Step,Step_Interval,Sweep_length,tau=tau


Episode=N_elements(data(*,0))
Sample=N_elements(data(0,*))+1.0000
Sample_interval=Sweep_length/Sample


n=Round((holding*1.0-First_Step)/Step_Interval)+1.000
temp=Step_Interval*(n)+First_Step
;if temp gt holding then begin
;  n=n-2
;endif


window,/free,title='Tau Calculation',xsize=1900,ysize=1000,xpos=0,ypos=0
plot,median(reform(data(1,*)), 1/Sample_interval),yrange = [min(data(1,*))-1, max(data(n,*))+1], xrange =[95/Sample_interval,300/Sample_interval]

Tau = fltarr(n-1)
For j=1, n-1 do begin
  TauStart = 101/Sample_interval
  TauEnd = 300/Sample_interval
  DataTau = reform(data(j,*))
  DataTau =median(DataTau, 1/Sample_interval)
  DataTau = DataTau(TauStart:TauEnd)

  TauEnd = where( DataTau eq min(DataTau) )
  TauEnd =TauEnd[0]
  TauEnd = TauEnd+TauStart

  DataTau = reform(data(j,*))
  DataTau =median(DataTau, 1/Sample_interval)
  Y=DataTau
  ;Y= data(j,*)
  oplot, DataTau
  X= data(0,*)
  Y=Y(TauStart:TauEnd)
  X=X(TauStart:TauEnd)
  y=y+300
  sz=size(X)
  X=float(indgen(sz[1]))
  A0=max(y)-min(y)
  A1=-0.002
  A2=min(y)



  ;Plot,x,y,yrange=[min(y),max(y)]
  weights = 1/Y

  A = [A0,A1,A2]
  ; Print,A

  yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='gfunct', STATUS=STATUS, ITMAX = 10000)
  ;print,'STATUS', STATUS, '    Sigma',Sigma



  A0=A(0)
  A1=A(1)
  A2=A(2)
  A = [A0,A1,A2]
  yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='gfunct')

  oplot,x+TauStart,yfit-300,color='0000FF'x
  
  ;print,A[2]-200

 print,'Tau',(-1)/A[1]*Sample_interval
  Tau[j-1]=(-1)/A[1]*Sample_interval
Endfor

wait,0.2
tau=tau(sort(tau))
;tau = tau( tau(where(tau) ne 0))

if N_elements(tau) lt 5 then begin

  finaltau= mean(tau)

endif else begin


  Taubox= CREATEBOXPLOTDATA(tau, MEAN_VALUES=means,  CI_VALUES=CI,OUTLIER_VALUES=outlier,SUSPECTED_OUTLIER_VALUES=SP_outlier)
  print,Taubox
  print,CI
  print,outlier
  print,SP_outlier
  print,means

  if N_elements(outlier) eq 0 and  N_elements(SP_outlier) eq 0 then begin

    FinalTau= mean(tau)

  endif


  if N_elements(outlier) ne 0 and  N_elements(SP_outlier) ne 0 then begin

    FinalTau= [total(tau)-total(outlier)-total(SP_outlier)]/$
      [ N_elements(tau)-0.5*N_elements(outlier)-0.5* N_elements(SP_outlier)]

  endif


  if N_elements(outlier) ne 0 and  N_elements(SP_outlier) eq 0 then begin

    FinalTau= [total(tau)-total(outlier)]/$
      [ N_elements(tau)-0.5*N_elements(outlier)]

  endif

  if N_elements(outlier) eq 0 and  N_elements(SP_outlier) ne 0 then begin

    FinalTau= [total(tau)-total(SP_outlier)]/$
      [ N_elements(tau)-0.5* N_elements(SP_outlier)]

  endif
endelse

print,'Tau(ms)',FinalTau
tau = FinalTau

cursor,xn,yn, /DOWN
wdelete

Return
END