PRO gfunct, X, A, F, pder
;Junzhan Data
  bx = EXP(A[1] * X)
  F = A[0] * bx + A[2]

IF N_PARAMS() GE 4 THEN $
    pder = [[bx], [A[0] * X * bx], [replicate(1.0, N_ELEMENTS(X))]]
END



Pro rebound_fit, data,Sweep_length,rebound = rebound


flag1:

    Episode=N_elements(data(*,0))
    Sample=N_elements(data(0,*))+0.0000
    ;print,sample
    Sample_interval=Sweep_length/Sample
   ; print,Sample_interval
    line_a=reform(data(1,*))
    line_a=smooth(line_a,1/Sample_interval)
window,/free,title='Fit the rebound range to calculate the rebound',xsize=1900,ysize=1000,xpos=0,ypos=0
    plot,line_a
    Base_Start=100
    Plots, Base_Start,[-180,100]
    Base_end=100/Sample_interval-100
    Plots, Base_End,[-180,100]
    rebound_Start=700/Sample_interval
    Plots, rebound_Start,[-180,100]
    
    cursor,xn,yn, /DOWN
    rebound_End=xn
    
    Plots,rebound_End,[-180,100]   
    
   
    Y= data(1,*)
    X= data(0,*)
    Y=Y(rebound_Start:rebound_end)
    X=X(rebound_Start:rebound_end)
    y=y+200
    sz=size(X)
    X=float(indgen(sz[1])) 
    A0=min(y)-max(y)
    A1=-0.001
    A2=max(y)
    
    ;Plot,x,y,yrange=[min(y),max(y)]    
    weights = 1/Y
    A = [A0,A1,A2]
    ;Print,A                   
    yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='gfunct')  
    ;Print,A
    A0=A(0)
    A1=A(1)
    A2=A(2)
    A = [A0,A1,A2]
    yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='gfunct')  
    ;Print,A
    
    oplot,x+rebound_Start,yfit-200,color=255
    ;print,A[2]-200
    
    wait,1
    plot,line_a[17400:rebound_end]
    oplot,x+rebound_Start-17399,yfit-200,color=255
     
    
    
    ;rebound(i)=max(line_a(rebound_Start:rebound_end))-mean(line_a(Base_Start:Base_end)) 
    rebound=A[2]-200-mean(line_a(Base_Start:Base_end)) 
    print,A[2]-200,mean(line_a(Base_Start:Base_end)) ,rebound
    
    cursor,x1,y1,3,/device
    if !mouse.button eq 2 then begin
      wdelete
      goto,flag1;
    endif

    print, 'rebound(mV) ',rebound 


 
cursor,xn,yn, /DOWN 
wdelete

return
End