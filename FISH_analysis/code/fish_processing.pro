pro Fish_Processing, Dgene, Threshold,DgeneSubBg = DgeneSubBg, DgeneBgImg = DgeneBgImg
  ;read *.LSM files taken by Zeiss confocal
  ;File_name=Dialog_Pickfile(PATH='Y:\FISH2\Ass1 Vglut2\',FILTER = '*.LSM'); pick up the file
  ;File_name = 'Y:\FISH2\Gad1 GlyT2\Slide1_S2.lsm'

  ;  Print,file_name
  ;  if File_name eq '' then begin
  ;  infor=DIALOG_MESSAGE('Select 1 file!Please run again',Title='Warning!')
  ;  endif else begin
  ;
  ;  img=lsm_read(file_name)
  ;  Dapi = Congrid(REVERSE(img[*,*,0],2),1024,1024) > 0;blue channel,scale data to 1024x1024 pixels
  ;  Cmarker = Congrid(REVERSE(img[*,*,1],2),1024,1024) > 0;common marker gene, green channel
  ;  Dgene = Congrid(REVERSE(img[*,*,1],2),1024,1024)> 0;detected gene,red channel
  print,Threshold
  Dgeneraw = Dgene
  Dgene = EDGE_DOG(Dgene,Threshold = Threshold, RADIUS1=4.0, RADIUS2=7.0,ZERO_CROSSINGS=[0,255])
  lar = label_region(Dgene)
  m=MAX(lar)
  Dgene = Dgeneraw
  Dgene(where(lar))=0
  DgeneBgImg = smooth(Dgene,30)
  DgeneSubBg = Dgeneraw - Dgene > 0
  ;DgeneSubBg =  WATERSHED(DgeneSubBg)
  ;  window,0, title = 'Detecion',xsize = 3072, ysize = 1024,xpos= -1024,ypos=0
  ;  DEVICE, DECOMPOSED=0
  ;  loadct,3
  ;  Tvscl, DgeneSubBg,0
  ;,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  ; DgeneSubBg = Dgeneraw - DgeneBgImg > 0
  ; DgeneSubBgRaw = DgeneSubBg
  ; DgeneSubBg = EDGE_DOG(DgeneSubBg,Threshold = 5.0, RADIUS1=3.0, RADIUS2=5.0,ZERO_CROSSINGS=[0,255])
  ; lar = label_region(DgeneSubBg)
  ; m=MAX(lar)
  ;  DgeneSubBg = DgeneSubBgRaw
  ;  DgeneSubBg(where(lar))=0
  ;  DgeneSubBg = DgeneSubBgRaw - DgeneSubBg > 0
  ;
  ;''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  lar = label_region(DgeneSubBg); label region of gad expression
  m=MAX(lar)
  for i=1, m  do begin
    k=size(lar(where(lar eq i)))
    if k(1) lt 10 then DgeneSubBg(where(lar eq i))=0 ; remove small dots
    if k(1) gt 200 then DgeneSubBg(where(lar eq i))=0 ; remove small dots
  endfor
  ;    loadct,3
  ; Tvscl, DgeneSubBg ,1
  lar = label_region(DgeneSubBg); label region of gad expression
  m=MAX(lar)

  for i=1, m do begin


    k=size(lar(where(lar eq i)))

    indices = where(lar eq i)
    results = Find_Boundary(indices, XSize=1024, YSize=1024, Perimeter=length, Area=area)
    Circularity = 4.0*!pi*area/(length^2)

    if  Circularity lt 0.4 then begin
      DgeneSubBg(where(lar eq i))=0
    endif
    
  endfor
  ;loadct,3
  ;Tvscl, DgeneSubBg,2
  ;
  ; wait, 3
  ; wdelete

  ;endelse
end