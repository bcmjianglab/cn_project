PRO FISH_final_spatial
  ;read *.LSM files taken by Zeiss confocal
  File_name=Dialog_Pickfile(PATH='..t\',FILTER = '*.lsm'); pick up the file
  Print,file_name
  if File_name eq '' then begin
    infor=DIALOG_MESSAGE('Select 1 file!Please run again',Title='Warning!')
  endif else begin
 print,'Processing......'

    filename=file_basename(file_name,'.lsm')
    DIR=file_DIRname(file_name,/MARK_DIRECTORY)
    DIR = DIR +'Intensity_analysis_spatial\'
    FILE_MKDIR,DIR
    img_out = DIR + Filename + '.tif'
    Analysis_out = DIR +Filename + '.pdf'
    data_out= DIR +filename +'.csv'
   
   
   ok = QUERY_TIFF(file_name,s)
if s.NUM_IMAGES gt 2 then begin
    images = bytarr(3,3072,3072)
    IF (ok) THEN BEGIN
      FOR i=0,s.NUM_IMAGES-1 DO BEGIN
        if i mod 2 eq 0 then begin
          img = READ_TIFF(file_name,IMAGE_INDEX=i)
          j = i/2
          x = 0 + 682*round(j/4)
          y = 0 + 682*(j mod 4)
          images[0,y:y+1023,x:x+1023] = img[0, *, *]
          images[1,y:y+1023,x:x+1023] = img[1, *, *]
          images[2,y:y+1023,x:x+1023] = img[2, *, *]
          ;wait,2
          ;tv,img[2, *, *]
        endif
      ENDFOR
    ENDIF
  img = images
endif else begin
  img = read_tiff(file_name)
endelse

    ;img = read_tiff(file_name)
    
    bch = bytarr(3072,3072)
    gch = bytarr(3072,3072)
    rch = bytarr(3072,3072)

    for i = 0, 3071 do begin
      for j = 0 , 3071 do begin
        bch(i,j) = img(0,i,j);;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        gch(i,j) = img(1,i,j)
        rch(i,j) = img(2,i,j);;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      endfor
    endfor

    
    Dapi = Congrid(REVERSE(smooth(bch,4),2),1024,1024) > 0;blue channel,scale data to 1024x1024 pixels
    GreenCh = Congrid(REVERSE(smooth(gch,4),2),1024,1024) > 0;common marker gene, green channel
    RedCh = Congrid(REVERSE(smooth(rch,4),2),1024,1024) > 0;detected gene,red channel
;   

RedCh = BYTSCl(RedCh,max = 150)
GreenCh =  BYTSCl(GreenCh,max = 200)
;Dapi =  BYTSCl(Dapi,max = 255)
Threshold = 15

img = bytarr(3,1024,1024)
img[1,*,*] = GreenCh;G
img[2,*,*] = Dapi;B
img[0,*,*] = RedCh;R

    
    RedChRaw = RedCh
    GreenChRaw = GreenCh
    merged = FLTARR(1024, 1024)

    for i = 0 , 1023  do begin
      for j = 0, 1023 do begin
        if RedCh(i,j) ge GreenCh(i,j) then begin
          merged(i,j) = RedCh(i,j)*1.000
        endif else begin
          merged(i,j) = GreenCh(i,j)*1.000
        endelse
      endfor
    endfor
    MergeRaw = Merged
    print,'Processing Merged'
    FISH_Processing, Merged,Threshold, DgeneSubBg = Subbg, DgeneBgImg = bg
    AllCell = Subbg

    print,'Processing Red Channel'
    FISH_Processing, RedCh,Threshold,DgeneSubBg = RedChSubbg, DgeneBgImg = RedChbg
    RedChSubbg = RedChSubbg
    RedChbg = RedChbg
    
    print,'Processing Green Channel'
    FISH_Processing, GreenCh,Threshold,DgeneSubBg = GreenChSubbg, DgeneBgImg = GreenChbg
    GreenChSubbg = GreenChSubbg
    GreenChbg = GreenChbg

    xroi, img,img[0,*,*],img[1,*,*],img[2,*,*],Regions_Out=thisROI, /Block, TOOLS= 'Polygon Draw'
    mask = thisROI -> ComputeMask(Dimensions=Size(Dapi, /Dimensions), Mask_Rule=2)
    Result = thisROI -> ComputeGeometry(CENTROID=CT)
    CTx = ct(0)
    CTy = ct(1)
    maxdis = 0
    mindis = 1000.0

    for i = 0, 1023 do begin
      if i le CTx then begin
      for j = 0, 1023 do begin
       ;if j lt CTy then begin;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        if mask(i,j) eq 255 then begin
          dis = sqrt((i-CTx)*(i-CTx)+(j-CTy)*(j-CTy))
          if dis ge maxdis then begin
            maxdis = dis
            maxX = i
            maxY = j
          endif
        endif
       ;endif;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      endfor
      endif
    endfor

    for i = 0, 1023 do begin
      for j = 0, 1023 do begin
        if mask(i,j) eq 0 then begin
          dis = sqrt((i-CTx)*(i-CTx)+(j-CTy)*(j-CTy))
          if dis le mindis then begin
            mindis = dis
            minX = i
            minY = j
          endif
        endif
      endfor
    endfor
    ;print,CTx,CTy,maxX,maxy,minX,miny
    window, xsize = 1024, ysize = 1024
    DEVICE, DECOMPOSED=0
    loadct,7
    tv,MergeRaw
    plots,ctx,cty,Psym=1,  SYMSIZE=2, /DEVICE
    plots,maxX,maxy,Psym=1,  SYMSIZE=2,/DEVICE
    WAIT,1
    
    DetaX=maxX-ctx
    DetaY=maxy-cty
    If DetaX EQ 0 then begin
      If DetaY GT 0 then begin
        Ang=0
      Endif else begin
        Ang=180
      Endelse
    Endif else begin
      Ang=90.0*abs(DetaX)/(abs(DetaX)+abs(DetaY))
      If DetaX GT 0 then begin
        If DetaY GE 0 then begin
          Ang=360-Ang
        Endif else begin
          Ang=180+Ang
        Endelse
      Endif else begin
        If DetaY GE 0 then begin
          Ang=Ang
        Endif else begin
          Ang=180-Ang
        Endelse
      Endelse
    Endelse
    print,'Rotate Angle ',ang

ctx =511
cty = 511
;ang =0

    temp = RedChSubbg
    outside = RedChSubbg - mask > 0; outside of CN
    RedChSubbg = temp - outside; keep CN area only

    temp = GreenChSubbg
    outside = GreenChSubbg - mask > 0; outside of CN
    GreenChSubbg = temp - outside; keep CN area only

    temp = AllCell
    outside = AllCell - mask > 0; outside of CN
    AllCell = temp - outside; keep CN area only



        if ang lt 90 then begin
          AllCell =  rot(AllCell,Ang,1,ctx,cty,missing=0)
          RedChRaw =  rot(RedChRaw,Ang,1,ctx,cty,missing=0)
          RedChbg =  rot(RedChbg,Ang,1,ctx,cty,missing=0)
          GreenChRaw =  rot(GreenChRaw,Ang,1,ctx,cty,missing=0)
          GreenChbg =  rot(GreenChbg,Ang,1,ctx,cty,missing=0)
          Dapi =  rot(Dapi,Ang,1,ctx,cty,missing=0)
          
          MergeRaw = rot(MergeRaw,Ang,1,ctx,cty,missing=0)
        endif else begin
          if ang lt 180 then begin
          AllCell =  REVERSE(rot(AllCell,Ang,1,ctx,cty,missing=0))
          RedChRaw =  REVERSE(rot(RedChRaw,Ang,1,ctx,cty,missing=0))
          RedChbg =  REVERSE(rot(RedChbg,Ang,1,ctx,cty,missing=0))       
          GreenChRaw =  REVERSE(rot(GreenChRaw,Ang,1,ctx,cty,missing=0))
          GreenChbg =  REVERSE(rot(GreenChbg,Ang,1,ctx,cty,missing=0))
          Dapi =  REVERSE(rot(Dapi,Ang,1,ctx,cty,missing=0))
          MergeRaw =  REVERSE(rot(MergeRaw,Ang,1,ctx,cty,missing=0))
          endif else begin
            print,'No rotation'
            goto,flag1
          endelse
        endelse
    tv,MergeRaw
    wait,1
    wdelete
    larA = label_region(AllCell)
    An=MAX(larA)
    CellR = fltarr(An)
    CellG = fltarr(An)
    CellRb = fltarr(An)
    CellGb = fltarr(An)

    for i = 1, An do begin
      k=size(larA(where(larA eq i)))
      CellR(i-1) = mean(RedChRaw(where(larA eq i)))
      CellRb(i-1) = mean(RedChbg(where(larA eq i)))
      CellG(i-1) = mean(GreenChRaw(where(larA eq i)))
      CellGb(i-1) = mean(GreenChbg(where(larA eq i)))
    endfor
    cellRb(where(cellRb eq 0)) = mean(cellRb)
    cellGb(where(cellGb eq 0)) = mean(cellGb)

    x1 = (cellR-CellRb)/CellRb
    sc1 = max(x1(where(x1 ne max(x1))))*1.0
    x1 = x1/max(x1(where(x1 ne max(x1))))*10.0
    
   
    

    y1 = (cellG-CellGb)/CellGb
    sc2 = max(y1(where(y1 ne max(y1))))*1.0
    y1 = y1/max(y1(where(y1 ne max(y1))))*10.0


    ;print,sc1,sc2

    p = plot(x1,y1,linestyle=6,symbol = "Circle",SYM_SIZE = 0.5,SYM_FILLED = 1,SYM_COLOR = 'Red', SYM_FILL_COLOR = 'Red',AXIS_STYLE = 1, xrange =  [0.0001,10],yrange = [0.0001,10]);, xlog = 1,ylog = 1)
    p.save,Analysis_out

    Raw = AllCell
    openw,file_id,data_out,WIDTH=60000,/get_lun
    printf,file_id, ' Cell_ID' ,' Red', ' Red_Bg', ' Green',' Green_Bg', ' Scaled_Red', ' Scaled_Green',' AnteriortoPostX', ' VentraltoDorsalY',' Area';,xcutoff, ycutoff

    for i = 1, An do begin
      k=size(larA(where(larA eq i)))
      CellR = mean(RedChRaw(where(larA eq i)))
      CellRb = mean(RedChbg(where(larA eq i)))
      CellG = mean(GreenChRaw(where(larA eq i)))
      CellGb = mean(GreenChbg(where(larA eq i)))
      RC = (cellR-CellRb)/CellRb & RC = RC/sc1*10.0
      GC = (cellG-CellGb)/CellGb & GC = GC/sc2*10.0
      
      
      rrr = (RedChRaw(where(larA eq i)) - CellRb)/CellRb
      ggg = (GreenChRaw(where(larA eq i)) - CellGb)/CellGb
      r = 0
      p = 1
      if n_elements(rrr) ge 2 then begin
      c = R_CORRELATE(rrr,ggg, /KENDALL)
      r = c[0]
      p = c[1]
      endif
      pos=median(where(larA eq i))
      X= pos-fix(pos/1024)*1024l-511.5
      Y= fix(pos/1024)
      

        printf,file_id, i , CellR, CellRb, CellG, CellGb,RC,GC,x,y,k(1), r,p
        AllCell(where(larA eq i)) = 255



    endfor
    wait,1
    p.close
    free_lun,file_id
    
AllCell(where(AllCell lt 255)) = 0

img = bytarr(3,1024,1024)
img[1,*,*] = GreenChRaw;G
img[2,*,*] = AllCell;B
img[0,*,*] = RedChRaw;R
img = image(img, DIMENSIONS=[512, 512])
wait,2
img.save,img_out,RESOLUTION=300,BORDER=0, /TRANSPARENT
img.close

  result = FILE_TEST(DIR + 'fish_final_spatial.pro')
  if result eq 1 then file_delete,DIR + 'fish_final_spatial.pro'

endelse 
flag1:


Print,'DONE.....................................................Please selet next image '

END
