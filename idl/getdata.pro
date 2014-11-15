;;; BR notes.
;;
;; 1.  Make stacked image generate a single object (?)  Right now there are a bunch of sources scattered.
;;
;; Later -- 
;; put in seeing (just add sig2_seeing to Gaussians?)
;; make sure profiles are normalized to the same integral over all values of re (I think we need to divide them by r_dev^2 (?)
;; put in extinction.
;; do our values for A and Aerr agree with the values in the catalog?


;; things to add 
; 1.  how do we put in seeing (just add sig2_seeing to Gaussians?)
; 2.  how do we put in extinction?

; Don't forget about Galactic extinction!
; Call this function for 1 object at a time.
; Return Poisson errors for the object also?

;; here's an example CMASS object.
;>>> dd['RA'][0], dd['DEC'][0], dd['R_DEV'][0]
;(129.17618844774728, 48.946499245149603, array([  6.97442301e-05,   4.30699396e+00,   3.84225416e+00,
;         2.94124365e+00,   4.47876596e+00], dtype=float32))


function WHERE_XYZ, Array_expression, Count, XIND=xind, YIND=yind,ZIND=zind
  ; works for 1, 2 or 3 dimensional arrays
  ;
  ; Returns the 1D indices (same as WHERE)
  ;
  ; ARGUMENTS
  ;  - same as WHERE (see WHERE)
  ;
  ; KEYWORDS
  ; - Optionally returns X, Y, Z locations through:
  ;
  ; XIND: Output keyword, array of locations along the first dimension
  ; YIND: Output keyword, array of locations along the second dimension (if present)
  ; ZIND: Output keyword, array of locations along the third dimension (if present)
  ;
  ; If no matches where found, then XIND returns -1
  ;
  index_array=where(Array_expression, Count)
  dims=size(Array_expression,/dim)
  xind=index_array mod dims[0]
  case n_elements(dims) of
    2: yind=index_array / dims[0]
    3: begin
        yind=index_array / dims[0] mod dims[1]
        zind=index_array / dims[0] / dims[1]
      end
    else:
  endcase
  return, index_array
end

function  N2diso, rgrid, sig2inv
        ; Not too sure what imGauss is for ??
        ; Should be working now
        ; n2diso, make_array(5,1),make_array(1,5), 5,5,0.1,1
        ; Can't use this since we don't run getdata anymore 
        ;  n2diso, xgrid,ygrid,xpix[0],ypix[0],0.1,1 

        ;exparg = sig2inv*(xgrid-mx)^2+sig2inv*(ygrid-my)^2
        exparg = sig2inv*rgrid^2
        result= (sig2inv/(2*!pi))*exp(-0.5*exparg)
        return, result
;       stop
 end

; Modular testing works, but can't detect N2diso
;deVauc2d, make_array(5,5),1,1
;pro deVauc2d,xgrid,ygrid,Ie,re
function deVauc2d, rgrid,Ie,re
    ;a = [0.00139,0.00941,0.04441,0.16162,0.48121,1.20357,2.54182,4.46441,6.22820,6.15393] 
    ;v = [0.00087,0.00296,0.00792,0.01902,0.04289,0.09351,0.20168,0.44126,1.01833,2.74555]
    ; Values correct verified with sum a_m =21.2900
    mdev10 = [[0.00139,0.00087],[0.00941,0.00296],[0.04441,0.00792],[0.16162,0.01902],[0.48121,0.04289],[1.20357,0.09351],[2.54182,0.20168],[4.46441,0.44126],[6.22820,1.01833],[6.15393,2.74555]]
    xi = abs(rgrid/re)
    shape = size(rgrid,/dimensions)
    result = make_array(shape[0],shape[1])
    ; use /null so that array unmodified if condition not met
    xx = where_xyz(result GT 8., cnt, XIND=xind, YIND=yind)
    ;stop
    ; Maybe can reuse these from stackimage procedure ? 
    xgrid = intarr(2*shape[0]+1,2*shape[1]+1)
    ygrid = intarr(2*shape[0]+1,2*shape[1]+1)
    for i = 0,9 do $
        result = result + mdev10[0,i]*N2diso(xi,mdev10[1,i])
    if cnt GT 0 THEN $
      result[xind,yind] = 0.
    result = Ie*result
    return, result
 ;stop
end


pro fitimage,imdata,immodel,invvar,A,invvarA,bgd,invvarbgd
  ;; 
  ;; taking equations directly from NR, with immodel replacing x_i.
  s = 0.0
  sx = 0.0
  sy = 0.0
  sxx = 0.0
  sxy = 0.0
  syy = 0.0

  shape = size(imdata,/dimensions)
  shapem = size(immodel,/dimensions)
  if shape[0] ne shapem[0] or shape[1] ne shapem[1] then begin
    print,'shapes dont match in fitimage.  quitting!'
    stop
  endif
  ;; i'm sure there's a vectorized faster way.
  for i = 0,shape[0]-1 do begin
    for j = 0,shape[1]-1 do begin
      s += invvar[i,j]
      sx += immodel[i,j]*invvar[i,j]
      sy += imdata[i,j]*invvar[i,j]
      sxx += immodel[i,j]*immodel[i,j]*invvar[i,j]
      sxy += imdata[i,j]*immodel[i,j]*invvar[i,j]
    endfor
  endfor 
  delta = s*sxx - sx*sx
  bgd = (sxx*sy - sx*sxy)/delta
  invvarbgd = delta/sxx
  A = (s*sxy - sx*sy)/delta
  invvarA = delta/s
  end

;pro stackimage, im, iminvvar  ;; , ra, dec, r_dev
pro stackimage

  nstack = 1

  dd = mrdfits("/home/bareid/imagingsys/cmass-dr12v4-NS-Reid-DORIS.fits",1,hdr)

  ;; this is only correct for things in chunk >= 12.
  setenv,'PHOTO_RESOLVE=/clusterfs/riemann/raid006/bosswork/boss/resolve/2010-05-23'

  ;; chunks 1 <= chunk <= 4
  ;; setenv,'PHOTO_RESOLVE=/clusterfs/riemann/raid006/bosswork/boss/resolve/2009-06-14'
  ;; chunks 5 <= chunk <= 11
  ;; setenv,'PHOTO_RESOLVE=/clusterfs/riemann/raid006/bosswork/boss/resolve/2009-11-16'
  ;; generalize later!


  ;;setenv,'PHOTO_RESOLVE=/clusterfs/riemann/raid006/bosswork/boss/resolve/2013-07-29/'
  ;; I think we should use the targeting PHOTO_RESOLVE results.

  xx = where(dd.ichunk GE 12,cnt)

  ;; let's create the x,y grids only a single time for the maximum value of r_dev.
  ; this is in pixel units; just a test value for now.

  ;; see BR_2dmodel.ipynb; 99% of galaxies are < this R_DEV value.
  r_devmax = [74.9,74.9,14.75,11.97,15.31]

  ;for ifilter = 1,3 do begin 
  for ifilter = 3,3 do begin 
  
    ;; construct small arrays around the point we're interested in.
    ;; different sizes for each band.
  
    npix = floor(r_devmax[ifilter] * 8.) + 1
    xgrid = intarr(2*npix+1,2*npix+1)
    ygrid = intarr(2*npix+1,2*npix+1)
  
  
    xgrid1d = indgen(2*npix+1) - npix 
    ygrid1d = indgen(2*npix+1) - npix 
  
    nx = n_elements(xgrid1d)
    ny = n_elements(ygrid1d)
  
    for i=0,nx-1 do begin
      for j=0,ny-1 do begin
        xgrid[i,j] = xgrid1d[i]
        ygrid[i,j] = ygrid1d[j]
      endfor
    endfor
  
    ;; for a test we're just going to stack images and model and see they match! 
    for indx=0,nstack-1 do begin
      ra = dd[indx].ra
      dec = dd[indx].dec
      r_dev = dd[indx].r_dev[ifilter]
      if r_dev GT r_devmax[ifilter] then begin
        continue ;; skip weirdos.
      end
      ; Read in this image
      infile = sdss_name('idR', dd[indx].run, dd[indx].camcol, dd[indx].field, $
                          rerun=dd[indx].rerun, filter=ifilter)
      sdss_readimage, infile, image, invvar, /calib
      astrans = sdss_astrom(dd[indx].run, dd[indx].camcol, dd[indx].field, $
       rerun=dd[indx].rerun, filter=ifilter)
      astrans_eq2xy, astrans, ra, dec, xpix=xpix, ypix=ypix
      ;; these are arrays, we want scalars for each object.
      xpix = xpix[0]
      ypix = ypix[0]

      ;; which pixel is center closest to (pixcen) and what's the offset from that center (cenoffset).
      pixcen = [floor(xpix+0.5), floor(ypix+0.5)]
      cenoffset = [xpix - floor(xpix+0.5), ypix - floor(ypix+0.5)]

      print,'pixcen',pixcen
      print,'cenoffset',cenoffset

      immin = pixcen - npix
      immax = pixcen + npix
      rgmin = [0,0]
      rgmax = [2*npix,2*npix]
      for ii=0,1 do begin
        if immin[ii] LT 0 then begin
          rgmin[ii] = -immin[ii] ;; relevant rgrid region is smaller.
          immin[ii] = 0
        endif
        if immax[ii] GT n_elements(image[*,0])-1 then begin
          rgmax[ii] = 2*npix -(immax[ii] - n_elements(image[*,0]-1))  ;; relevant rgrid region is smaller.
          immax[ii] = n_elements(image[*,0])-1
        endif
      endfor

      if cenoffset[0] GT 0.501 or cenoffset[0] LT -0.501 or cenoffset[1] GT 0.501 or cenoffset[1] LT -0.501 then begin
        print,'screwed up offsets',xoffset,yoffset
      stop
      end

      print,'why??'
      print,immin
      print,immax
      print,rgmin
      print,rgmax
      print,'end why??'


      rgrid = sqrt((xgrid - cenoffset[0])^2 + (ygrid - cenoffset[1])^2)
      rgrid = rgrid[rgmin[0]:rgmax[0], rgmin[1]:rgmax[1]]

      imcut = image[immin[0]:immax[0],immin[1]:immax[1]]
      invvarcut = image[immin[0]:immax[0],immin[1]:immax[1]]
      t1 = size(rgrid,/dimensions)
      t2 = size(imcut,/dimensions)
      if t1[0] ne t2[0] or t1[1] ne t2[1] then begin
        print,'ack size mismatch',t1,t2
      endif
      immod = deVauc2d(rgrid,1.0,r_dev)
      fitimage,imcut,immod,invvarcut,A,invvarA,bgd,invvarbgd
      print,'results of fit',A,invvarA,bgd,invvarbgd
      print,'blah',dd[indx].modelflux[ifilter], dd[indx].modelflux_ivar[ifilter],A,invvarA
      print,'signal to noise should agree roughly',(dd[indx].modelflux[ifilter])^2*dd[indx].modelflux_ivar[ifilter],A^2*invvarA
      ;; make sure these are the right size.

      xtmp = where_xyz(imcut EQ max(imcut),cnt,XIND=xind,YIND=yind)
      ytmp = where_xyz(immod EQ max(immod),cnt,XIND=xind2,YIND=yind2)
      ztmp = where_xyz(rgrid EQ min(rgrid),cnt,XIND=xind3,YIND=yind3)
      print,'same image centers?'
      print,'mod:',xind2[0],yind2[0],immod[xind2[0],yind2[0]]
      print,'data:',xind[0],yind[0],imcut[xind[0],yind[0]]
      print,'rgrid:',xind3[0],yind3[0],rgrid[xind3[0],yind3[0]]
      ;; ok, looks like peak in imcut is not showing up at the same pixel every time.
      ;; that makes sense; i thought we were fixing that with the imdata offsets.

      ;; for later inspection; just testing!
      mwrfits,imcut,'imdata'+string(indx)+'.fits',/create
      mwrfits,invvarcut,'imdatainvvar'+string(indx)+'.fits',/create
      mwrfits,immod,'immod'+string(indx)+'.fits',/create

    endfor ;end loop over BOSS galaxies.
  endfor ; end loop over ifilter.



stop
end
