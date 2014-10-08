;Changes that I made, I will have to manually merge with your new stackimages code later.


; Don't forget about Galactic extinction!
; Call this function for 1 object at a time.
; Return Poisson errors for the object also?

;; here's an example CMASS object.
;>>> dd['RA'][0], dd['DEC'][0], dd['R_DEV'][0]
;(129.17618844774728, 48.946499245149603, array([  6.97442301e-05,   4.30699396e+00,   3.84225416e+00,
;         2.94124365e+00,   4.47876596e+00], dtype=float32))

; make sure to index xpix,ypix value or else the array operation will not work
 function N2diso, xgrid,ygrid,mx,my,sig2inv,imGauss
        ; Should be working now
	; n2diso, make_array(5,1),make_array(1,5), 5,5,0.1,1
        ; Can't use this since we don't run getdata anymore 
	;  n2diso, xgrid,ygrid,xpix[0],ypix[0],0.1,1 
	exparg = sig2inv*(xgrid-mx)^2+sig2inv*(ygrid-my)^2
        result= (sig2inv/(2*!pi))*exp(-0.5*exparg)    
	return, result 
	;stop
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
    xx = where(result GT 8.)
    ; Maybe can reuse these from stackimage procedure ? 
    xgrid = intarr(2*shape[0]+1,2*shape[1]+1)
    ygrid = intarr(2*shape[0]+1,2*shape[1]+1)
    for i = 0,9 do $
	result = result + mdev10[0,i]*N2diso(xi,mdev10[1,i])
    if n_elements(xx) GT 0 THEN $
      result[xx] = 0.
    result = Ie*result
    return, result
 stop
end
 
pro getdata ;; , ra, dec, r_dev

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

  ;; eventually loop over all hte galaxies in dd.
  indx = xx[0]
  ra = dd[indx].ra
  dec = dd[indx].dec
  print,'working on',ra,dec,dd[indx].ichunk

   ; Build a little model of this galaxy image...

   ; Loop over each of the SDSS bands
   ;for ifilter=1,3 do begin
   for ifilter=1,3 do begin
      r_dev = dd[indx].r_dev[ifilter]
      ; Read in this image
      infile = sdss_name('idR', dd[indx].run, dd[indx].camcol, dd[indx].field, $
       rerun=dd[indx].rerun, filter=ifilter)
      print,'beth file:',infile
      sdss_readimage, infile, image, invvar, /calib

      ; Find the center of this object on this image,
      ; assuming it's a de Vauc profile...
      astrans = sdss_astrom(dd[indx].run, dd[indx].camcol, dd[indx].field, $
       rerun=dd[indx].rerun, filter=ifilter)
      astrans_eq2xy, astrans, ra, dec, xpix=xpix, ypix=ypix
      print,'the xy coordinate of this gal is ',xpix,ypix
      ;; next need to figure out how image and invvar are 
    endfor

  ;;for each image, we want to grab relevant subset of pixels
  npix = r_dev * 8.
  
  xmin = floor(xpix - r_dev*8.) - 1
  xmax = floor(xpix + r_dev*8.) + 1

  ymin = floor(ypix - r_dev*8.) - 1
  ymax = floor(ypix + r_dev*8.) + 1

  imcut = image[xmin:xmax,ymin:ymax]
  invvarcut = invvar[xmin:xmax,ymin:ymax]

  nx = n_elements(imcut[*,0])
  ny = n_elements(imcut[0,*])

  xgrid = intarr(nx,ny)
  ygrid = intarr(nx,ny)

  ;xgrid1d = indgen(nx) + xmin
  ;ygrid1d = indgen(ny) + ymin
  ; Be careful of xmin,ymin here, they are array of length 1. 
  ; If no indexing for scalar value, then no elementwise operation occur and the array gets flatten to len = 1
  xgrid1d = make_array(nx,/float)+xmin[0]
  ygrid1d = make_array(ny,/float)+ymin[0]
  for i=0,nx-1 do begin
    for j=0,ny-1 do begin
      xgrid[i,j] = xgrid1d[i]
      ygrid[i,j] = ygrid1d[j]
    endfor
  endfor

stop

  ;; print,n_elements(imcut[*,0]) ;; this is the dimensions of the array.

  ;; testing -- is there a galaxy there?
  atv,image
  atvplot,xpix,ypix,ps=1
  ;; this doesn't work yet! hmmm...

  ; Measure the effective noise from a matched aperture on this
  ; image, the same way that PHOTO does...

  ;; this amounts to summing over invar weighted by our model. just need to figure out how to index invvar.

stop
end
