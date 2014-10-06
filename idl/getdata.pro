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

;function N2d():

;pro N2diso,xgrid,ygrid,mx,my,sig2inv,imGauss
;  exparg = (xgrid - mx)*(xgrid - mx)
;  return, 1
;  end
 
;function deVauc2d,xgrid,ygrid,Ie,re
  ;; Doris
;  return, 2
;  end

;pro fitimage,imdata,immodel,A,B
;  A = 5.
;  B = 7.


;pro stackimage, im, iminvvar  ;; , ra, dec, r_dev
pro stackimage, imdata 

  nstack = 2

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
    rgrid = intarr(2*npix+1,2*npix+1) 
  
    imdata = fltarr(2*npix+1,2*npix+1)
    immodel = fltarr(2*npix+1,2*npix+1)

    imdata[*,*] = 0.0d0
  
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
      print,indx,xpix,ypix


      xoffset = xpix - floor(xpix+0.5)
      yoffset = ypix - floor(ypix+0.5)

      xpixcen = floor(xpix+0.5)
      ypixcen = floor(ypix+0.5)
      print, 'yo',xoffset, yoffset,xpixcen,ypixcen

      xmin = max([xpixcen - npix,0])
      ymin = max([ypixcen - npix,0])
      xmax = min([xpixcen + npix,n_elements(image[*,0])-1])
      ymax = min([ypixcen + npix,n_elements(image[0,*])-1])

      if xoffset GT 0.501 or xoffset LT -0.501 or yoffset GT 0.501 or yoffset LT -0.501 then begin
        print,'screwed up offsets',xoffset,yoffset
      stop
      end

      rgrid = sqrt((xgrid - xoffset)^2 + (ygrid - yoffset)^2)
      ;; now send this to the model; we should make hte models accept rgrid instead of xgrid,ygrid?
      ;; now get the data stack
      print,xmin,ymin,xmax,ymax
      print,xpixcen,ypixcen,n_eleements(imcut[*,0]),
      imcut = image[xmin:xmax,ymin:ymax]
      invvarcut = invvar[xmin:xmax,ymin:ymax]
      
      ;; make sure these are the right size.
      print,'same?',n_elements(imcut[*,0]), n_elements(imcut[0,*]), n_elements(imdata[*,0]), n_elements(imdata[0,*])

      xdmin = 0
      xdmax = 2*npix
      if xpixcen - npix LT 0 then begin
        xdmin = npix - xpixcen
      end
      if xpixcen + npix GT n_elements(image[*,0]) then begin
        xdmax = 2*npix - (n_elements(image[*,0]) - xpixcen - npix)  
      end

      ydmin = 0
      ydmax = 2*npix
      if ypixcen - npix LT 0 then begin
        ydmin = npix - ypixcen
      end
      if ypixcen + npix GT n_elements(image[0,*]) then begin
        ydmax = 2*npix - (n_elements(image[0,*]) - ypixcen - npix)  
      end

      print,'this should be',xmax-xmin, xdmax-xdmin, ymax-ymin, ydmax-ydmin
      print,'beth this is broken!'
      ;; this is allowed even when the arrays aren't hte same size!  ackk!!!
      imdata += imcut
      mwrfits,imdata,'imstacktest.fits',/create
    endfor ;end loop over BOSS galaxies.
  endfor ; end loop over ifilter.


stop
end
