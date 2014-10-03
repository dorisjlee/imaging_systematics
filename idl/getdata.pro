; Don't forget about Galactic extinction!
; Call this function for 1 object at a time.
; Return Poisson errors for the object also?

;; here's an example CMASS object.
;>>> dd['RA'][0], dd['DEC'][0], dd['R_DEV'][0]
;(129.17618844774728, 48.946499245149603, array([  6.97442301e-05,   4.30699396e+00,   3.84225416e+00,
;         2.94124365e+00,   4.47876596e+00], dtype=float32))

;function N2d():

function N2diso,xgrid,ygrid,mx,my,sig2inv
  ;; Doris
  return, 1
  end
 
function deVauc2d,xgrid,ygrid,Ie,re
  ;; Doris
  return, 2
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
   for ifilter=3 do begin
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

  xgrid1d = indgen(nx) + xmin
  ygrid1d = indgen(ny) + ymin

  for i=0,nx-1 do begin
    for j=0,ny-1 do begin
      xgrid[i,j] = xgrid1d[i]
      ygrid[i,j] = ygrid1d[j]
    endfor
  endfor



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
