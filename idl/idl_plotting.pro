; IDL Version 7.0 (linux x86_64 m64)
; Journal File for dorislee@ln000
; Working directory: /clusterfs/riemann/home/dorislee/imaging_systematics/idl
; Date: Wed Oct 15 19:14:18 2014
 
nx = 128L
ny = 100L
X = FINDGEN(nx)
Y = REPLICATE(1.0, nx)
; Define input function parameters:
aAxis = nx/6.
bAxis = ny/10.
h = 0.5*nx
k = 0.6*ny
tilt = 30*!PI/180
A = [ 5., 10., aAxis, bAxis, h, k, tilt]
 
; Create an ellipse:
xprime = (X - h)*cos(tilt) - (Y - k)*sin(tilt)
yprime = (X - h)*sin(tilt) + (Y - k)*cos(tilt)
U = (xprime/aAxis)^2 + (yprime/bAxis)^2
 
; Create gaussian Z with random noise:
Zideal = A[0] + A[1] * EXP(-U/2)
Z = Zideal + RANDOMN(seed, nx, ny)
B = []  ; clear out the variable
; % Syntax error.
 
; Make about 20% of the points be "bad" data.
bad = WHERE(RANDOMU(1, nx, ny) gt 0.8)
Z[bad] = 999
 
; Create the mask of the bad data points.
mask = REPLICATE(1, nx, ny)
mask[bad] = 0
 
;***** Fit the function *****
yfit = GAUSS2DFIT(Z, B, /TILT, MASK=mask)
; % Keyword MASK not allowed in call to: GAUSS2DFIT
 
; Report results:
PRINT, 'Should be: ', STRING(A, FORMAT='(6f10.4)')
;Should be:      5.0000   10.0000   21.3333   10.0000   64.0000   60.0000     0.5236
PRINT, 'Is: ', STRING(B, FORMAT='(6f10.4)')
; % STRING: Variable is undefined: B.
 
; Create an array with our fitted results
xprime = (X - B[4])*cos(B[6]) - (Y - B[5])*sin(B[6])
; % Variable is undefined: B.
yprime = (X - B[4])*sin(B[6]) + (Y - B[5])*cos(B[6])
; % Variable is undefined: B.
Ufit = (xprime/B[2])^2 + (yprime/B[3])^2
; % Variable is undefined: B.
Zfit = B[0] + B[1] * EXP(-Ufit/2)
; % Variable is undefined: B.
 
; Plot the results. The black dots are missing data.
im = IMAGE(BYTSCL(Z, MAX=20), RGB_TABLE=40, LAYOUT=[1,1,1])
; % Syntax error.
 
; Contour plot of the fit.
c = CONTOUR(Zfit, /OVERPLOT, C_THICK=[4], COLOR='red')
; % Syntax error.
year = [1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006]
landings = [2203, 1905, 1573, 854, 1077, 880, 572, 626, 614]
graph = plot(year, landings, TITLE='Horseshoe Crab Landings in Thousands', $
; % Syntax error.
   XTITLE='Year', YTITLE='Landings in Thousands', $
   thick=3, color='black')
graph = plot(year, landings)
; % Variable is undefined: PLOT.
