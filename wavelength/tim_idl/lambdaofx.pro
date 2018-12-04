pro lambdaofx,xx,mm,d,sinalp,fl,a0,y0,r0,gltype,priswedge,lamcen,lam
; This routine computes wavelength lam(nm) as a function of x-coordinate xx
; and order number mm, for a cross-dispersed echelle.
; On input,
;  xx(nx) = x-coordinate vector, measured from CCD center, in mm
;  mm(nord) = vector containing order numbers for which lam is to be calculated
;  d = grating groove spacing (mm)
;  sinalp = nominal sin of incidence angle on grating
;  a0 = correction to sinalp, resulting in x-coord displacement of spectrum
;     on chip, relative to nominal position.
;  fl = focal length of the spectrograph camera (mm)
;  y0 = y-coordinate on CCD (mm) at which gamma angle = 0.
;  r0 = rotation angle of CCD (radian) relative to main dispersion along x
;  gltype = prism glass type
;  priswedge = prism vertex angle (degrees)
;  lamcen = nominal wavelength (micron) for which prism out-and-back 
;           deviation is zero
;
; On return,
;  lam(nx,nord) = computed wavelength vs x,order
;
; Method is to combine the diffraction equation with geometrical relations
; for image scale and rotation.

; constants
!except=2
radian=180d0/!dpi

; get sizes of things
nx=n_elements(xx)
nord=n_elements(mm)

; make output arrays, some working arrays
y0m=dblarr(nord)       ; y(mm), ordered the same as mm
xo=dblarr(nx,nord)     ; x-coordinate in SG natural coords, before det rotation
lam=dblarr(nx,nord)
dlamda0=dblarr(nx,nord)
dlamdfl=dblarr(nx,nord)
dlamdy0=dblarr(nx,nord)
dlamdr0=dblarr(nx,nord)
mmo=rebin(reform(mm,1,nord),nx,nord)

; fill y0m
; compute prism half-deflection angle for center wavelength
glass_index,gltype,[lamcen],nn
nnmid=nn[0]
delta=asin(nnmid*sin(priswedge/(2d0*radian))) - priswedge/(2d0*radian)
; compute central wavelengths for each order
dd=d*1000.                   ; groove spacing (microns)
lamc=dd*2.*sinalp/mm         ; central wavelength for each order (microns)

; get refractive indices, center y-coordinates for each order
glass_index,gltype,lamc,nnc
y0m=4d0*fl*delta*(nnc-nnmid)    ; 2 passes through full deflection angle

; compute gamma angle, sines and cosines of same
gammar=asin((y0m-y0)/fl)
cosgam=cos(gammar)
singam=sin(gammar)
cosgamo=rebin(reform(cosgam,1,nord),nx,nord)
singamo=rebin(reform(singam,1,nord),nx,nord)

; compute unrotated x coordinates
for i=0,nord-1 do xo[*,i]=xx+y0m[i]*sin(r0)

; compute sin(beta)
;sinbet=xo/fl + sinalp + a0
alp=asin(sinalp + a0)
bet=alp+atan(xo/fl)
sinbet=sin(bet)
sinalpo=sin(alp)

; compute lambda (mm)
lam=(d/mmo)*(sinalpo+sinbet)*cosgamo

; convert lam to nm
lam=lam*1.d6

end
