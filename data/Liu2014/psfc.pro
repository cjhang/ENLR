pro psfc,lo3=lo3,loiii=loiii,sdsserr=sdsserr,calib=calib
;make PSF comparison figure

ctrx=34.08-1 & ctry=25.51-1
obj=['J0233','J0304','J0311','J0412','J0753','J0809','J0847','J0909','J0924','J0935','J1144','J2214']
dec=['-0743','+0022','-0707','-0511','+3153','+0743','+2940','+3459','+0642','+5348','+1043','+2115']
psf=[2.9,3.2,3.95,2.7,3.5,3.3,4.1,4.3,3.3,4.8,3.7,3.1]
;psf=[2.9,3.2,4.1,2.7,2.7,3.3,4.1,4.3,3.3,4.8,3.7,3.1]
readcol,'~/iraf/T1.dat',a,a,a,a,a,a,a,kpc,format='a,a,a,a,a,a,a,f'
Rct_5sig=[8.5,12.3,13.2,9.9,7.8,12.4,8.3,11.7,11.5,11.5,10.3,9.7]
Rmax_ct=Rct_5sig/kpc
R_5sig=[16.6,10.9,11.,13.6,8.6,12.6,12.,13.3,9.1,13.6,14.1,10.7]
Rmax=R_5sig/kpc
pk=[30.15, 8.121, 13.98, 239.9, 23.34, 14.7, 10.16, 9.501, 16.34, 16.89, 20.92, 24.37]

;;;;;;;;;;;;; Flux Cal using SDSS ;;;;;;;;;;;;;

flx=[2216.,363.,704.,6671.,1456.,1021.,917.,1009.,744.,2044.,2398.,990.]
pk=[29.58, 8.544, 11.54, 193.1, 23.6, 14.17, 8.445, 9.021, 15.71, 16.01, 20.18, 18.56] ;1-gauss
pk=[30.15, 8.121, 13.98, 239.9, 23.34, 14.7, 10.16, 9.501, 16.34, 16.89, 20.92, 24.37]
readcol,'~/iraf/T1.dat',fname,radec,raa,decc,zapp,Lo3,e_Lo3,kpc,Dl,xx,Rint,format=('a,a,a,a,f,f,f,f,f,f,f')
sdsscal=zapp*0.
sdsserr=zapp*0.

for i=0,11 do begin
  readcol,'~/iraf/resp_'+fname[i]+'.dat',wv,response,format=('F,F')
  readcol,'~/iraf/sdss_simu_'+fname[i]+'.dat',wv,spec,format=('F,F')
  readcol,'~/iraf/sdss/int_'+fname[i]+'.dat',x0,y0,invar,format=('F,F,F')
  x0=x0/(1+6.4328D-5+2.94981D-2/(146-1D8/x0^2d)+2.554D-4/(41.-1D8/x0^2d))/(1.+zapp[i])
  quadterp,x0,y0,wv,ycal
  quadterp,x0,invar,wv,invarcal
  cal=[4980,5050]
  use=[where(wv ge cal[0] and wv le cal[1])]
  sdsscal[i]=total(ycal[use])/total(spec[use]*response[use])
  use0=[where(x0 ge cal[0] and x0 le cal[1])]
  use0b=[where(x0 ge 4985 and x0 le 5050)]
;  flx0=10d^Lo3[i]/(1d-17*4.*!dpi*(Dl[i]*3.08568d24)^2)/(2926./1873.)
; print,sqrt(total((invar[use0])^2))/total(y0[use0])*100
  uncer1=sqrt(total((invarcal[use0])^2))/total(ycal[use0])
  uncer2=sqrt(total((1/invarcal[use0])^2))/total(ycal[use0])
  sdsserr[i]=max([uncer1,uncer2])

  print,'Object:  ',obj[i]
  print,'calibration:  ',sdsscal[i]
  print,'uncertainty:  ',strtrim(string(sdsserr[i]*100),2),'%'

endfor

aper_uncorr=32.04/41.6
pk*=sdsscal/299792.458*5006.7*sqrt(2*!dpi)/0.1^2*1e-3*aper_uncorr
calib=sdsscal/299792.458*5006.7*sqrt(2*!dpi)/0.1^2*1e-3*aper_uncorr
peak=string(pk,format='(f4.2)')
loiii=alog10(flx*sdsscal*aper_uncorr/299792.458*5006.7*sqrt(2*!dpi)*1d-17*4.*!dpi*(Dl*3.08568d24)^2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

eta=obj & err=obj

s1=[0.9,1.2,1.2,1.3,1.1,1.3,1.3,1.3,1.2,1.1,1.3,1.0]
s2=[1.5,2.0,2.0,2.5,2.5,2.2,2.5,2.5,2.0,2.5,2.5,2.5]

for i=0,11 do begin
  fits_read,'~/iraf/dat/'+obj[i]+'/o3lp.fits',tab,htab,exten=1
  x = alog10(tbget(htab,tab,'sma')*0.1)
  y = alog10(tbget(htab,tab,'intens'))
  use=where(x ge alog10(s1[i]) and x le alog10(s2[i]))
  sixlin,x[use],y[use],a,siga,b,sigb
  plot,x,y,title=obj[i],xrange=[-0.5,1]
  oplot,x[use],x[use]*b[0]+a[0]
  eta[i]=strtrim(string(-b[0],format='(f4.2)'),2) 
  err[i]=strtrim(string(sigb[0],format='(f4.2)'),2) 
;  wait,2
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	star=readfits('~/iraf/psf/star.fits',h)-2e-4
;        star[28:42,0:11]=!values.f_nan & star[27:40,39:48]=!values.f_nan
;        star[28:42,0:11]=-10 & star[27:40,39:48]=-10
        nx=(size(star))[1]  &  ny=(size(star))[2]
        dist_circle,dist,[nx,ny],ctrx,ctry,/double

cdelt=0.1456

	set_plot,'ps'                 ;plot on PS file
	device,filename='psfc.ps',xsize=11.0,ysize=10,/inches;,xoffset=1.5,yoffset=1
	!p.thick=3
	!x.margin=[8,2]
	!y.margin=[3,1]
	!x.thick=2
	!y.thick=2
  	!p.charthick=2

	multiplot,[3,4], gap=0. ;, mXtitle='xx',/square $
			;mYtitle='Normalized Surface Brightness'
for i=0,11 do begin

   line=readfits('~/iraf/dat/'+obj[i]+'/nonpar/flx.fits',h)
   cont=readfits('~/iraf/psf/cont/'+obj[i]+'ct.fits',h)

   tmp=findgen(1000)/100.
   bin=0.1
   x=findgen(round(2.5/bin))*bin+0.5*bin
   y0=x*0.
   yl=y0 & yc=y0
   for j=0,n_elements(x)-1 do begin
     use=where(dist*0.1 ge j*bin and dist*0.1 le (j+1)*bin)
;     y0[j]=median(star[use])
     RESISTANT_Mean,star[use], 1, mean, meansig, num & y0[j]=mean
;     y0[j]=median(star[use])
     RESISTANT_Mean,line[use], 2, mean, meansig, num & yl[j]=mean
;     y0[j]=median(star[use])
     RESISTANT_Mean,cont[use], 2, mean, meansig, num & yc[j]=mean
   endfor
   if (i eq 9 or i eq 10) then xrange=[0,2.4999] else xrange=[0,2.5]           
   plot,x*psf[i]*cdelt/0.508,y0,xrange=xrange, yrange=[0.001,3],/ylog,/xstyle,/ystyle,linestyle=1
   goodc=where(x lt Rmax_ct[i])
   goodl=where(x lt Rmax[i])
   oplot,x[goodl],yl[goodl]/max(yl),linestyle=0
   oplot,x[goodc],yc[goodc]/max(yc),linestyle=2

   xyouts,0.65,0.97,'SDSS '+obj[i]+dec[i],charsize=1.3
   xyouts,0.10,0.0043,textoidl('I_0=')+peak[i],charsize=1.2
   xyouts,0.10,0.0019,textoidl('\eta=')+eta[i]+textoidl('\pm')+err[i],charsize=1.2

   oplot,[1,1]*Rint[i]/kpc[i],[1e-4,0.35],linestyle=3

   multiplot
endfor

xyouts,0.23,.005,'R (arcsec)',font=-1, charsize=1.5, charthick=3, /normal, align=0.5
xyouts,0.53,.005,'R (arcsec)',font=-1, charsize=1.5, charthick=3, /normal, align=0.5
xyouts,0.83,.005,'R (arcsec)',font=-1, charsize=1.5, charthick=3, /normal, align=0.5
xyouts,0.024,0.45,textoidl('log [I_R / I_0]'),orientation=90, charsize=1.6, charthick=3,/normal

multiplot,/reset

!p.multi=[0,1,1]               ;go back to Xwindow
device,/close_file
set_plot,'X'

spawn,'gv psfc.ps &'
spawn,'cp psfc.ps ~/Desktop/AAStex/'

end
