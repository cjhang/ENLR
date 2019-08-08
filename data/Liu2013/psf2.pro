pro psf2
;make PSF comparison figure

obj=['J0149-0048','J0210-1001','J0319-0019','J0319-0058', $
		 'J0321+0016','J0759+1339','J0841+2042','J0842+3625', $
		 'J0858+4417','J1039+4512','J1040+4745', $
		 'J0224+2750','J0807+4946','J1101+4004']
psf=[3.4,4.0,2.5,3.9,3.7,4.1,3.25,3.15,4.0,3.7,5.0,3.7,3.9,4.3]
;peak=[1286.01,	2489.36,1051.09,1068.59, 745.516,3160.18,5312.75,4776.56,2097.67,4153.9,2594.53,1951.23,1000.]
;peak=['1.29','2.49','1.05','1.07','0.75','3.16','5.31','4.78','2.10','4.15','2.59','2.40','2.11','1.95']
peak=['1.58','3.19','1.37','1.38','0.96','3.93','6.72','6.13','2.64','5.20','3.31','2.40','2.66','2.36']
eta=['5.70','3.00','3.00','4.30','3.02','3.39','3.84','3.93','4.06','4.00','3.53','3.38','3.49','4.64']
err=['0.14','0.02','0.20','0.11','0.08','0.02','0.02','0.02','0.02','0.02','0.01','0.06','0.19','0.12']
R_5sig=[0.970425,0.757695,1.44483,0.807992,1.19048,0.867303,0.799419,1.00604,1.19771,0.746951,1.15230,1.14816,0.886986,0.812586]
	fits_read,'/home/liu/idl/flux/orig/star_ell.fits',tab,htab,exten=1
    	sma0 = tbget(htab,tab,'sma')*1d
	intens0 = tbget(htab,tab,'intens')*1d

cdelt=0.145779

	set_plot,'ps'                 ;plot on PS file
	device,filename='psf2.ps',xsize=10.0,ysize=12,/inches;,xoffset=1.5,yoffset=1
	!p.thick=3
	!x.margin=[10,5]
	!y.margin=[5,3]
	!x.thick=2
	!y.thick=2
  	!p.charthick=2

	multiplot,[2,7], gap=0. ;, mXtitle='xx',/square $
			;mYtitle='Normalized Surface Brightness'
	for i=0,13 do begin
		fits_read,'/home/liu/idl/flux/orig/'+obj[i]+'ell.fits',tab,htab,exten=1
        sma = tbget(htab,tab,'sma')*1d
        intens = tbget(htab,tab,'intens')*1d
        ellip = tbget(htab,tab,'ellip')*1d
	    ellip[where(finite(ellip) ne 1)]=0.
	  sig=psf[i]*cdelt/sqrt(8*alog(2))
	  tmp=findgen(1000)/100.
	  if (i eq 13) then xrange=[0,2.0] else xrange=[0,2.0-1e-5]
	  if (i eq 12) then yrange=[1e-2,1] else yrange=[1.001e-2,1]
		plot,tmp,exp(-(tmp/sig)^2/2.),linestyle=1,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,/ylog,/nodata,charsize=1.3
		oplot, sma0*0.1*psf[i]*cdelt/0.6,intens0/max(intens0),linestyle=1
;		oplot, sma0*0.1,intens0/max(intens0),linestyle=1
		oplot, sma*0.1,intens/max(intens)
;		oplot, 0.1*smb[sort(smb)],intens/max(intens),linestyle=3
		oplot, 0.1*sma*(1-ellip),intens/max(intens),linestyle=3
	  if (obj[i] ne 'J0224+2750') then xyouts,0.85,0.48,'SDSS '+obj[i],charsize=1.4 else $
		xyouts,0.85,0.48,'3C67/'+obj[i],charsize=1.4
;		peak[i]=max(intens[0:3])
		xyouts,0.07,0.028,textoidl('I_0=')+peak[i],charsize=1.3
		xyouts,0.07,0.015,textoidl('\eta=')+eta[i]+textoidl('\pm')+err[i],charsize=1.3
		
		fits_read,'im_cont/p'+obj[i]+'ct.fits',tab,htab,exten=1
	        sma = tbget(htab,tab,'sma')*1d
        	intens = tbget(htab,tab,'intens')*1d
        	ellip = tbget(htab,tab,'ellip')*1d
		ellip[where(finite(ellip) ne 1)]=0.
		
		print,obj[i]

		use=where(sma*0.1 lt R_5sig[i])
		oplot, sma[use]*0.1,intens[use]/max(intens),linestyle=2
	
;;;		oplot, 0.1*sma*(1-ellip),intens/max(intens),linestyle=4

	  multiplot
		print,obj[i],max(intens)
endfor
xyouts,0.33,.037,'R (arcsec)',font=-1, charsize=1.5, charthick=3, /normal, align=0.5
xyouts,0.73,.037,'R (arcsec)',font=-1, charsize=1.5, charthick=3, /normal, align=0.5
xyouts,0.07,0.45,textoidl('log [I_R / I_0]'),orientation=90, charsize=1.5, charthick=3,/normal

multiplot,/reset

!p.multi=[0,1,1]               ;go back to Xwindow
device,/close_file
set_plot,'X'

spawn,'gv psf2.ps &'
spawn,'cp psf2.ps ~/Desktop/AAStex/'

end
