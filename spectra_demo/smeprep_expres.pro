function smeprep_expres, fitsSpectrum, outFile=outFile, segPath=segPath, segFiles=segFiles, $
					freeElems=freeElems, vmacPro=vmacPro, atmogrid_file=atmogrid_file, $
					verbose=verbose, stop_on_err=stop_on_err
		
	;
	; EXPRES Spectra (reduced) are fits files containing the original raw header
	;	plus the reduced spectrum, continuum normalized spectrum, associated wavelengths
	;	both pre/post barycentric correction and uncertainties, along with additional
	;	reduction information all in a FITS Binary table.
	;
	; The function returns a structure containing 'status' and 'file'
	;	where status: 0 is a succesful completion and 'file' will contain the generated .inp file.
	;	If status != 0, then file may not exist and 'msg' will be present explaining the
	;	reason for the failure.
	;
	
	; Constants
	c_light = double(2.9979246d5) ; km/s
	solar_abund, solab, elems
	expres_ipres = 137500 ; will look in header for RESOLUTN and if it isn't found use this
	
	; -----------------------------------------------------------------
	; Set up all of the parameter defaults for keyword arguments
	;
	IF ~ keyword_set(segPath) THEN segPath='/home/bizard/astronomy/projects/spocs/segf/'
	IF ~ keyword_set(segFiles) THEN BEGIN
		segFiles = segpath + [ $
				'5164_5190c', $		; SPOCS I: RE-retune with clean vald3 lines
				'5190_5207', $		; Ca, Cr, Y
				'5232_5262', $		; Nd, Mn, Cr, V, Sc, Co, Ca
				'6000_6015', $		; SPOCS I - All SPOCS I segments have been retuned against
				'6015_6030', $		;	...		Wallace 2011 atlas and contains lines from VALD 3
				'6030_6050', $		;  ...		
				'6050_6070', $		;  ...
				'6100_6120no_li', $ ;  ...
				'6121_6140b', $		;  ...
				'6143_6160', $		;  ...
				'6160_6180', $		;  ...
				'6295_6305', $		; O, Mg
				'6311_6320', $		; O, Mg - gap from prior segment as it crossed orders at HIRES
				'6579_6599no_li', $ ; C, plus Cr, Ti, Ni, Si
				'6688_6702', $		; Al, Ca
				'6703_6711', $		; Li, CN, Ca
				'6711_6718', $		; Ca, Cr, Zr (1 or two of each)
				'7440_7470b', $		; N plus S
				'7697_7702', $		; K (nlte)
				'7769_7800b' ] $ 	; O triplet plus a few Ca and 1 faint Nd
				+ '_jmbf.out'
	ENDIF
	IF ~ keyword_set(freeElems) THEN freeElems = ['C','N','O','Na','Mg','Al','Si','Ca','Ti','V','Cr','Mn','Fe','Ni','Y']
	IF ~ keyword_set(atmogrid_file) THEN atmogrid_file = 'atlas12.sav' ; 'marcs2012.sav' ; normal is atlas12.sav (kurucz)
	IF n_elements(vmacPro) THEN vmac_pro = vmacPro ELSE vmac_pro = ''
	IF ~ keyword_set(verbose) THEN verbose = 0
	IF ~ keyword_set(outFile) THEN BEGIN
		outfile = (stregex(fitsSpectrum,'([^\/]+)\.fits$',/extract,/subexpr))[1] + '.inp'
		; should generally just be /some/full/path/_init.inp
	ENDIF
	IF strlowcase(strmid(outFile,3,4,/REVERSE_OFFSET)) ne '.inp' THEN outFile = outFile + ".inp"
	if ~ keyword_set(stop_on_err) THEN stop_on_err = 0
	; -------------------------------------------------------------
	; Ensure that we understand all of the abundances passed
	;
	free_abund = []
	FOR i=0, n_elements(freeElems)-1 DO free_abund = [free_abund,WHERE(STREGEX(elems,"^"+freeElems[i]+"[ ]*$",/fold_case,/boolean) EQ 1)]
	IF N_ELEMENTS(free_abund) NE N_ELEMENTS(freeElems) THEN RETURN, {status: -1, file: outfile, msg: "Unable to locate all specified elements!"}

	; -------------------------------------------------------------
	; Open the FITS file and read in the reduced data and the 
	;	header from the original raw file (in first extension).
	;
	CATCH, err ; Catch any error reading the FITS file
	
	IF err NE 0 THEN BEGIN
		IF verbose THEN BEGIN
			PRINT, 'Error index: ', err
			PRINT, 'Error message: ', !ERROR_STATE.MSG
		ENDIF
		
		; Can't open the file, so can't continue
		IF stop_on_err THEN stop
		; If we were able to get the table and tblHeader, we can do without the meterTable
		IF keyword_set(gotTable) THEN BEGIN
			IF verbose THEN BEGIN
				PRINT,"Meter table seems to be missing from " + file
			ENDIF
		ENDIF ELSE BEGIN
			RETURN, {status: err, file: outFile, msg: STRING('Unable to open spectrum at: ',fitsSpectrum)}
		ENDELSE
		CATCH, /CANCEL
	ENDIF
	
	; Get the primary FITS header, we don't need the image in the extension, just the header
	junk = mrdfits(fitsSpectrum, 0, header, /SILENT)
	
	; Read a specific type of extraction from the file by naming the extension or just get
	;	the first extension regardless of how it was extracted:
	; spectrum = mrdfits(fitsSpectrum,'optimal extraction',tblHeader, /SILENT)
	table = mrdfits(fitsSpectrum,1,tblHeader, /SILENT)
	
	gotTable = 1 ; in case this file doesn't have the meter spectra but everything else is ok
	
	; Read in the exposure meter extension...this will have the flux weighted barycentric correction
	;
	meterTable = mrdfits(fitsSpectrum,2,meterHeader, /SILENT)
	
	gotMeterSpectra = 1
	CATCH, /CANCEL ; end of catch block for reading FITS file
	
	; -------------------------------------------------------------
	; Check to see ensure that the required columns are here
	;
	required = ['WAVELENGTH','UNCERTAINTY','SPECTRUM','CONTINUUM','BARY_WAVELENGTH']
	cols = tag_names(table)
	FOR i=0, n_elements(required)-1 DO BEGIN
		idx = where(cols eq required[i], count)
		IF count eq 0 THEN BEGIN
			return, {status: -2, file: outFile, msg: STRING("Missing column ",required[i]," in FITS file.")}
		ENDIF
	ENDFOR
	
	; -------------------------------------------------------------
	; Make sure that the wavelengths are in Angstroms or Nanometers
	;	We only check the uncorrected wavelengths and assume that
	;	the barycentric corrected wavelengths are in the same units
	;	if not, then this will need to be updated to check both.
	;
	waveunits = 'unknown'
	num_cols = sxpar(tblHeader,'TFIELDS',count=nFound)
	IF nFound eq 0 THEN BEGIN
		return, {status: -1, file: outFile, msg: STRING("Column definitions not found in FITS Table.")}
	ENDIF
	FOR i=1, num_cols DO BEGIN
		colname = sxpar(tblHeader,STRING('TTYPE',i,format="(A,I0)"),count=nFound)
		IF nFound eq 0 THEN CONTINUE
		if colname ne 'wavelength' THEN CONTINUE
		waveunits = sxpar(tblHeader,STRING('TUNIT',i,format="(A,I0)"),count=nFound)
		if nFound eq 0 THEN waveunits = 'unknown'
		BREAK
	ENDFOR
	
	wavefactor = 1d0 ; for angstroms we are already in the right units
	IF waveunits eq 'nanometers' THEN BEGIN
		wavefactor = 10d0
	ENDIF ELSE IF waveunits ne 'angstroms' THEN BEGIN
		return, {status: -3, file: outFile, msg: STRING("Unknown units in wavelength table: ",waveunits)}
	ENDIF
	
	; -------------------------------------------------------------
	; Get the shape of the arrays and then reform them to the same
	;	shape so that we can access them normally.
	;
	shape = size(table.wavelength)
	nCols = shape[1]
	nOrds = shape[2]
	wavelengths = reform(table.wavelength * wavefactor, nCols, nOrds)
	bary_wavelengths = reform(table.bary_wavelength * wavefactor, nCols, nOrds)
	raw_spectrum = reform(table.spectrum, nCols, nOrds)
	uncertainty = reform(table.uncertainty, nCols, nOrds)
	continuum = reform(table.continuum, nCols, nOrds)
	; tellurics = reform(table.tellurics, nCols, nOrds)
	
	; The spectrum is not continuum normalized, so we need to fix that using the continuum
	;	column before continuuing
	;
	;	We wrap it in a catch block in case the continuum is 0 anywhere
	;
	CATCH, err ; Catch any error dividing spectrum by continuum
	
	IF err NE 0 THEN BEGIN
		IF verbose THEN BEGIN
			PRINT, 'Error index: ', err
			PRINT, 'Error message: ', !ERROR_STATE.MSG
		ENDIF
		
		; Can't open the file, so can't continue
		IF stop_on_err THEN STOP
		return, {status: err, file: outFile, msg: STRING('Unable to continuum normalize the spectrum for: ',fitsSpectrum,': ',!ERROR_STATE.MSG)}
		CATCH, /CANCEL
	ENDIF
	
	spectrum = raw_spectrum / continuum
	uncertainty = uncertainty / continuum
	
	CATCH, /CANCEL
	
	; -------------------------------------------------------------
	; The barycentric radial velocity is in the meter extension header.
	;	The barycentric corrected wavelengths are in the primary table.
	;	We can get both, preferring the flux weighted meter barycentric correction
	;	but accepting one calculated from the difference in the two wavelength scales.
	;
	;	'Observer towards star' == positive in blueshifted direction
	;
	
	IF gotMeterSpectra THEN BEGIN
		; The meter table currently uses Heirarchical FITS keywords, so we have to
		;	find the key in a funny way
		heirKeyData = stregex(meterHeader,'HIERARCH (.*) = (.*)',/extract,/subexpr)
		ourKey = (where(heirKeyData[1,*] eq 'wtd_mdpt_bc'))[0] ; this is 'z', we want bc in m/s
		bary_vel = double(heirKeyData[2,ourKey]) * c_light
		bc_vel_type = "Weighted Midpoint Bary Vel"
	ENDIF ELSE BEGIN
		bv_order = 40
		bv_pixel = 3300
		max_bary_vel = 40.0
		bary_vel = c_light * ((bary_wavelengths[bv_pixel,bv_order] / wavelengths[bv_pixel,bv_order])-1.0d0) ; km/s
	
		IF abs(bary_vel) GT max_bary_vel THEN BEGIN
			IF verbose THEN BEGIN
				PRINT,"Abnormally high barycentric velocity shift: ", bary_vel, " km/s", format="(A,F0.1,A)"
			ENDIF
			IF stop_on_err THEN STOP
			RETURN, {status: -4, file: outFile, msg: STRING('Bad barycentric velocity: ',bary_vel, " km/s", format="(A,F0.1,A)")}
		ENDIF
		bc_vel_type = "Diff btwn bary corrected & obs waves"
	ENDELSE
	
	; -----------------------------------------------------------------
	; Get the peculiar velocity of the star in the barycentric frame
	; 	We use cross-correlation with Mg I b segment to find radial velocity
	;	of star with respect to BARYCENTER.  This is in traditional redshfit
	;	notation, so the star toward observer is NEGATIVE.
	;
	;	We wrap it in a catch block in case something goes wrong in the cross-correlation
	;
	CATCH, err ; Catch any error calculating the total velocity offset
	
	IF err NE 0 THEN BEGIN
		IF verbose THEN BEGIN
			PRINT, 'Error index: ', err
			PRINT, 'Error message: ', !ERROR_STATE.MSG
		ENDIF
		
		; Can't open the file, so can't continue
		IF stop_on_err THEN stop
		return, {status: err, file: outFile, msg: STRING('Unable to cross correlate spectrum for: ',fitsSpectrum,': ',!ERROR_STATE.MSG)}
		CATCH, /CANCEL
	ENDIF
	
	obs_waves = wavelengths
	vactoair,obs_waves  ; EXPRES Wavelengths are stored in vacuum wavelengths
	
	total_rv = get_expres_rv(obs_waves,spectrum,stop_on_err=stop_on_err)
	IF total_rv eq -9999 THEN MESSAGE, STRING('Unable to determine Total RV!')
	peculiar_rv = bary_vel + total_rv
		
	IF verbose THEN BEGIN
		; Doing this as a sanity check
		tmp2 = bary_wavelengths
		vactoair,tmp2
		peculiar_from_bary = get_expres_rv(tmp2, spectrum, stop_on_err=stop_on_err)

		print,"Peculiar RV: ",peculiar_rv," Barycentric Correction: ",bary_vel," Total (earth-centric) RV: ",total_rv, format="(A,F0.4,A,F0.4,A,F0.4)"
		print,"Peculiar RV from barycentric wavelengths (Should match peculiar): ",peculiar_from_bary,format="(A,F0.4)"
	ENDIF
	CATCH, /CANCEL
	
	; ------------------ CHOOSE WAVELENGTH SCALE TO USE --------------------
	; We will be using the observed wavelength scale with the total RV correction applied
	;	Alternately, we could use the barycentric corrected wavelengths with just the
	;	peculiar RV correction applied.
	;
	waves = obs_waves * (1.0d0 - (total_rv/c_light))
	
	; ------------- COLLECT SEGMENTS, BUILD LINE/SPECIES LISTS -------------
	; 
	; This is the meat of the prep...we look at our pre-defined line list segments
	;	and attempt to find observed spectrum over the same wavelength ranges.
	;	We then build up our observed spectral chunks and apply a telluric mask.
	;	For EXPRES, it may be better to use the Telluric mask that is already
	;	a part of the reduction process for this observation sinc it is much more
	;	sophisticated, however, for now we will just use what we have been using.
	;
	nSegments = n_elements(segFiles)
	nSegFiles = nSegments
	
	; It is possible for a single file to contain more than one segment.
	for i=0, nSegFiles-1 do begin
		restore, segFiles[i] ; sme
		segsInFile = n_elements(sme.wran[0,*])
		nSegments = nSegments + segsInFile - 1 ; count segments beyond first segment in file
	endfor
	
	; Initialize an array for wavelength ranges and pixel offsets for each of the segments.
	wran = fltarr(2,nSegments)
	wind = lonarr(nSegments)
	segBase = 0
	
	; The first (and last) valid wavelength of each order (i.e. where there is sepctrum) is different,
	; and does not always correspnod to the wavelengths in the zeroth and final columns.  Therefore,
	; we need to get the *real* start/end points for each order.
	;
	; Additionally, the EXPRES spectra have a lot of overlap for most orders, but the edges of  the
	;	orders are often low S/N.  We will trim a bit of each order (below relative order 80)
	;	to ensure we are just getting the good stuff.
	;
	oTrimLimit = 81 ; Trim orders 80 and below
	oTrimLeft = 3 ; Angstroms
	oTrimRight = 3 ; Angstroms
	
	pixRngs = intarr(2,nOrds)
	orderWaveMids = fltarr(nOrds)
	for oIdx=0, nOrds-1 do begin
		oMask = where(FINITE(spectrum[*,oIdx]))
		IF (oIdx lt oTrimLimit) THEN BEGIN
			startWave = waves[oMask[0],oIdx] + oTrimLeft
			endWave = waves[oMask[n_elements(oMask)-1],oIdx] - oTrimRight
			
			; Adjust the mask to be oTrimLeft and oTrimRight inside the valid spectral range for the order
			oMask = where(waves[*,oIdx] ge startWave AND waves[*,oIdx] le endWave)
		ENDIF
		
		pixRngs[*,oIdx] = [oMask[0],oMask[n_elements(oMask)-1]]
		orderWaveMids[oIdx] = 0.5 * (waves[pixRngs[0,oIdx],oIdx] + waves[pixRngs[1,oIdx],oIdx]) ; wavelength midpoint of each order
	endfor
	
	fullWaveRange = minmax(waves[where(FINITE(spectrum))])
	
	; Find the wavelength borders of the spectral orders
	;
	orderMins = min(waves[pixRngs,*],dimension=1,max=orderMaxs)	; min & max wavelength for each order
	wranges = reform([orderMins,orderMaxs],[nOrds,2]) 	; wranges[oIdx,*] gives min,max for each order oIdx
	
	curSeg = -1
	
	CATCH, err ; Catch any error in prepping the segments
	
	IF err NE 0 THEN BEGIN
		IF verbose THEN BEGIN
			PRINT, 'Error index: ', err
			PRINT, 'Error message: ', !ERROR_STATE.MSG
		ENDIF
		
		; Can't extract segments properly so exit
		IF stop_on_err THEN stop
		return, {status: err, file: outFile, msg: STRING("Error prepping segment (",curSeg,"): ",!ERROR_STATE.MSG)}
		CATCH, /CANCEL
	ENDIF
	
	; margin padding
	; Look for orders where our wavelength range has at least PAD Å of spectrum that overlaps
	; the wavelength range we are looking for.  If not, then we sometimes get just a couple of
	; pixels that are low S/N anyway.
	;
	margin_pad = 1.0 ; Å
	
	;Loop through segment FILES, extracting observed spectrum and atomic data.
	for i=0, nsegFiles-1 do begin
		curSeg = i
		restore, segFiles[i]
		segsInFile = n_elements(sme.wran[0,*])
		
		; Loop through segments _IN_ File and add them all to a single segment
		; !!! Have to set choose the right orders here !!!
		incompleteChunk = 0
		extraSegs = 0
		FOR j=0, segsInFile-1 DO BEGIN
			;Set wavelength range to model in segment.
			wlo = max([floor(sme.wran[0,j]),fullWaveRange[0]])
			whi = min([ceil(sme.wran[1,j]),fullWaveRange[1]])
			
			; Find all orders that COMPLETELY CONTAIN this segment
			;	If we don't find any, then we will instead get all overlapping ones
			orderIdxs = WHERE( (wranges[*,0] le wlo-margin_pad and $ ; blue edge of order less than blue edge of segment
								wranges[*,1] ge wlo+margin_pad) and $ ; red edge of order greater than blue edge of segment
								(wranges[*,0] le whi-margin_pad and $ ; blue edge of order less than red edge of segment
								wranges[*,1] ge whi+margin_pad),nOrdersFound ) ; blue edge of order greater than red edge of segment
			
			IF (nOrdersFound eq 0) THEN BEGIN
				; Find all of the orders that OVERLAP this segment
				orderIdxs = WHERE( (wranges[*,0] le wlo-margin_pad and wranges[*,1] ge wlo+margin_pad) or $
									(wranges[*,0] le whi-margin_pad and wranges[*,1] ge whi+margin_pad),nOrdersFound )
			ENDIF
			
			IF nOrdersFound eq 0 THEN BEGIN					;check for no spectrum
				STOP, 'no spectrum for wavelength range in line data file: ' + segFiles[i]
				PRINT,i
				CONTINUE
			ENDIF
			
			orderIdx = orderIdxs[0]
			pixMin = pixRngs[0,orderIdx]
			pixMax = pixRngs[1,orderIdx]
			quarterPixel = 0.25 * (waves[pixMin+1,orderIdx] - waves[pixMin,orderIdx])
			; as we add more spectral coverage from orders, we update this (second term is 1/4 pixel in angstroms)
			wlo_current = wlo - quarterPixel
			
			oWave = []
			oSOB = []
			oUOB = []
			oMOB = []
			oWinds = 0
			
			;stop,"Found ",nOrdersFound," orders for seg ",i,".",j
			
			for subIdx=0, nOrdersFound-1 do begin
				if wlo_current ge whi-margin_pad then break ; we've filled the whole segment
				orderIdx = orderIdxs[subIdx]
				pixMin = pixRngs[0,orderIdx]
				pixMax = pixRngs[1,orderIdx]
				quarterPixel = 0.25 * (waves[pixMin+1,orderIdx] - waves[pixMin,orderIdx])
				
				
				; don't run off the bottom of the order
				; we shift by a quarter of a pixel down and use GT rather than using GE to avoid
				;  overlapping pixels when we stitch two orders together.
				;
				wlo_current = wlo_current > waves[pixMin,orderIdx] - quarterPixel
				
				; Extract observed spectrum in the specified wavelength range
				waveIdxs = where(waves[*,orderIdx] gt wlo_current and waves[*,orderIdx] le whi and waves[*,orderIdx] le waves[pixMax,orderIdx], nwhr)
				if ( nwhr eq 0 ) then begin
				  stop, 'Logic error - no spectrum for valid wavelength range'
				endif
				
				wave = reform(waves[waveIdxs,orderIdx])
				if n_elements(wave) lt 10 then continue; don't bother grabbing less than 10 pixels
				
				whi_prev = max([oWave,0])
				whi_current = max(wave)
				
				; Add up all of the pixels from multiple-order segments
				oWinds += nwhr
				
				oWave = [oWave,wave]
				oSOB = [oSOB,reform(spectrum[waveIdxs,orderIdx])]
				oUOB = [oUOB,reform(uncertainty[waveIdxs,orderIdx])]
				
				; wavelengths of masking segment in this range
				sunWhr = where(sme.wave ge wlo_current and sme.wave le whi_current, nsunWhr)
				;stop, i, j, subIdx, nsunWhr, n_elements(wave)
				
				; not really 'solar' mask.  Instead applies 'high' resolution mask
				;	to 'low' resolution observation
				;
				apply_solar_mask, sme.wave[sunWhr], sme.mob[sunWhr], wave, tmpMOB
				
				; Get the telluric mask, corrected for radial velocity
				telluric = get_telluric_mask(wave,bary_vel,total_rv)
				
				; modify the mask by masking out tellurics
				tmpMOB = tmpMOB * telluric
				oMOB = [oMOB,tmpMOB]
				
				wlo_current = whi_current
			endfor
			
			if (i eq 0 and j eq 0) then begin
				wind[segBase + j + extraSegs] = oWinds - 1
				smeWave = oWave
				sob  = oSOB
				uob  = oUOB
				mob  = oMOB
			endif else begin
				wind[segBase + j + extraSegs] = wind[segBase + j -1] + oWinds
				smeWave = [smeWave, oWave]
				sob  = [sob, oSOB]
				uob  = [uob, oUOB]
				mob  = [mob, oMOB]
			endelse
			wran[*,segBase+j] = [min(oWave), max(oWave)]			;save range for output
		endfor
		
		;Extract atomic data from the file.
		if ( i eq 0 ) then begin			;first pass, create arrays
		  species = sme.species
		  atomic  = sme.atomic
		  lande   = sme.lande
		  depth   = sme.depth
		  lineref = sme.lineref
		endif else begin					;later passes (IN 1 FILE), append to arrays
		  species = [species, sme.species]
		  atomic  = transpose([transpose(atomic),  transpose(sme.atomic)])
		  lande   = [lande,   sme.lande  ]
		  depth   = [depth,   sme.depth  ]
		  lineref = [lineref, sme.lineref]
		endelse
		
		segBase = segBase + segsInFile
	endfor
	
	CATCH, /CANCEL
	
	nlin = n_elements(species)				;total number of lines
	IF verbose THEN BEGIN
		print,nlin," Lines in the spectral model",format="(I0,A)"
	ENDIF
	
	; ---------------- Create Input File ------------------
	
	; Loop through SME masked/tuned segments
	
	; Grab bits of spectrum in same wavelength range as chunk
		; Match to line/atomic data in tuned segments
	
	; Apply solar mask ?
	
	; Get other neccessary SME lines/depths/references/species/etc
	
	;
	; Set initial parameters for SME
	;
	version = 3.0
	id = systime()
	obs_name = strjoin(strsplit(sxpar(header,'OBJECT'),' ',/EXTRACT),'') ; name without spaces
	obs_type = -1 ; ??
	obsnm = sxpar(header, 'OBS_ID')
	obsdate = sxpar(header, 'DATE-OBS')
	iptype = 'gauss'	; instrumental profile
	ipres = sxpar(header,'RESOLUTN',count=resFound)		; resolution
	
	IF resFound NE 1 THEN BEGIN
		if verbose THEN PRINT,"WARNING: No IP resolution found, using default"
		ipres = expres_ipres
	ENDIF
	
	; Global Free Parameters
	; first solve with VMAC free, VSINI=0, then solve for VSINI using VMAC relation for VMAC
	glob_free = ['TEFF','GRAV','FEH', 'VMAC'] 
	; Initial values (solar...could guess better since we know these are cool)
	teff = 5777
	grav = 4.44
	feh = 0.0
	vsini = 0.0 ; let all broadening be absorbed by vmac
	jbsme_guess_b, obs_name, teff=teff, logg=logg, vmac=vmac
	
	vmic = 0.85		; fixed
	vmac = vmac		; 3.558==solar, all rotational broadening will initially be absorbed here
	
	abund = sme.abund ; from one of the line list files ?
	
	; Individual Abundances: array elements are Atomic Number - 1
	ab_free = intarr(99) * 0
	ab_free[free_abund] = 1
	
	; Radial Velocity - We have already 'zeroed out' the RV
	vrad = replicate(0.0, nSegments)
	vrad_flag = 0 ; solve for vrad in each segment
	
	; continuum scaling 
	;	-3=fixed global/normalized flux (scalar)
	;	-2=fixed global/absolute flux   (scalar)
	;	-1=free global					(scalar)
	;	 0=fit scalar to each segment	(array[sme.nseg])
	;	 1=fit line to each segment		(array[2,sme.nseg])
	;
	cscale_flag = 1
	cscale = replicate( 1.0, 2, nSegments )
	
	; Gamma 6 factor
	gam6 = 2.5
	
	; ?? Other SME settings ??
	nseg = nSegments
	accwi = 0.003	;
	accrt = 0.001	;
	clim = 0.01		;
	maxiter = 100	;
	chirat = 0.002	;
	
	nmu = 7			; number of mu angles for macroturbulence
	mu = rotate(sqrt(0.5*(2*dindgen(nmu)+1)/nmu), 5) ; actual mu angles
	
	atmo_pro = 'interp_atmo_grid' ; 'interpkrz2'
	atm_type = 'RHOX'

	; Save new SME structure:
	sme = { $
		version:     version,     $
		id:          id,          $
		teff:        teff,        $
		grav:        grav,        $
		feh:         feh,         $
		vmic:        vmic,        $
		vmac:        vmac,        $
		vsini:       vsini,       $
		vrad:        vrad,        $
		vrad_flag:   vrad_flag,   $ 
		cscale:      cscale,      $
		cscale_flag: cscale_flag, $
		gam6:        gam6,        $
		accwi:       accwi,       $
		accrt:       accrt,       $
		clim:        clim,        $
		maxiter:     maxiter,     $
		chirat:      chirat,      $
		nmu:         nmu,         $
		nseg:        nseg,        $
		abund:       abund,       $
		species:     species,     $
		atomic:      atomic,      $
		lande:       lande,       $
		depth:       depth,       $
		lineref:     lineref,     $
		wran:        wran,        $
		mu:          mu,          $
		glob_free:   glob_free,   $
		ab_free:     ab_free,     $
		atmo_pro:    atmo_pro,    $
		atmogrid_file: atmogrid_file, $
		atm_type:	 atm_type,	  $
		wave:        smeWave,     $
		wind:        wind,        $
		sob:         sob,         $
		uob:         uob,         $
		mob:         mob,         $
		obs_name:    obs_name,    $
		obs_type:    obs_type,    $
		iptype:      iptype,      $
		ipres:       ipres,       $
		obsnm:       obsnm,       $
		obsdate:	 obsdate,	  $
		baryvel:	 bary_vel,	  $
		totalrv:	 total_rv,	  $
		starrv:      peculiar_rv, $
		bcvel_type:  bc_vel_type  $
	}
	
	save, sme, file=outFile
	print, "Saved " + outFile
	return, {status: 0, file: outFile, msg: "OK"}	
END