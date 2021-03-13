#!/usr/bin/env python3

import pathlib

####### Import ##########

uvfits = [str(s) for s in pathlib.Path('uvfits').glob('*.uvfits')]
ms = [s.replace('uvfits', 'ms') for s in uvfits]
for uv, vis in zip(uvfits, ms): 
     importuvfits(fitsfile = uv, vis = vis)
vis = 'ATA-20210221.ms'
concat(vis = ms, concatvis = vis)

####### Flagging ##########

# manual flagging of some data
# plotms(vis = vis, xaxis = 'time', yaxis = 'amp', coloraxis = 'field', avgchannel = '64', avgtime = '10', correlation = 'XX,YY')
# manual flags
flagmanager(vis = vis, mode = 'restore', versionname = 'manual')

####### Calibration ##########

setjy(vis = vis, field = '3C286', standard = 'Perley-Butler 2017',
          model = '3C286_C.im', usescratch = False, scalebychan = True)
setjy(vis = vis, field = '3C147', standard = 'Perley-Butler 2017',
          model = '3C147_C.im', usescratch = False, scalebychan = True)

# Bandpass calibration on 3C84

# preliminary phase calibration for bandpass calibration
gaincal(vis = vis, caltable = f'3C84.G0',
            spw = '0:4~60',
            field = f'3C84', refant = '"1h"',  
            gaintype = 'G', calmode = 'p',  
            solint = 'inf', combine = 'scan,obs',  
            minblperant = 1)

# we skip delay calibration because only baselines
# where refant is present are used for delay calibration

# bandpass calibration
bandpass(vis = vis, caltable = 'cal.B0',
         field = '3C84', refant = '"1h"',
         bandtype = 'B', solint = 'inf,4ch',
         combine = 'scan,obs', 
         gaintable = ['3C84.G0'])

# X-Y phase offset calibration on 3C84
gaincal(vis = vis, caltable = f'3C84.G1',
            spw = '0:4~60',
            field = f'3C84', refant = '"1h"',  
            gaintype = 'G', calmode = 'p',  
            solint = 'inf', combine = 'scan,obs',  
            minblperant = 1,
            gaintable = ['cal.B0'])


# timeranges to combine scans
timeranges = {'3C84': ['02:00:00~04:00:00', '04:00:00~05:50:00', '06:00:00~08:20:00'],
              '3C286' : ['06:00:00~08:00:00', '08:00:00~09:00:00', '09:00:00~11:00:00'],
              '3C147' : ['07:00:00~08:10:00', '08:30:00~09:30:00', '09:30:00~10:30:00'],
              '3C345' : ['10:00:00~12:00:00'],
              '2202+422' : ['14:00:00~17:00:00', '17:00:00~20:00:00']}

# Common XY phase calibrations
rmtables('cal.T1')
for cal in timeranges:
    for timerange in timeranges[cal]:
        gaincal(vis = vis, caltable = f'cal.T1', field = cal, 
                spw = '0:4~60', solint = 'inf', combine = 'scan,obs', 
                refant = '"1h"', gaintype = 'T', calmode = 'p', 
                gaintable = ['cal.B0', '3C84.G1'], minblperant = 1,
                timerange = timerange,
                append = True)

# Amplitude calibration
rmtables('cal.G2')
for cal in timeranges:
    for timerange in timeranges[cal]:
        extra = {'solnorm' : False} if cal == '3C286' \
          else {'preavg' : 10} if cal in ['3C84'] \
          else {'preavg' : 10, 'minsnr' : 1} if cal in ['3C345', '2202+422'] \
          else {}
        gaincal(vis = vis, caltable = f'cal.G2', field = cal, 
                spw = '0:4~60', solint = 'inf', combine = 'scan,obs', 
                refant = '"1h"', gaintype = 'G', calmode = 'a', 
                gaintable = ['cal.B0', '3C84.G1', 'cal.T1'],
                gainfield = ['', '', cal],
                minblperant = 1,
                timerange = timerange,
                append = True, **extra)

# Flux calibration
rmtables('cal.flux2')
fluxscale(vis = vis,
              caltable = 'cal.G2', fluxtable = 'cal.flux2',
              reference = '3C286',
              transfer = ['3C84', '3C147', '3C345', '2202+422'],
              incremental = False)

# Apply
for cal in timeranges:
    applycal(vis = vis, field = cal, gaintable = ['cal.B0', '3C84.G1', 'cal.T1', 'cal.flux2'],
             gainfield = ['', '', cal, cal], calwt = False)

applycal(vis = vis, field = 'CasA', gaintable = ['cal.B0', '3C84.G1', 'cal.T1', 'cal.flux2'],
             calwt = False)
applycal(vis = vis, field = 'CygA', gaintable = ['cal.B0', '3C84.G1', 'cal.T1', 'cal.flux2'],
             calwt = False)

statwt(vis = vis, datacolumn = 'data')

####### Data reduction ##########

# Extract and reduce CasA and CygA fields
for field in ['CasA', 'CygA']:
    mstransform(vis = vis, outputvis = f'ATA-20210221-{field}.ms',
                field = field, datacolumn = 'corrected',
                chanaverage = True, chanbin = 4,
                timeaverage = True, timebin = '10s')

####### Imaging ##########

tclean(vis = 'ATA-20210221-CasA.ms', imagename = 'CasA',
           mask = 'CasA-mask.mask',
           specmode = 'mfs', 
           niter = 300, gain = 0.1, threshold = '2Jy', 
           imsize = 256, 
           cell = '6arcsec', stokes = 'I',
           deconvolver = 'multiscale', scales=[0, 5, 15, 45], smallscalebias=0.9, 
           weighting = 'briggs', robust = 0.5, 
           pbcor = False, 
           savemodel = 'modelcolumn')

tclean(vis = 'ATA-20210221-CygA.ms', imagename = 'CygA',                                                                                                                                 
           mask = 'CygA-mask.mask', 
           niter = 1000, gain = 0.1, threshold = '1.5Jy',  
           imsize = 128,  
           cell = '6arcsec', stokes = 'I', 
           deconvolver = 'multiscale', scales=[0, 2, 5, 10], smallscalebias=0.9,  
           weighting = 'briggs', robust = 0.5,  
           pbcor = False,  
           savemodel = 'modelcolumn')

####### Plots ##########

plotms(vis = vis, coloraxis = 'field',
           avgchannel = '64', avgtime = '25',
           xaxis = 'time', yaxis = 'amp',
           correlation = 'XX,YY',
           plotfile = '/tmp/amp_uncal.png')

plotms(vis = vis, field = 'CasA,CygA', coloraxis = 'field',
           avgchannel = '64', avgtime = '1000',
           xaxis = 'uwave', yaxis = 'vwave',
           plotfile = '/tmp/uv_coverage.png')

plotms(vis = 'cal.B0', xaxis = 'frequency', yaxis = 'gainphase',
           coloraxis = 'antenna1', xconnector = 'line',
           plotfile = '/tmp/bandpass_phase.png')

plotms(vis = 'cal.B0', xaxis = 'frequency', yaxis = 'gainamp',
           coloraxis = 'antenna1', xconnector = 'line',
           plotfile = '/tmp/bandpass_amp.png')

plotms(vis = 'cal.T1', xaxis = 'time', yaxis = 'gainphase',
           xconnector = 'line', coloraxis = 'field',
           plotfile = '/tmp/phase_cal.png')

plotms(vis = 'cal.flux2', xaxis = 'time', yaxis = 'gainamp',
           xconnector = 'line', coloraxis = 'field',
           plotfile = '/tmp/gain_cal.png')

plotms(vis = vis, xaxis = 'real', yaxis = 'imag', avgchannel = '64',
           avgtime = '25', correlation = 'XX,YY', coloraxis = 'field',
           field = '3C84,3C286,3C147,3C345,2202+422',
           xdatacolumn = 'corrected', ydatacolumn = 'corrected',
           plotfile = '/tmp/calibrators_IQ.png')

plotms(vis = vis, xaxis = 'uvwave', yaxis = 'amp', avgchannel = '64',  
           avgtime = '25', correlation = 'XX,YY', coloraxis = 'field',  
           field = 'CasA',  
           ydatacolumn = 'corrected', 
           plotfile = '/tmp/CasA-uvwave.png')

plotms(vis = vis, xaxis = 'uvwave', yaxis = 'amp', avgchannel = '64',  
           avgtime = '25', correlation = 'XX,YY', coloraxis = 'field',  
           field = 'CygA',  
           ydatacolumn = 'corrected', 
           plotfile = '/tmp/CygA-uvwave.png')
