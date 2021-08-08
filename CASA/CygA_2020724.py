#!/usr/bin/env python3
#
# This script was done for CASA 6.2.0
#

import pathlib

dirs = ['20210723_3C286', '20210724_cygA', '20210725_cygA']
vis = ['20210723_3C286/20210723.ms', '20210724_cygA/20210724.ms',
       '20210725_cygA/20210725.ms']

def snr(img, box='100,100,230,230'):
    return imstat(img)['max'][0]/imstat(img, box=box)['rms'][0]

####### Import ##########

for v,d in zip(vis, dirs):
    uvfits = [str(s) for s in pathlib.Path(d).glob('*.uvfits')]
    ms = [s.replace('uvfits', 'ms') for s in uvfits]
    for uv, w in zip(uvfits, ms):
        importuvfits(fitsfile=uv, vis=w)
    concat(vis=ms, concatvis=v)
    for vv in ms:
        rmtables(vv)
    tb.open(v, nomodify=False)
    a = tb.getcol('OBSERVATION_ID')
    a[:] = 0
    tb.putcol('OBSERVATION_ID', a)
    tb.close()

####### Flagging ##########

for v in vis:
    flagdata(vis=v, mode='tfcrop', datacolumn='data', action='apply',
             flagbackup=False, freqcutoff=3.5, timecutoff=3.0)

# Zenith keyhole
flagdata(vis=vis[1], mode='manual', scan='697~707', flagbackup=False)
flagdata(vis=vis[2], mode='manual', scan='1033~1041', flagbackup=False)

# Not visible
flagdata(vis=vis[1], mode='manual', scan='689~696', flagbackup=False)
flagdata(vis=vis[2], mode='manual', scan='1017~1023', flagbackup=False)

# Other problems
flagdata(vis=vis[1], mode='manual', scan='201~216', flagbackup=False)
flagdata(vis=vis[1], mode='manual', scan='225~232', flagbackup=False)
flagdata(vis=vis[2], mode='manual', scan='283~288', flagbackup=False)
flagdata(vis=vis[2], mode='manual', scan='545~548', flagbackup=False)
flagdata(vis=vis[2], mode='manual', scan='551~552', flagbackup=False)
flagdata(vis=vis[2], mode='manual', scan='568', flagbackup=False)
flagdata(vis=vis[2], mode='manual', scan='537~544', flagbackup=False)

####### Calibration ##########

for v in vis:
    setjy(vis=v, field='3C286', standard='Perley-Butler 2017',
          model='3C286_X.im', usescratch=False, scalebychan=True)

refant = '"2j"'

# Bandpass calibration on 3C84

timerange = [None,
             '2021/07/25/15:10:00~17:00:00',
             '2021/07/26/15:15:00~17:00:00']
for j in range(1, len(vis)):
# preliminary phase calibration for bandpass calibration
    gaincal(vis=vis[j], caltable=f'{dirs[j]}/3C84.G0', spw='0:4~60',
            field='3C84', refant=refant,
            timerange=timerange[j],
            gaintype='G', calmode='ap',
            solint='inf', combine='scan',
            minblperant=1)
    bandpass(vis=vis[j], caltable=f'{dirs[j]}/3C84.B0',
             field='3C84', refant=refant,
             timerange=timerange[j],
             bandtype='B', solint='inf',
             combine='scan',
             gaintable=[f'{dirs[j]}/3C84.G0'],
             minblperant=1, minsnr=2)
    # X-Y phase offset calibration on 3C84
    gaincal(vis=vis[j], caltable=f'{dirs[j]}/3C84.G1',
            field='3C84', refant=refant,
            timerange=timerange[j],
            gaintype='G', calmode='p',
            solint='inf', combine = 'scan',
            minblperant=1,
            gaintable=[f'{dirs[j]}/3C84.B0'])


# the X-Y and bandpass from dirs[1] look better (there was
# much greater SNR), so we use those

bp_cal = f'{dirs[1]}/3C84.B0'
xy_cal = f'{dirs[1]}/3C84.G1'

# timeranges to combine scans
timeranges = [
    {},
    {'J2007+404':
         ['07:20:00~08:00:00', '08:10:00~09:00:00',
          '09:00:00~09:45:00', '09:45:00~10:40:00',
          '10:40:00~11:20:00', '11:20:00~12:10:00',
          '12:10:00~13:00:00', '13:00:00~13:50:00',
          '13:50:00~15:00:00'],
      '3C84':
          ['8:50:00~9:10:00', '9:40:00~10:00:00',
           '10:30:00~10:50:00', '11:15:00~11:40:00',
           '12:05:00~12:20:00', '12:55:00~13:10:00',
           '13:45:00~14:00:00', '14:35:00~14:50:00',
           '15:00:00~15:42:00', '15:42:00~15:58:00',
           '15:58:00~17:00:00'],
    },
    {'J2007+404':
         ['1:00:00~2:05:00', '2:05:00~2:55:00',
          '2:55:00~4:00:00', '4:00:00~4:55:00',
          '4:55:00~5:45:00', '5:45:00~6:40:00',
          '6:40:00~7:40:00', '7:40:00~9:15:00',
          '9:15:00~10:00:00', '10:00:00~10:50:00',
          '10:50:00~11:40:00', '11:40:00~12:30:00',
          '12:30:00~13:15:00', '13:15:00~14:00:00',
          '14:00:00~15:00:00'],
      '3C84':
         ['8:20:00~8:40:00', '9:10:00~9:20:00',
          '10:00:00~10:10:00', '10:45:00~11:00:00',
          '11:30:00~11:50:00', '12:20:00~12:40:00',
          '13:10:00~13:30:00', '14:00:00~14:20:00',
          '14:50:00~15:00:00', '15:10:00~15:37:00',
          '15:37:00~15:53:00', '15:53:00~17:00:00'],
      '3C286':
         ['24:20:00~25:30:00', '25:50:00~26:20:00',
          '26:40:00~27:20:00', '27:40:00~28:20:00',
          '28:40:00~29:20:00', '29:40:00~30:00:00',
          '30:30:00~30:50:00', '31:20:00~31:50:00'],
    },
    ]

# Common XY phase calibrations
for j in range(1, len(dirs)):
    rmtables(f'{dirs[j]}/J2007+404.T0')
    for timerange in timeranges[j]['J2007+404']:
        gaincal(vis=vis[j], caltable=f'{dirs[j]}/J2007+404.T0',
                field='J2007+404',
                solint='inf', combine='scan', refant=refant,
                gaintype='T', calmode='p',
                gaintable=[bp_cal, xy_cal], minblperant=1,
                timerange=timerange, append=True)

    rmtables(f'{dirs[j]}/3C84.T0')
    for timerange in timeranges[j]['3C84']:
        gaincal(vis=vis[j], caltable=f'{dirs[j]}/3C84.T0', field='3C84',
                solint='inf', combine='scan', refant=refant,
                gaintype='T', calmode='p',
                gaintable=[bp_cal, xy_cal], minblperant=1,
                timerange=timerange, append=True)

for j in range(len(dirs)):
    rmtables(f'{dirs[j]}/3C286.T0')
    gaincal(vis=vis[j], caltable=f'{dirs[j]}/3C286.T0', field='3C286',
            solint='600' if j == 0 else '1200',
            combine='scan', refant=refant,
            gaintype='T', calmode='p',
            gaintable=[bp_cal, xy_cal], minblperant=1,
            timerange='18:00:00~33:00:00' if j != 2 else '18:00:00~24:20:00')

j = 2
for timerange in timeranges[j]['3C286']:
    gaincal(vis=vis[j], caltable=f'{dirs[j]}/3C286.T0', field='3C286',
            solint='inf', combine='scan', refant=refant,
            gaintype='T', calmode='p',
            gaintable=[bp_cal, xy_cal], minblperant=1,
            timerange=timerange, append=True)

# Amplitude calibration
for j in range(1, len(dirs)):
    rmtables(f'{dirs[j]}/amplitude.G0')
    for timerange in timeranges[j]['J2007+404']:
        gaincal(vis=vis[j], caltable=f'{dirs[j]}/amplitude.G0',
                field='J2007+404',
                solint='inf', combine='scan', refant=refant,
                gaintype='G', calmode='a',
                gaintable=[bp_cal, xy_cal,
                           f'{dirs[j]}/J2007+404.T0'],
                minblperant=1,
                timerange=timerange, append=True)
    for timerange in timeranges[j]['3C84']:
        # We force a flux of 37 Jy on 3C84
        # This is done because there is a lot of time elapsed between
        # the scans of 3C286 and the other calibrators in this observation,
        # so it is not possible to transfer the fluxes accurately.
        # The flux of 37 Jy has been transferred from vis[2]
        additional = {} if j != 1 else {'smodel': [37, 0, 0, 0]}
        gaincal(vis=vis[j], caltable=f'{dirs[j]}/amplitude.G0', field='3C84',
                solint='inf', combine='scan', refant=refant,
                gaintype='G', calmode='a',
                gaintable=[bp_cal, xy_cal,
                           f'{dirs[j]}/3C84.T0'],
                minblperant=1,
                timerange=timerange, append=True,
                **additional)

for j in range(len(dirs)):
    gaincal(vis=vis[j], caltable=f'{dirs[j]}/amplitude.G0', field='3C286',
            solint='600' if j == 0 else '1200',
            combine='scan', refant=refant,
            gaintype='G', calmode='a',
            gaintable=[bp_cal, xy_cal, f'{dirs[j]}/3C286.T0'], minblperant=1,
            timerange='18:00:00~33:00:00' if j != 2 else '18:00:00~24:20:00',
            append=True)

j = 2
for timerange in timeranges[j]['3C286']:
    gaincal(vis=vis[j], caltable=f'{dirs[j]}/amplitude.G0', field='3C286',
            solint='inf', combine='scan', refant=refant,
            gaintype='G', calmode='a',
            gaintable=[bp_cal, xy_cal,
                       f'{dirs[j]}/3C286.T0'],
            minblperant=1,
            timerange=timerange, append=True)


# Flux calibration

j = 1
rmtables(f'{dirs[j]}/flux.G0')
fluxscale(vis=vis[j], caltable=f'{dirs[j]}/amplitude.G0',
          fluxtable=f'{dirs[j]}/flux.G0', reference='3C84',
          transfer=['J2007+404'],
          incremental=False)
j = 2
rmtables(f'{dirs[j]}/flux0.G0')
fluxscale(vis=vis[j], caltable=f'{dirs[j]}/amplitude.G0',
          fluxtable=f'{dirs[j]}/flux0.G0', reference='3C286',
          transfer=['3C84', 'J2007+404'],
          timerange='2021/07/25/00:00:00~2021/07/26/03:36:00',
          incremental=False)
rmtables(f'{dirs[j]}/flux1.G0')
fluxscale(vis=vis[j], caltable=f'{dirs[j]}/amplitude.G0',
          fluxtable=f'{dirs[j]}/flux1.G0', reference='3C286',
          transfer=['3C84', 'J2007+404'],
          timerange='2021/07/26/03:36:00~2021/07/27/00:00:00',
          incremental=False)

# Apply
for j in range(len(vis)):
    clearcal(vis=vis[j])

# Cygnus A
j = 1
applycal(vis=vis[j], field='CygA',
         gaintable=[f'{dirs[j]}/flux.G0', f'{dirs[j]}/J2007+404.T0',
                    bp_cal, xy_cal],
         gainfield=['J2007+404', 'J2007+404', '', ''],
         calwt=False)
j = 2
applycal(vis=vis[j], field='CygA',
         gaintable=[f'{dirs[j]}/flux0.G0', f'{dirs[j]}/J2007+404.T0',
                    bp_cal, xy_cal],
         gainfield=['J2007+404', 'J2007+404', '', ''],
         timerange='2021/07/25/00:00:00~2021/07/26/03:36:00',
         calwt=False)
applycal(vis=vis[j], field='CygA',
         gaintable=[f'{dirs[j]}/flux1.G0', f'{dirs[j]}/J2007+404.T0',
                    bp_cal, xy_cal],
         gainfield=['J2007+404', 'J2007+404', '', ''],
         timerange='2021/07/26/03:36:00~2021/07/27/00:00:00',
         calwt=False)

# calibrators
for j in range(2):
    for calfield in ['3C286', 'J2007+404', '3C84']:
        if j == 0 and calfield != '3C286':
            continue
        applycal(vis=vis[j], field=calfield,
                 gaintable=[(f'{dirs[j]}/flux.G0' if j != 0
                             else f'{dirs[j]}/amplitude.G0'),
                             f'{dirs[j]}/{calfield}.T0',
                             bp_cal, xy_cal],
                 gainfield=[calfield, calfield, '', ''],
                 calwt=False)

j = 2
for calfield in ['3C286', 'J2007+404', '3C84']:
    applycal(vis=vis[j], field=calfield,
             gaintable=[(f'{dirs[j]}/flux0.G0' if j != 0
                        else f'{dirs[j]}/amplitude.G0'),
                        f'{dirs[j]}/{calfield}.T0',
                            bp_cal, xy_cal],
             gainfield=[calfield, calfield, '', ''],
             calwt=False,
             timerange='2021/07/25/00:00:00~2021/07/26/03:36:00')
    applycal(vis=vis[j], field=calfield,
             gaintable=[(f'{dirs[j]}/flux1.G0' if j != 0
                        else f'{dirs[j]}/amplitude.G0'),
                        f'{dirs[j]}/{calfield}.T0',
                            bp_cal, xy_cal],
             gainfield=[calfield, calfield, '', ''],
             calwt=False,
             timerange='2021/07/26/03:36:00~2021/07/27/00:00:00')

# Extract and concatenate calibrated visibilities
for j in range(len(vis)):
    for field in ['CygA', '3C286', 'J2007+404', '3C84']:
        if j == 0 and field != '3C286':
            continue
        rmtables(f'{dirs[j]}/{field}.ms')
        rmtables(f'{dirs[j]}/{field}.ms.flagversions')
        mstransform(vis=vis[j], outputvis=f'{dirs[j]}/{field}.ms',
                    field=field, datacolumn='corrected',
                    antenna='*&', correlation='XX,YY')

for field in ['CygA', '3C286', 'J2007+404', '3C84']:
    rmtables(f'{field}.ms')
    use_dirs = dirs if field == '3C286' else dirs[1:]
    concat(vis=[f'{d}/{field}.ms' for d in use_dirs], concatvis=f'{field}.ms')
    statwt(f'{field}.ms', datacolumn='data')

# Clean

# 3C84

tclean(vis='3C84.ms', imagename='3C84',
       mask='masks/3C84.mask',
       niter=100, gain=0.1, threshold='0.1Jy',
       imsize=512,
       cell='4arcsec', stokes='I',
       deconvolver='multiscale', scales=[0, 2, 5, 10],
       smallscalebias=0.9,
       weighting='briggs', robust=0.5,
       pbcor=True, savemodel='modelcolumn')

gaincal('3C84.ms', caltable='3C84.G0',
        solint='300', refant='"2j"', minblperant=1,
        gaintype='G', calmode='ap', combine='scan', minsnr=1)

applycal('3C84.ms', gaintable='3C84.G0', calwt=False)

tclean(vis='3C84.ms', imagename='3C84_selfcal1',
       mask='masks/3C84.mask',
       niter=1000, gain=0.1, threshold='0.1Jy',
       imsize=512,
       cell='4arcsec', stokes='I',
       deconvolver='multiscale', scales=[0, 2, 5, 10],
       smallscalebias=0.9,
       weighting='briggs', robust=0.5,
       pbcor=True, savemodel='modelcolumn')

# 3C286

tclean(vis='3C286.ms', imagename='3C286',
       mask='masks/3C286.mask',
       niter=100, gain=0.1, threshold='0.01Jy',
       imsize=512,
       cell='4arcsec', stokes='I',
       deconvolver='multiscale', scales=[0, 2, 5, 10],
       smallscalebias=0.9,
       weighting='briggs', robust=0.5,
       pbcor=True, savemodel='modelcolumn')

gaincal('3C286.ms', caltable='3C286.G0',
        solint='1200', refant=refant, minblperant=1,
        gaintype='G', calmode='ap', combine='scan', minsnr=1)

applycal('3C286.ms', gaintable='3C286.G0', calwt=False)

tclean(vis='3C286.ms', imagename='3C286_selfcal1',
       mask='masks/3C286.mask',
       niter=1000, gain=0.1, threshold='0.01Jy',
       imsize=512,
       cell='4arcsec', stokes='I',
       deconvolver='multiscale', scales=[0, 2, 5, 10],
       smallscalebias=0.9,
       weighting='briggs', robust=0.5,
       pbcor=True, savemodel='modelcolumn')

# J2007+404

tclean(vis='J2007+404.ms', imagename='J2007+404',
       mask='masks/J2007+404.mask',
       niter=100, gain=0.1, threshold='0.005Jy',
       imsize=512,
       cell='4arcsec', stokes='I',
       deconvolver='multiscale', scales=[0, 2, 5, 10],
       smallscalebias=0.9,
       weighting='briggs', robust=0.5,
       pbcor=True, savemodel='modelcolumn')

gaincal('J2007+404.ms', caltable='J2007+404.G0',
        solint='3600', refant=refant, minblperant=1,
        gaintype='G', calmode='ap', combine='scan', minsnr=1)

applycal('J2007+404.ms', gaintable='J2007+404.G0', calwt=False)

tclean(vis='J2007+404.ms', imagename='J2007+404_selfcal1',
       mask='masks/J2007+404.mask',
       niter=1000, gain=0.1, threshold='0.005Jy',
       imsize=512,
       cell='4arcsec', stokes='I',
       deconvolver='multiscale', scales=[0, 2, 5, 10],
       smallscalebias=0.9,
       weighting='briggs', robust=0.5,
       pbcor=True, savemodel='modelcolumn')

# CygA

tclean(vis='CygA.ms', imagename='CygA',
       mask='masks/CygA.mask',
       niter=1000, gain=0.1, threshold='0.1Jy',
       imsize=512,
       cell='4arcsec', stokes='I',
       deconvolver='multiscale', scales=[0, 2, 5, 10],
       smallscalebias=0.9,
       weighting='briggs', robust=0.5,
       pbcor=True, savemodel='modelcolumn')

gaincal(vis='CygA.ms', caltable='CygA.selfcal1',
        solint='600', refant=refant, combine='scan',
        gaintype='G', calmode='p', minblperant=1, minsnr=1)

applycal(vis='CygA.ms', gaintable=['CygA.selfcal1'],
         interp=['nearest'])

tclean(vis='CygA.ms', imagename='CygA_selfcal1',
       mask='masks/CygA.mask',
       niter=3000, gain=0.1, threshold='0.07Jy',
       imsize=512,
       cell='4arcsec', stokes='I',
       deconvolver='multiscale', scales=[0, 2, 5, 10],
       smallscalebias=0.9,
       weighting='briggs', robust=0.5,
       pbcor=True, savemodel='modelcolumn')

gaincal(vis='CygA.ms', caltable='CygA.selfcal2',
        solint='600', refant=refant, combine='scan',
        gaintype='G', calmode='ap', minblperant=1, minsnr=1)

applycal(vis='CygA.ms', gaintable=['CygA.selfcal2'],
         interp=['nearest'])

tclean(vis='CygA.ms', imagename='CygA_selfcal2',
       mask='masks/CygA.mask',
       niter=3000, gain=0.1, threshold='0.07Jy',
       imsize=512,
       cell='4arcsec', stokes='I',
       deconvolver='multiscale', scales=[0, 2, 5, 10],
       smallscalebias=0.9,
       weighting='briggs', robust=0.5,
       pbcor=True, savemodel='modelcolumn')

# Image stats

for im in [
        '3C84.image', '3C84_selfcal1.image',
        '3C286.image', '3C286_selfcal1.image',
        'J2007+404.image', 'J2007+404_selfcal1.image',
        'CygA.image', 'CygA_selfcal1.image', 'CygA_selfcal2.image']:
    print(im, f'SNR: {snr(im):.2f}')
