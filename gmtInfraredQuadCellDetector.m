function wfs = gmtInfraredQuadCellDetector(tel,loopRate,pixelScaleInArcsec,gs)

wfs = shackHartmann(1,2);
wfs.camera.exposureTime = 1/loopRate;
wfs.lenslets.throughput = 0.55;
wfs.camera.quantumEfficiency = 0.85;
wfs.camera.readOutNoise      = 5;
wfs.camera.pixelScale   = pixelScaleInArcsec;
d = tel.D/wfs.lenslets.nLenslet;
wfs.lenslets.fieldStopSize ...
                        = pixelScaleInArcsec*constants.arcsec2radian*d/gs.wavelength;