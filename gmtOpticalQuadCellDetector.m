function wfs = gmtOpticalQuadCellDetector(tel,loopRate,pixelScaleInArcsec,gs)

wfs = shackHartmann(1,2);
wfs.camera.exposureTime = 1/loopRate;
wfs.lenslets.throughput = 0.51;
wfs.camera.quantumEfficiency = 0.45;
wfs.camera.readOutNoise      = 0.1;
wfs.camera.pixelScale   = pixelScaleInArcsec;
d = tel.D/wfs.lenslets.nLenslet;
wfs.lenslets.fieldStopSize ...
                        = pixelScaleInArcsec*constants.arcsec2radian*d/gs.wavelength;