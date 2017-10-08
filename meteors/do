ABmag_band_sed.perl -m 19 -b r > a
src_ray  -b 5330.09708737864 7107 2 -f 0.000277262362052392 -r 4500 -t 15 -s 0 0 -o 100000000 | tran_ray -r 2.5e-4 2.5e-4 0 >> a
src_ray  -b 5330.09708737864 7107 2 -f 0.000277262362052392 -r 4500 -t 15 -s 0 0 -o 50000000 | tran_ray -r 5e-4 5e-4 0 >> a
cat a | diffuser -k -s 0.6 | lsst_mirrors.perl -m 1008 -b r -doparts BTC | plane_collapse -c 0 0 1 0 > tmp
cat tmp | print_rays | niceplot.perl -c "finite sources" -c "1) D=50km" -c "2) D=100km" -c "3) D=infinite" > t.t
cat tmp | tran_ray -t -3 3 0 | image -l 8 -d 800 -e -o tmp.fits
