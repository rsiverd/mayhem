# Warnings emitted during a test of Z90/Z82/Z81 calibration building
# identified a number of FITS files in the Vandy /archive/engineering
# 'mirror' that were incomplete. This is likely the result of running
# out of disk space without blocking the auto-download script. This
# implies a bug in the fetching script which I still need to track
# down and fix. In the meantime, I am identifying affected files
# using fitsverify and removing them so that the fetcher can retrieve
# valid replacements. Relevant commands are reported below for posterity.

# -----------------------------------------------------------------------
# CORRUPT FILE DISCOVERY (krang):
krang
cdvae
cd lsc/nres01

# Bad fits.fz files from April/May:
ls 2019{04,05}??/raw/*.fz > z_fits_files.txt
wcl z_fits_files.txt
cat z_fits_files.txt | xargs -n10 fitsverify -e -q | tee z_fits_status.txt
grep -v 'OK:' z_fits_status.txt | cut -d' ' -f3 | cut -d, -f1 # inspect
grep -v 'OK:' z_fits_status.txt | cut -d' ' -f3 | cut -d, -f1 | xargs -n10 rm -v

# Bad fits.fz files from earlier (just in case):
ls 20190{1..3}??/raw/*.fz > z_fits_files.txt
wcl z_fits_files.txt
cat z_fits_files.txt | xargs -n10 fitsverify -e -q | tee z_fits_status.txt
grep -v 'OK:' z_fits_status.txt | cut -d' ' -f3 | cut -d, -f1 # inspect
grep -v 'OK:' z_fits_status.txt | cut -d' ' -f3 | cut -d, -f1 | xargs -n10 rm -v

rm z_fits_files.txt z_fits_status.txt

# -----------------------------------------------------------------------
# Repeat for other sites (using newly-written script):
krang
cdvae

find cpt/nres??/2019????/raw -type f -name "*.fits.fz" | sort > z_fits_files.lsc.txt
~/NRES/mayhem/bin/remove-invalid-fits-files.sh z_fits_files.lsc.txt
rm z_fits_files.lsc.txt

find elp/nres??/2019????/raw -type f -name "*.fits.fz" | sort > z_fits_files.elp.txt
~/NRES/mayhem/bin/remove-invalid-fits-files.sh z_fits_files.elp.txt
rm z_fits_files.elp.txt

find cpt/nres??/2019????/raw -type f -name "*.fits.fz" | sort > z_fits_files.cpt.txt
~/NRES/mayhem/bin/remove-invalid-fits-files.sh z_fits_files.cpt.txt
rm z_fits_files.cpt.txt

find tlv/nres??/2019????/raw -type f -name "*.fits.fz" | sort > z_fits_files.tlv.txt
~/NRES/mayhem/bin/remove-invalid-fits-files.sh z_fits_files.tlv.txt
rm z_fits_files.tlv.txt


rm z_fits_files.???.txt z_fits_files.txt

# -----------------------------------------------------------------------
# Check for bad tarballs:
krang
cdvae
cd lsc/nres01

du -k 2019????/specproc/*gz | awk '{print $2,$1}' > z_file_sizes.txt
awk '{ if ($2 < 900) print }' z_file_sizes.txt # inspect
awk '{ if ($2 < 900) print $1 }' z_file_sizes.txt | xargs -n10 rm -v

du -k ????????/specproc/*gz | awk '{print $2,$1}' > z_file_sizes.txt
awk '{ if ($2 < 900) print }' z_file_sizes.txt # inspect
awk '{ if ($2 < 900) print $1 }' z_file_sizes.txt | xargs -n10 rm -v

# Everything in one pass:
cdvae
du -k ???/nres??/*/specproc/*gz | awk '{print $2,$1}' > z_file_sizes.txt
awk '{ if ($2 < 900) print }' z_file_sizes.txt # inspect
awk '{ if ($2 < 900) print $1 }' z_file_sizes.txt | xargs -n10 rm -v

# -----------------------------------------------------------------------
# FILE REDOWNLOAD (splinter):
sssp
cd NRES/mayhem

# April/May:
./bin/fetch-nres-data.py -o /net/shredder.phy.vanderbilt.edu/hd1/siverd/lco_arch_eng --old-first -s lsc --start 2019-04-01 --end 2019-05-15

# earlier:
./bin/fetch-nres-data.py -o /net/shredder.phy.vanderbilt.edu/hd1/siverd/lco_arch_eng --old-first -s lsc --start 2019-01-01 --end 2019-04-01

./bin/fetch-nres-data.py -o /net/shredder.phy.vanderbilt.edu/hd1/siverd/lco_arch_eng --old-first -s elp --start 2019-04-10 --end 2019-04-20
./bin/fetch-nres-data.py -o /net/shredder.phy.vanderbilt.edu/hd1/siverd/lco_arch_eng --old-first -s cpt --start 2019-04-10 --end 2019-04-20
./bin/fetch-nres-data.py -o /net/shredder.phy.vanderbilt.edu/hd1/siverd/lco_arch_eng --old-first -s cpt --proc --start 2019-04-10 --end 2019-04-20

./bin/fetch-nres-data.py -o /net/shredder.phy.vanderbilt.edu/hd1/siverd/lco_arch_eng --old-first -s tlv --start 2019-04-10 --end 2019-04-20

