cd b02_S2_Q1
../SIR_TE -P 500000 -T 2048 -D 0 -W 500 -K 2 -Z 60. -E 2. -d 5.5e-5 -e 2.e-6 -b 1.50 -B0 1.50 -B1 0.25 -m 0. -a 0. -g 0.077 -S 2 -Q 1 -Y 365. -R 400 > b02.out &
cd ..
cd b02_S2_Q0
../SIR_TE -P 500000 -T 2048 -D 0 -W 500 -K 2 -Z 60. -E 2. -d 5.5e-5 -e 2.e-6 -b 1.50 -B0 1.50 -B1 0.25 -m 0. -a 0. -g 0.077 -S 2 -Q 0 -Y 365. -R 400 > b02.out &

