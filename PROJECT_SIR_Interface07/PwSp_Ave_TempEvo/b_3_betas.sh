cd b1.175_S2_Q0
../SIR_TE -P 500000 -T 2048 -D 0 -W 500 -K 2 -Z 60. -E 2. -d 5.5e-5 -e 1.e-5 -b 1.175 -B0 1.175 -B1 0.25 -m 0. -a 0. -g 0.077 -S 2 -Q 0 -Y 365. -R 400 > b1.175.out &
cd ..
cd b055_S2_Q0
../SIR_TE -P 500000 -T 2048 -D 0 -W 500 -K 2 -Z 60. -E 2. -d 5.5e-5 -e 1.e-5 -b 0.55 -B0 0.55 -B1 0.25 -m 0. -a 0. -g 0.077 -S 2 -Q 0 -Y 365. -R 400 > b055.out &
cd ..
cd b04_S2_Q0
../SIR_TE -P 500000 -T 2048 -D 0 -W 500 -K 2 -Z 60. -E 2. -d 5.5e-5 -e 1.e-5 -b 0.4 -B0 0.4 -B1 0.25 -m 0. -a 0. -g 0.077 -S 2 -Q 0 -Y 365. -R 400 > b04.out &
